#GSASII - phase data display routines
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
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import numpy as np
import numpy.linalg as nl

VERY_LIGHT_GREY = wx.Colour(235,235,235)
WHITE = wx.Colour(255,255,255)

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

def UpdatePhaseData(self,item,data,oldPage):

    Atoms = []
    if self.dataDisplay:
        self.dataDisplay.Destroy()
    PhaseName = self.PatternTree.GetItemText(item)
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataFrame.SetLabel('Phase Data for '+PhaseName)
    self.dataFrame.CreateStatusBar()
    self.dataDisplay = G2gd.GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())

    def SetupGeneral():
        generalData = data['General']
        atomData = data['Atoms']
        generalData['AtomTypes'] = []
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['vdWRadii'] = []
        generalData['AtomMass'] = []
        generalData['Color'] = []
        generalData['Mydir'] = self.dirname
        colType = 1
        colSS = 7
        if generalData['Type'] =='macromolecular':
            colType = 4
            colSS = 10
        for atom in atomData:
            atom[colType] = atom[colType].lower().capitalize()              #force to standard form
            if generalData['AtomTypes'].count(atom[colType]):
                generalData['NoAtoms'][atom[colType]] += atom[colSS-1]*float(atom[colSS+1])
            elif atom[colType] != 'UNK':
                Info = G2elem.GetAtomInfo(atom[colType])
                generalData['AtomTypes'].append(atom[colType])
                generalData['BondRadii'].append(Info['Drad'])
                generalData['AngleRadii'].append(Info['Arad'])
                generalData['vdWRadii'].append(Info['Vdrad'])
                generalData['AtomMass'].append(Info['Mass'])
                generalData['NoAtoms'][atom[colType]] = atom[colSS-1]*float(atom[colSS+1])
                generalData['Color'].append(Info['Color'])

    def UpdateGeneral():
        
        ''' default dictionary structure for "General" phase item: (taken from GSASII.py)
        'General':{
            'Name':PhaseName
            'Type':'nuclear'
            'SGData':SGData
            'Cell':[False,10.,10.,10.,90.,90.,90,1000.]
            'Histogram list':['',]
            'Pawley dmin':1.0}
            'Atoms':[]
            'Drawing':{}
        })
        '''
        phaseTypes = ['nuclear','modulated','magnetic','macromolecular','Pawley']
        SetupGeneral()
        generalData = data['General']
        
        def OnPhaseName(event):
            oldName = generalData['Name']
            generalData['Name'] = NameTxt.GetValue()
            self.G2plotNB.Rename(oldName,generalData['Name'])
            self.dataFrame.SetLabel('Phase Data for '+generalData['Name'])
            self.PatternTree.SetItemText(item,generalData['Name'])
                        
        def OnPhaseType(event):
            if not generalData['AtomTypes']:             #can change only if no atoms!
                generalData['Type'] = TypeTxt.GetValue()
                dataDisplay.Destroy()           #needed to clear away bad cellSizer, etc.
                UpdateGeneral()
                if generalData['Type'] == 'Pawley':
                    if self.dataDisplay.FindPage('Atoms'):
                        self.dataDisplay.DeletePage(self.dataDisplay.FindPage('Atoms'))
                        self.dataDisplay.DeletePage(self.dataDisplay.FindPage('Draw Options'))
                        self.dataDisplay.DeletePage(self.dataDisplay.FindPage('Draw Atoms'))
                        self.dataDisplay.AdvanceSelection()
                    if not self.dataDisplay.FindPage('Pawley reflections'):      
                        self.dataDisplay.AddPage(G2gd.GSGrid(self.dataDisplay),'Pawley reflections')
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
                if ObjId == 0:
                    cell[1] = cell[2] = value
                    Obj.SetValue("%.5f"%(cell[1]))
                else:
                    cell[3] = value
                    Obj.SetValue("%.5f"%(cell[3]))
            elif laue in ['mmm']:
                cell[ObjId+1] = value
                Obj.SetValue("%.5f"%(cell[ObjId+1]))
            elif laue in ['2/m'+'a']:
                if ObjId != 3:
                    cell[ObjId+1] = value
                    Obj.SetValue("%.5f"%(cell[ObjId+1]))
                else:
                    cell[4] = value
                    Obj.SetValue("%.3f"%(cell[4]))
            elif laue in ['2/m'+'b']:
                if ObjId != 3:
                    cell[ObjId+1] = value
                    Obj.SetValue("%.5f"%(cell[ObjId+1]))
                else:
                    cell[5] = value
                    Obj.SetValue("%.3f"%(cell[5]))
            elif laue in ['2/m'+'c']:
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
            dataDisplay.Destroy()           #needed to clear away bad cellSizer, etc.
            UpdateGeneral()
                        
        def OnPawleyVal(event):
            try:
                dmin = float(pawlVal.GetValue())
                if 0.25 <= dmin <= 10.:
                    generalData['Pawley dmin'] = dmin
            except ValueError:
                pass
            pawlVal.SetValue("%.2f"%(generalData['Pawley dmin']))          #reset in case of error            
                                    
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
        
        General.DestroyChildren()
        dataDisplay = wx.Panel(General)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'General phase data:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        nameSizer = wx.BoxSizer(wx.HORIZONTAL)
        nameSizer.Add(wx.StaticText(dataDisplay,-1,' Phase name: '),0,wx.ALIGN_CENTER_VERTICAL)
        NameTxt = wx.TextCtrl(dataDisplay,-1,value=generalData['Name'],style=wx.TE_PROCESS_ENTER)
        NameTxt.Bind(wx.EVT_TEXT_ENTER,OnPhaseName)
        NameTxt.Bind(wx.EVT_KILL_FOCUS,OnPhaseName)
        nameSizer.Add(NameTxt,0,wx.ALIGN_CENTER_VERTICAL)
        nameSizer.Add(wx.StaticText(dataDisplay,-1,'  Phase type: '),0,wx.ALIGN_CENTER_VERTICAL)
        TypeTxt = wx.ComboBox(dataDisplay,-1,value=generalData['Type'],choices=phaseTypes,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        TypeTxt.Bind(wx.EVT_COMBOBOX, OnPhaseType)
        nameSizer.Add(TypeTxt,0,wx.ALIGN_CENTER_VERTICAL)
        nameSizer.Add(wx.StaticText(dataDisplay,-1,'  Space group: '),0,wx.ALIGN_CENTER_VERTICAL)
        SGTxt = wx.TextCtrl(dataDisplay,-1,value=generalData['SGData']['SpGrp'],style=wx.TE_PROCESS_ENTER)
        SGTxt.Bind(wx.EVT_TEXT_ENTER,OnSpaceGroup)
        nameSizer.Add(SGTxt,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(nameSizer,0)
        mainSizer.Add((5,5),0)
        cell = generalData['Cell']
        laue = generalData['SGData']['SGLaue']
        if laue == '2/m':
            laue += generalData['SGData']['SGUniq']
        for cellGUI in cellGUIlist:
            if laue in cellGUI[0]:
                useGUI = cellGUI
        cellList = []
        cellSizer = wx.FlexGridSizer(2,useGUI[1]+1,5,5)
        cellRef = wx.CheckBox(dataDisplay,label='Refine unit cell:')
        cellSizer.Add(cellRef,0,wx.ALIGN_CENTER_VERTICAL)
        cellRef.Bind(wx.EVT_CHECKBOX, OnCellRef)
        cellRef.SetValue(cell[0])
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
        mainSizer.Add(cellSizer,0)
        mainSizer.Add((5,5),0)
        
        if len(generalData['AtomTypes']):
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
            mainSizer.Add(denSizer)
            mainSizer.Add((5,5),0)
            
            elemSizer = wx.FlexGridSizer(7,len(generalData['AtomTypes'])+1,1,1)
            elemSizer.Add(wx.StaticText(dataDisplay,label='Elements'),0,wx.ALIGN_CENTER_VERTICAL)
            for elem in generalData['AtomTypes']:
                typTxt = wx.TextCtrl(dataDisplay,value=elem,style=wx.TE_READONLY)
                typTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(typTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='No. per cell'),0,wx.ALIGN_CENTER_VERTICAL)
            for elem in generalData['AtomTypes']:
                numbTxt = wx.TextCtrl(dataDisplay,value='%.1f'%(generalData['NoAtoms'][elem]),
                    style=wx.TE_READONLY)
                numbTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(numbTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='Atom weight'),0,wx.ALIGN_CENTER_VERTICAL)
            for wt in generalData['AtomMass']:
                wtTxt = wx.TextCtrl(dataDisplay,value='%.3f'%(wt),style=wx.TE_READONLY)
                wtTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(wtTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='Bond radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['BondRadii']:
                bondRadii = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                bondRadii.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(bondRadii,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='Angle radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['AngleRadii']:
                elemTxt = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='van der Waals radii'),0,wx.ALIGN_CENTER_VERTICAL)
            for rad in generalData['vdWRadii']:
                elemTxt = wx.TextCtrl(dataDisplay,value='%.2f'%(rad),style=wx.TE_READONLY)
                elemTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                elemSizer.Add(elemTxt,0,wx.ALIGN_CENTER_VERTICAL)
            elemSizer.Add(wx.StaticText(dataDisplay,label='Default color'),0,wx.ALIGN_CENTER_VERTICAL)
            for color in generalData['Color']:
                colorTxt = wx.TextCtrl(dataDisplay,value='',style=wx.TE_READONLY)
                colorTxt.SetBackgroundColour(color)
                elemSizer.Add(colorTxt,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add(elemSizer)
            
        elif generalData['Type'] == 'Pawley':
            pawlSizer = wx.BoxSizer(wx.HORIZONTAL)
            pawlSizer.Add(wx.StaticText(dataDisplay,label=' Pawley dmin: '),0,wx.ALIGN_CENTER_VERTICAL)
            pawlVal = wx.TextCtrl(dataDisplay,value='%.2f'%(generalData['Pawley dmin']),style=wx.TE_PROCESS_ENTER)
            pawlVal.Bind(wx.EVT_TEXT_ENTER,OnPawleyVal)        
            pawlVal.Bind(wx.EVT_KILL_FOCUS,OnPawleyVal)
            pawlSizer.Add(pawlVal,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add(pawlSizer)

        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)

    def FillAtomsGrid():

        self.dataFrame.setSizePosLeft([700,300])
        generalData = data['General']
        atomData = data['Atoms']
        Items = [G2gd.wxID_ATOMSEDITINSERT, G2gd.wxID_ATOMSEDITDELETE, G2gd.wxID_ATOMSREFINE, 
            G2gd.wxID_ATOMSMODIFY, G2gd.wxID_ATOMSTRANSFORM, G2gd.wxID_ATONTESTINSERT]
        if atomData:
            for item in Items:    
                self.dataFrame.AtomsMenu.Enable(item,True)
        else:
            for item in Items:
                self.dataFrame.AtomsMenu.Enable(item,False)            
            
        AAchoice = ": ,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE,HOH,UNK"
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,X,XU,U,F,FX,FXU,FU",
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,4', #x,y,z,frac
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+":I,A",
            wg.GRID_VALUE_FLOAT+':10,4',                                                            #Uiso
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,                         #Uij - placeholders
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING]
        colLabels = ['Name','Type','refine','x','y','z','frac','site sym','mult','I/A','Uiso','U11','U22','U33','U12','U13','U23']
        if generalData['Type'] == 'magnetic':
            colLabels += ['Mx','My','Mz']
            Types[2] = wg.GRID_VALUE_CHOICE+": ,X,XU,U,M,MX,MXU,MU,F,FX,FXU,FU,FM,FMX,FMU,"
            Types += [
                wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
        elif generalData['Type'] == 'macromolecular':
            colLabels = ['res no','residue','chain'] + colLabels
            Types = [wg.GRID_VALUE_STRING,
                wg.GRID_VALUE_CHOICE+AAchoice,
                wg.GRID_VALUE_STRING] + Types
        elif generalData['Type'] == 'modulated':
            Types += []
            colLabels += []

        def RefreshAtomGrid(event):

            def chkUij(Uij,CSI):
                return Uij

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
                    dlg = wx.MultiChoiceDialog(self,'Select','Refinement controls',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelections()
                        parms = ''
                        for x in sel:
                            parms += choice[x][0]
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'I/A':
                    choice = ['Isotropic','Anisotropic']
                    dlg = wx.SingleChoiceDialog(self,'Select','Thermal Motion',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel][0]
                    dlg.Destroy()
                elif Atoms.GetColLabelValue(c) == 'Type':
                    choice = generalData['AtomTypes']
                    dlg = wx.SingleChoiceDialog(self,'Select','Atom types',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel]
                        noSkip = False
                        Atoms.ClearSelection()
                        for row in range(Atoms.GetNumberRows()):
                            if parms == atomData[row][c]:
                                Atoms.SelectRow(row,True)
                elif Atoms.GetColLabelValue(c) == 'residue':
                    choice = []
                    for r in range(Atoms.GetNumberRows()):
                        if str(atomData[r][c]) not in choice:
                            choice.append(str(atomData[r][c]))
                    choice.sort()
                    dlg = wx.SingleChoiceDialog(self,'Select','Residue',choice)
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
                    dlg = wx.SingleChoiceDialog(self,'Select','Residue no.',choice)
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
                    dlg = wx.SingleChoiceDialog(self,'Select','Chain',choice)
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
                        if parms != atomData[r][c] and Atoms.GetColLabelValue(c) == 'I/A':
                            if parms == 'A':                #'I' --> 'A'
                                Uiso = atomData[r][us]
                                sytsym = atomData[r][ss]
                                CSI = G2spc.GetCSuinel(sytsym)
                                atomData[r][ui:ui+6] = Uiso*np.array(CSI[3])
                                atomData[r][us] = ''
                                Atoms.SetCellRenderer(r,us,wg.GridCellStringRenderer())
                                Atoms.SetCellStyle(r,us,VERY_LIGHT_GREY,True)
                                for i in range(6):
                                    ci = ui+i
                                    Atoms.SetCellRenderer(r,ci,wg.GridCellFloatRenderer(10,4))
                                    Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                                    if CSI[2][i]:
                                        Atoms.SetCellStyle(r,ci,WHITE,False)
                            else:                           #'A' --> 'I'
                                Uij = atomData[r][ui:ui+6]
                                atomData[r][us] = (Uij[0]+Uij[1]+Uij[2])/3.0
                                atomData[r][ui:ui+6] = [0,0,0,0,0,0]
                                Atoms.SetCellRenderer(r,us,wg.GridCellFloatRenderer(10,4))
                                Atoms.SetCellStyle(r,us,WHITE,False)
                                for i in range(6):
                                    ci = ui+i
                                    Atoms.SetCellRenderer(r,ci,wg.GridCellStringRenderer())
                                    Atoms.SetCellValue(r,ci,'')
                                    Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                        Atoms.SetCellValue(r,c,parms)
            elif Atoms.GetColLabelValue(c) in ['Name']:
                value = Atoms.GetCellValue(r,c)
                atomData[r][c] = value
                ID = [atomData[r][-1],]
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByIDs(data['Drawing'],atomData[r],ID)
            elif Atoms.GetColLabelValue(c) in ['x','y','z']:
                atomData[r][c] = float(Atoms.GetCellValue(r,c))
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
                Atoms.SetCellValue(r,SScol,Sytsym)
                Atoms.SetCellValue(r,Mulcol,str(Mult))
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
            elif Atoms.GetColLabelValue(c) == 'I/A':
                atomData[r][c] = Atoms.GetCellValue(r,c)
                if atomData[r][c] == 'I':
                    Uij = atomData[r][c+2:c+8]
                    atomData[r][c+1] = (Uij[0]+Uij[1]+Uij[2])/3.0
                    atomData[r][c+2:c+8] = [0,0,0,0,0,0]
                    Atoms.SetCellRenderer(r,c+1,wg.GridCellFloatRenderer(10,4))
                    Atoms.SetCellStyle(r,c+1,WHITE,False)
                    for i in range(6):
                        ci = i+colLabels.index('U11')
                        Atoms.SetCellRenderer(r,ci,wg.GridCellStringRenderer())
                        Atoms.SetCellValue(r,ci,'')
                        Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                else:
                    Uiso = atomData[r][c+1]
                    atomData[r][c+1] = ''
                    CSI = G2spc.GetCSuinel(atomData[r][colLabels.index('site sym')])
                    atomData[r][c+2:c+8] = Uiso*np.array(CSI[3])
                    Atoms.SetCellRenderer(r,c+1,wg.GridCellStringRenderer())
                    Atoms.SetCellStyle(r,c+1,VERY_LIGHT_GREY,True)
                    for i in range(6):
                        ci = i+colLabels.index('U11')
                        Atoms.SetCellRenderer(r,ci,wg.GridCellFloatRenderer(10,4))
                        Atoms.SetCellStyle(r,ci,VERY_LIGHT_GREY,True)
                        if CSI[2][i]:
                            Atoms.SetCellStyle(r,ci,WHITE,False)
            elif Atoms.GetColLabelValue(c) in ['U11','U22','U33','U12','U13','U23']:
                atomData[r][c] = float(Atoms.GetCellValue(r,c))
                CSI = G2spc.GetCSuinel(atomData[r][colLabels.index('site sym')])
                iUij = CSI[0][c-colLabels.index('U11')]
                for i in range(6):
                    if iUij == CSI[0][i]:
                        atomData[r][i+colLabels.index('U11')] = value*CSI[1][i]
            elif Atoms.GetColLabelValue(c) == 'Uiso':
                atomData[r][c] = float(Atoms.GetCellValue(r,c))                
            data['Atoms'] = atomData

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

        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if Atoms.GetColLabelValue(c) == 'Type':
                PE = G2elem.PickElement(self)
                if PE.ShowModal() == wx.ID_OK:
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
                ID = [atomData[r][-1],]
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByIDs(data['Drawing'],atomData[r],ID)
            else:
                event.Skip()

        table = []
        rowLabels = []
        for i,atom in enumerate(atomData):
            if atom[colLabels.index('I/A')] == 'I':
                table.append(atom[:colLabels.index('U11')]+['','','','','',''])
            else:
                table.append(atom)
            rowLabels.append(str(i+1))
        atomTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        Atoms.SetTable(atomTable, True)
        Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, ChangeSelection)
        Atoms.Bind(wg.EVT_GRID_SELECT_CELL, AtomTypeSelect)
        Atoms.SetMargins(0,0)
        Atoms.AutoSizeColumns(True)
        colType = colLabels.index('Type')
        colSS = colLabels.index('site sym')
        colIA = colLabels.index('I/A')
        colU11 = colLabels.index('U11')
        colUiso = colLabels.index('Uiso')
        attr = wg.GridCellAttr()                                #to set Uij defaults
        attr.SetBackgroundColour(VERY_LIGHT_GREY)
        attr.SetReadOnly(True)
        for i in range(colU11,colU11+6):
            Atoms.SetColAttr(i,attr)
        for row in range(Atoms.GetNumberRows()):
            Atoms.SetReadOnly(row,colSS,True)                         #site sym
            Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            if Atoms.GetCellValue(row,colIA) == 'A':
                CSI = G2spc.GetCSuinel(atomData[row][colLabels.index('site sym')])
                Atoms.SetCellRenderer(row,colUiso,wg.GridCellStringRenderer())
                Atoms.SetCellStyle(row,colUiso,VERY_LIGHT_GREY,True)
                Atoms.SetCellValue(row,colUiso,'')
                for i in range(6):
                    ci = colU11+i
                    Atoms.SetCellRenderer(row,ci,wg.GridCellFloatRenderer(10,4))
                    if CSI[2][i]:
                        Atoms.SetCellStyle(row,ci,WHITE,False)

    def OnAtomAdd(event):
        AtomAdd(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
        
    def OnAtomTestAdd(event):
        try:
            drawData = data['Drawing']
            x,y,z = drawData['testPos']
            AtomAdd(x,y,z)
        except:
            AtomAdd(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
                
    def AtomAdd(x,y,z):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        atId = ran.randint(0,sys.maxint)
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','H','',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'nuclear':
            atomData.append(['UNK','H','',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'magnetic':
            atomData.append(['UNK','H','',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0,atId])
        SetupGeneral()
        if 'Drawing' in data:            
            DrawAtomAdd(data['Drawing'],atomData[-1])
            G2plt.PlotStructure(self,data)

    def OnAtomInsert(event):
        AtomInsert(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
        
    def OnAtomTestInsert(event):
        if 'Drawing' in data:
            drawData = data['Drawing']
            x,y,z = drawData['testPos']
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
                G2plt.PlotStructure(self,data)
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
            dlg = wx.MultiChoiceDialog(self,'Select','Refinement controls',choice)
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

    def AtomTransform(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = colLabels.index('U11')
            css = colLabels.index('site sym')
            atomData = data['Atoms']
            generalData = data['General']
            SGData = generalData['SGData']
            dlg = SymOpDialog(self,SGData,True)
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

    def SetupDrawingData():
        generalData = data['General']
        atomData = data['Atoms']
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        defaultDrawing = {'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,'backColor':[0,0,0],'depthFog':False,
            'Zclip':50.0,'cameraPos':50.,'radiusFactor':0.85,'showBadContacts':False,
            'bondRadius':0.1,'ballScale':0.33,'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,
            'unitCellBox':False,'showABC':True,'showSymElem':False,'selectedAtoms':[],
            'Rotation':[0.0,0.0,0.0,[]],'bondList':{},'testPos':[-.1,-.1,-.1]}
        try:
            drawingData = data['Drawing']
        except KeyError:
            data['Drawing'] = {}
            drawingData = data['Drawing']
        if not drawingData:                 #fill with defaults if empty
            drawingData = copy.copy(defaultDrawing)
            drawingData['Atoms'] = []
        cx,ct,cs = [0,0,0]
        if generalData['Type'] == 'nuclear':
            cx,ct,cs = [2,1,6]         #x, type & style
        elif generalData['Type'] == 'macromolecular':
            cx,ct,cs = [5,4,9]         #x, type & style
        elif generalData['Type'] == 'magnetic':
            cx,ct,cs = [2,1,6]         #x, type & style
#        elif generalData['Type'] == 'modulated':
#           ?????   for future
        if not drawingData['Atoms']:
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
        if generalData['Type'] == 'nuclear':
            if oldatom:
                atomInfo = [atom[:2]+oldatom[2:]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+['1',]+['vdW balls',]+
                    ['',]+[[255,255,255],]+atom[9:]+[[]]][0]
            ct,cs = [1,8]         #type & color
        elif generalData['Type'] == 'macromolecular':
            try:
                oneLetter = AA3letter.index(atom[1])
            except ValueError:
                oneLetter = -1
            atomInfo = [[atom[1].strip()+atom[0],]+
                [AA1letter[oneLetter]+atom[0],]+atom[2:5]+
                atom[6:9]+['1',]+['sticks',]+['',]+[[255,255,255],]+atom[12:]+[[]]][0]
            ct,cs = [4,11]         #type & color
        elif generalData['Type'] == 'magnetic':
            if oldatom:
                atomInfo = [atom[:2]+oldatom[3:]][0]
            else:
                atomInfo = [atom[:2]+atom[3:6]+['vdW balls',]+['',]+atom[9:]+[[]]][0]
            ct,cs = [1,8]         #type & color
#        elif generalData['Type'] == 'modulated':
#           ?????   for future
        atNum = generalData['AtomTypes'].index(atom[ct])
        atomInfo[cs] = list(generalData['Color'][atNum])
        return atomInfo
            
    def DrawAtomAdd(drawingData,atom):
        drawingData['Atoms'].append(MakeDrawAtom(atom))
        
    def DrawAtomsReplaceByIDs(drawingData,atom,IDs):
        atomData = drawingData['Atoms']
        indx = FindAtomIndexByIDs(atomData,IDs)
        for ind in indx:
            atomData[ind] = MakeDrawAtom(atom,atomData[ind])
            
    def UpdateDrawAtoms():
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
        cx,ct,cs = drawingData['atomPtrs']
        atomData = drawingData['Atoms']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',    #x,y,z
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,lines,vdW balls,sticks,balls & sticks,ellipsoids,polyhedra",
            wg.GRID_VALUE_CHOICE+": ,type,name,number",wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,]
        styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']
        labelChoice = [' ','type','name','number']
        colLabels = ['Name','Type','x','y','z','Sym Op','Style','Label','Color','I/A']
        if generalData['Type'] == 'macromolecular':
            colLabels = ['Residue','1-letter','Chain'] + colLabels
            Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING]+Types
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

                dlg = wx.MultiChoiceDialog(self,'Select',name,choice)
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
                    G2plt.PlotStructure(self,data)                    
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
                    dlg = wx.SingleChoiceDialog(self,'Select','Atom drawing style',styleChoice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = styleChoice[sel]
                        for r in range(len(atomData)):
                            atomData[r][c] = parms
                            drawAtoms.SetCellValue(r,c,parms)
                        FindBonds()
                        G2plt.PlotStructure(self,data)
                    dlg.Destroy()
                elif drawAtoms.GetColLabelValue(c) == 'Label':
                    dlg = wx.SingleChoiceDialog(self,'Select','Atom labelling style',labelChoice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = labelChoice[sel]
                        for r in range(len(atomData)):
                            atomData[r][c] = parms
                            drawAtoms.SetCellValue(r,c,parms)
                    dlg.Destroy()                    
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    dlg = wx.ColourDialog(self)
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
                    FindBonds()
                elif drawAtoms.GetColLabelValue(c) == 'Color':
                    dlg = wx.ColourDialog(self)
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
            G2plt.PlotStructure(self,data)
                    
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
            G2plt.PlotStructure(self,data)                    
                
        table = []
        rowLabels = []
        for i,atom in enumerate(drawingData['Atoms']):
            table.append(atom[:colLabels.index('I/A')+1])
            rowLabels.append(str(i+1))

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
        self.dataFrame.setSizePosLeft([600,300])
        
        FindBonds()
        G2plt.PlotStructure(self,data)

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
            dlg = wx.SingleChoiceDialog(self,'Select','Atom drawing style',styleChoice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parms = styleChoice[sel]
                for r in indx:
                    atomData[r][cs] = parms
                    drawAtoms.SetCellValue(r,cs,parms)
            dlg.Destroy()
            FindBonds()
            G2plt.PlotStructure(self,data)

    def DrawAtomLabel(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs = data['Drawing']['atomPtrs']
            styleChoice = [' ','type','name','number']
            if generalData['Type'] == 'macromolecular':
                styleChoice = [' ','type','name','number','residue','1-letter','chain']
            dlg = wx.SingleChoiceDialog(self,'Select','Atom label style',styleChoice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parms = styleChoice[sel]
                for r in indx:
                    atomData[r][cs+1] = parms
                    drawAtoms.SetCellValue(r,cs+1,parms)
            dlg.Destroy()
            G2plt.PlotStructure(self,data)
            
    def DrawAtomColor(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            generalData = data['General']
            atomData = data['Drawing']['Atoms']
            cx,ct,cs = data['Drawing']['atomPtrs']
            dlg = wx.ColourDialog(self)
            if dlg.ShowModal() == wx.ID_OK:
                color = dlg.GetColourData().GetColour()
                for r in indx:
                    atomData[r][cs+2] = color
                    attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                    attr.SetBackgroundColour(color)
                    drawAtoms.SetAttr(r,cs+2,attr)
                    data['Drawing']['Atoms'][r][cs+2] = color
            dlg.Destroy()
            
    def ResetAtomColors(event):
        generalData = data['General']
        atomData = data['Drawing']['Atoms']
        cx,ct,cs = data['Drawing']['atomPtrs']
        for atom in atomData:            
            atNum = generalData['AtomTypes'].index(atom[ct])
            atom[cs+2] = list(generalData['Color'][atNum])
        UpdateDrawAtoms()
        G2plt.PlotStructure(self,data)        
        
    def SetViewPoint(event):
        indx = drawAtoms.GetSelectedRows()
        if indx:
            atomData = data['Drawing']['Atoms']
            cx = data['Drawing']['atomPtrs'][0]
            data['Drawing']['viewPoint'] = [atomData[indx[0]][cx:cx+3],[indx[0],0]]
            drawAtoms.ClearSelection()                                  #do I really want to do this?
            G2plt.PlotStructure(self,data)
            
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
            dlg = SymOpDialog(self,SGData,False)
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
                            atom[cx+3] = str(((Opr+1)+100*Cent)*(1-2*Inv))+'+'+ \
                                str(int(Cell[0]))+str(int(Cell[1]))+str(int(Cell[2]))
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
            G2plt.PlotStructure(self,data)
            
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
            dlg = SymOpDialog(self,SGData,False)
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
                        atom[cx+3] = str(((Opr+1)+100*Cent)*(1-2*Inv))+'+'+ \
                            str(int(Cell[0]))+str(int(Cell[1]))+str(int(Cell[2]))
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
            G2plt.PlotStructure(self,data)
            
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
                                newAtom = atomB[:]
                                newAtom[cx:cx+3] = xyz
                                atomData.append(newAtom)
            data['Drawing']['Atoms'] = atomData
            UpdateDrawAtoms()
            G2plt.PlotStructure(self,data)
            
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
                    result = G2spc.GenAtom(XYZ,SGData,False,Uij)
                    for item in result:
                        atom = copy.copy(atomData[ind])
                        atom[cx:cx+3] = item[0]
                        atom[cx+3] = str(item[2])+'+000'
                        atom[cuij:cuij+6] = item[1]
                        Opp = G2spc.Opposite(item[0])
                        for xyz in Opp:
                            if noDuplicate(xyz,atomData):
                                atom[cx:cx+3] = xyz
                                atomData.append(atom[:])
                else:
                    result = G2spc.GenAtom(XYZ,SGData,False)
                    for item in result:
                        atom = copy.copy(atomData[ind])
                        atom[cx:cx+3] = item[0]
                        atom[cx+3] = str(item[1])+'+000'
                        Opp = G2spc.Opposite(item[0])
                        for xyz in Opp:
                            if noDuplicate(xyz,atomData):
                                atom[cx:cx+3] = xyz
                                atomData.append(atom[:])               
                data['Drawing']['Atoms'] = atomData
                
            UpdateDrawAtoms()
            G2plt.PlotStructure(self,data)
            
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
                    
    def FindBonds():                    #uses numpy & masks - very fast even for proteins!
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
            atom[-1] = []               #clear out old bonds
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
                dist = ma.masked_less(np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0)),0.5) #gets rid of self & disorder "bonds" < 0.5A
                sumR = atomA[3]+Radii
                IndB = ma.nonzero(ma.masked_greater(dist-data['Drawing']['radiusFactor']*sumR,0.))                 #get indices of bonded atoms
                i = atomA[0]
                for j in IndB[0]:
                    if Styles[i] == 'polyhedra':
                        atomData[i][-1].append(np.inner(Amat,Dx[j]))
                    elif Styles[j] != 'polyhedra' and j > i:
                        atomData[i][-1].append(Dx[j]*Radii[i]/sumR[j])
                        atomData[j][-1].append(-Dx[j]*Radii[j]/sumR[j])

    def DrawAtomsDelete(event):   
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
            atomData = data['Drawing']['Atoms']
            indx.reverse()
            for ind in indx:
                del atomData[ind]
            UpdateDrawAtoms()
            G2plt.PlotStructure(self,data)
        event.StopPropagation()
        
    def FindAtomIndexByIDs(atomData,IDs):
        indx = []
        for i,atom in enumerate(atomData):
            if atom[-2] in IDs:
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

        def OnZclip(event):
            drawingData['Zclip'] = Zclip.GetValue()
            ZclipTxt.SetLabel('Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
            G2plt.PlotStructure(self,data)
            
        def OnCameraPos(event):
            drawingData['cameraPos'] = cameraPos.GetValue()
            cameraPosTxt.SetLabel('Camera Position: '+'%.2f'%(drawingData['cameraPos']))
            ZclipTxt.SetLabel('Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
            G2plt.PlotStructure(self,data)

        def OnBackColor(event):
            drawingData['backColor'] = event.GetValue()
            G2plt.PlotStructure(self,data)

        def OnBallScale(event):
            drawingData['ballScale'] = ballScale.GetValue()/100.
            ballScaleTxt.SetLabel('Ball scale: '+'%.2f'%(drawingData['ballScale']))
            G2plt.PlotStructure(self,data)

        def OnVdWScale(event):
            drawingData['vdwScale'] = vdwScale.GetValue()/100.
            vdwScaleTxt.SetLabel('van der Waals scale: '+'%.2f'%(drawingData['vdwScale']))
            G2plt.PlotStructure(self,data)

        def OnEllipseProb(event):
            drawingData['ellipseProb'] = ellipseProb.GetValue()
            ellipseProbTxt.SetLabel('Ellipsoid probability: '+'%d%%'%(drawingData['ellipseProb']))
            G2plt.PlotStructure(self,data)

        def OnBondRadius(event):
            drawingData['bondRadius'] = bondRadius.GetValue()/100.
            bondRadiusTxt.SetLabel('Bond radius, A: '+'%.2f'%(drawingData['bondRadius']))
            G2plt.PlotStructure(self,data)

        def OnShowABC(event):
            drawingData['showABC'] = showABC.GetValue()
            G2plt.PlotStructure(self,data)

        def OnShowUnitCell(event):
            drawingData['unitCellBox'] = unitCellBox.GetValue()
            G2plt.PlotStructure(self,data)

        def OnShowBadContacts(event):
            drawingData['showBadContacts'] = showBadContacts.GetValue()

        def OnShowSymElem(event):
            drawingData['showSymElem'] = showSymElem.GetValue()

        def OnShowHyd(event):
            drawingData['showHydrogen'] = showHydrogen.GetValue()
            FindBonds()
            G2plt.PlotStructure(self,data)

        def OnSizeHatoms(event):
            try:
                value = max(0.1,min(1.2,float(sizeH.GetValue())))
            except ValueError:
                value = 0.5
            drawingData['sizeH'] = value
            sizeH.SetValue("%.2f"%(value))
            G2plt.PlotStructure(self,data)
            
        def OnRadFactor(event):
            try:
                value = max(0.1,min(1.2,float(radFactor.GetValue())))
            except ValueError:
                value = 0.85
            drawingData['radiusFactor'] = value
            radFactor.SetValue("%.2f"%(value))
            FindBonds()
            G2plt.PlotStructure(self,data)
            
        


        dataDisplay = wx.Panel(drawOptions)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Drawing controls:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        
        slopSizer = wx.BoxSizer(wx.HORIZONTAL)
        slideSizer = wx.FlexGridSizer(6,2)
        slideSizer.AddGrowableCol(1,1)

        cameraPosTxt = wx.StaticText(dataDisplay,-1,
            'Camera Position: '+'%.2f'%(drawingData['cameraPos']),name='cameraPos')
        slideSizer.Add(cameraPosTxt,0,wx.ALIGN_CENTER_VERTICAL)
        cameraPos = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['cameraPos'],name='cameraSlider')
        cameraPos.SetRange(10,500)
        cameraPos.Bind(wx.EVT_SLIDER, OnCameraPos)
        slideSizer.Add(cameraPos,1,wx.EXPAND|wx.RIGHT)
        
        ZclipTxt = wx.StaticText(dataDisplay,-1,'Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
        slideSizer.Add(ZclipTxt,0,wx.ALIGN_CENTER_VERTICAL)
        Zclip = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['Zclip'])
        Zclip.SetRange(1,99)
        Zclip.Bind(wx.EVT_SLIDER, OnZclip)
        slideSizer.Add(Zclip,1,wx.EXPAND|wx.RIGHT)
        
        vdwScaleTxt = wx.StaticText(dataDisplay,-1,'van der Waals scale: '+'%.2f'%(drawingData['vdwScale']))
        slideSizer.Add(vdwScaleTxt,0,wx.ALIGN_CENTER_VERTICAL)
        vdwScale = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['vdwScale']))
        vdwScale.Bind(wx.EVT_SLIDER, OnVdWScale)
        slideSizer.Add(vdwScale,1,wx.EXPAND|wx.RIGHT)

        ellipseProbTxt = wx.StaticText(dataDisplay,-1,'Ellipsoid probability: '+'%d%%'%(drawingData['ellipseProb']))
        slideSizer.Add(ellipseProbTxt,0,wx.ALIGN_CENTER_VERTICAL)
        ellipseProb = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=drawingData['ellipseProb'])
        ellipseProb.SetRange(1,99)
        ellipseProb.Bind(wx.EVT_SLIDER, OnEllipseProb)
        slideSizer.Add(ellipseProb,1,wx.EXPAND|wx.RIGHT)

        ballScaleTxt = wx.StaticText(dataDisplay,-1,'Ball scale: '+'%.2f'%(drawingData['ballScale']))
        slideSizer.Add(ballScaleTxt,0,wx.ALIGN_CENTER_VERTICAL)
        ballScale = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['ballScale']))
        ballScale.Bind(wx.EVT_SLIDER, OnBallScale)
        slideSizer.Add(ballScale,1,wx.EXPAND|wx.RIGHT)

        bondRadiusTxt = wx.StaticText(dataDisplay,-1,'Bond radius, A: '+'%.2f'%(drawingData['bondRadius']))
        slideSizer.Add(bondRadiusTxt,0,wx.ALIGN_CENTER_VERTICAL)
        bondRadius = wx.Slider(dataDisplay,style=wx.SL_HORIZONTAL,value=int(100*drawingData['bondRadius']))
        bondRadius.SetRange(1,25)
        bondRadius.Bind(wx.EVT_SLIDER, OnBondRadius)
        slideSizer.Add(bondRadius,1,wx.EXPAND|wx.RIGHT)
        
        slopSizer.Add(slideSizer,1,wx.EXPAND|wx.RIGHT)
        slopSizer.Add((10,5),0)
        slopSizer.SetMinSize(wx.Size(300,10))
        mainSizer.Add(slopSizer,0)
        mainSizer.Add((5,5),0)

        flexSizer = wx.FlexGridSizer(6,2,5,0)
        flexSizer.Add(wx.StaticText(dataDisplay,-1,'View Point:  '),0,wx.ALIGN_CENTER_VERTICAL)
        VP = drawingData['viewPoint'][0]
        viewPoint = wx.TextCtrl(dataDisplay,value='%.3f, %.3f, %.3f'%(VP[0],VP[1],VP[2]),
            style=wx.TE_READONLY,size=wx.Size(120,20),name='viewPoint')
        viewPoint.SetBackgroundColour(VERY_LIGHT_GREY)
        flexSizer.Add(viewPoint,0,wx.ALIGN_CENTER_VERTICAL)
        
        showABC = wx.CheckBox(dataDisplay,-1,label='Show test point?')
        showABC.Bind(wx.EVT_CHECKBOX, OnShowABC)
        showABC.SetValue(drawingData['showABC'])
        flexSizer.Add(showABC,0,wx.ALIGN_CENTER_VERTICAL)

        unitCellBox = wx.CheckBox(dataDisplay,-1,label='Show unit cell?')
        unitCellBox.Bind(wx.EVT_CHECKBOX, OnShowUnitCell)
        unitCellBox.SetValue(drawingData['unitCellBox'])
        flexSizer.Add(unitCellBox,0,wx.ALIGN_CENTER_VERTICAL)

        showBadContacts = wx.CheckBox(dataDisplay,-1,label='Show bad contacts?')
        showBadContacts.Bind(wx.EVT_CHECKBOX, OnShowBadContacts)
        showBadContacts.SetValue(drawingData['showBadContacts'])
        flexSizer.Add(showBadContacts,0,wx.ALIGN_CENTER_VERTICAL)

        showSymElem = wx.CheckBox(dataDisplay,-1,label='Show sym. elem.?')
        showSymElem.Bind(wx.EVT_CHECKBOX, OnShowSymElem)
        showSymElem.SetValue(drawingData['showSymElem'])
        flexSizer.Add(showSymElem,0,wx.ALIGN_CENTER_VERTICAL)

        showHydrogen = wx.CheckBox(dataDisplay,-1,label='Show hydrogens?')
        showHydrogen.Bind(wx.EVT_CHECKBOX, OnShowHyd)
        showHydrogen.SetValue(drawingData['showHydrogen'])
        flexSizer.Add(showHydrogen,0,wx.ALIGN_CENTER_VERTICAL)

        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(dataDisplay,-1,'Background color:'),0,wx.ALIGN_CENTER_VERTICAL)
        backColor = wcs.ColourSelect(dataDisplay, -1,colour=drawingData['backColor'],size=wx.Size(25,25))
        backColor.Bind(wcs.EVT_COLOURSELECT, OnBackColor)
        lineSizer.Add(backColor,0,wx.ALIGN_CENTER_VERTICAL)
        flexSizer.Add(lineSizer,0,)

        flexSizer.Add(wx.StaticText(dataDisplay,-1,'Hydrogen radius, A:  '),0,wx.ALIGN_CENTER_VERTICAL)
        sizeH = wx.TextCtrl(dataDisplay,-1,value='%.2f'%(drawingData['sizeH']),style=wx.TE_PROCESS_ENTER)
        sizeH.Bind(wx.EVT_TEXT_ENTER,OnSizeHatoms)
        sizeH.Bind(wx.EVT_KILL_FOCUS,OnSizeHatoms)
        flexSizer.Add(sizeH,0,wx.ALIGN_CENTER_VERTICAL)

        flexSizer.Add(wx.StaticText(dataDisplay,-1,'Bond search factor:  '),0,wx.ALIGN_CENTER_VERTICAL)
        radFactor = wx.TextCtrl(dataDisplay,value='%.2f'%(drawingData['radiusFactor']),style=wx.TE_PROCESS_ENTER)
        radFactor.Bind(wx.EVT_TEXT_ENTER,OnRadFactor)
        radFactor.Bind(wx.EVT_KILL_FOCUS,OnRadFactor)
        flexSizer.Add(radFactor,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(flexSizer,0,)

        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[1] += 26                           #compensate for status bar
        dataDisplay.SetSize(Size)
        self.dataFrame.setSizePosLeft(Size)

    def UpdateDData():
        UseList = data['Histograms']
        if UseList:
            self.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,True)
        else:
            self.dataFrame.DataMenu.Enable(G2gd.wxID_DATADELETE,False)            
        generalData = data['General']
        SGData = generalData['SGData']
        keyList = UseList.keys()
        keyList.sort()
        Indx = {}
        
        def OnShowData(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['Show'] = Obj.GetValue()
            UpdateDData()
        
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
            
        def OnCutoffVal(event):
            Obj = event.GetEventObject()
            try:
                cutoff = float(Obj.GetValue())
                if cutoff > 0:
                    UseList[Indx[Obj.GetId()]]['Cutoff'] = cutoff
            except ValueError:
                pass
            Obj.SetValue("%.3f"%(UseList[Indx[Obj.GetId()]]['Cutoff']))          #reset in case of error

        def OnSizeType(event):
            Obj = event.GetEventObject()
            UseList[Indx[Obj.GetId()]]['Size'][0] = Obj.GetValue()
            dataDisplay.Destroy()
            UpdateDData()
            
        def OnSizeRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            UseList[hist]['Size'][2][pid] = Obj.GetValue()
            
        def OnSizeVal(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            try:
                size = float(Obj.GetValue())
                if size > 0:
                    UseList[hist]['Size'][1][pid] = size
            except ValueError:
                pass
            Obj.SetValue("%.1f"%(UseList[hist]['Size'][1][pid]))          #reset in case of error
            
        def OnSizeAxis(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            axis = np.array(UseList[hist]['Size'][3])
            UseList[hist]['Size'][3][pid] = Obj.GetValue()
            new = np.array(UseList[hist]['Size'][3])
            if not np.any(new):
                UseList[hist]['Size'][3] += new-axis
            Obj.SetValue(UseList[hist]['Size'][3][pid])
                        
        def OnStrainType(event):
            Obj = event.GetEventObject()
            UseList[Indx[Obj.GetId()]]['Mustrain'][0] = Obj.GetValue()
            dataDisplay.Destroy()
            UpdateDData()
            
        def OnStrainRef(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            UseList[hist]['Mustrain'][2][pid] = Obj.GetValue()
            
        def OnStrainVal(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            try:
                strain = float(Obj.GetValue())
                if strain >= 0:
                    UseList[hist]['Mustrain'][1][pid] = strain
            except ValueError:
                pass
            Obj.SetValue("%.1f"%(UseList[hist]['Mustrain'][1][pid]))          #reset in case of error
            
        def OnStrainAxis(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            axis = np.array(UseList[hist]['Mustrain'][3])
            UseList[hist]['Mustrain'][3][pid] = Obj.GetValue()
            new = np.array(UseList[hist]['Mustrain'][3])
            if not np.any(new):
                UseList[hist]['Mustrain'][3] += new-axis
            Obj.SetValue(UseList[hist]['Mustrain'][3][pid])

        def OnMDRef(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            UseList[hist]['MDtexture'][1] = Obj.GetValue()
            
        def OnMDVal(event):
            Obj = event.GetEventObject()
            hist = Indx[Obj.GetId()]
            try:
                mdVal = float(Obj.GetValue())
                if mdVal > 0:
                    UseList[hist]['MDtexture'][0] = mdVal
            except ValueError:
                pass
            Obj.SetValue("%.3f"%(UseList[hist]['MDtexture'][0]))          #reset in case of error
            
        def OnMDAxis(event):
            Obj = event.GetEventObject()
            hist,pid = Indx[Obj.GetId()]
            axis = np.array(UseList[hist]['MDtexture'][2])
            UseList[hist]['MDtexture'][2][pid] = Obj.GetValue()
            new = np.array(UseList[hist]['MDtexture'][2])
            if not np.any(new):
                UseList[hist]['MDtexture'][2] += new-axis
            Obj.SetValue(UseList[hist]['MDtexture'][2][pid])
            
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
            Obj.SetValue("%.4f"%(UseList[Indx[Obj.GetId()]]['Extinction'][0]))          #reset in case of error
            
        def checkAxis(axis):
            if not np.any(np.array(axis)):
                return False
            return axis
            
        DData.DestroyChildren()
        dataDisplay = wx.Panel(DData)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Histogram data for '+PhaseName+':'),0,wx.ALIGN_CENTER_VERTICAL)
        for item in keyList:
            histData = UseList[item]
            mainSizer.Add(wx.StaticText(dataDisplay,-1,50*'_'))                
            mainSizer.Add((5,5),0)
            showData = wx.CheckBox(dataDisplay,label=' Show '+item)
            showData.SetValue(UseList[item]['Show'])
            Indx[showData.GetId()] = item
            showData.Bind(wx.EVT_CHECKBOX, OnShowData)
            mainSizer.Add(showData,0,wx.ALIGN_CENTER_VERTICAL)
            mainSizer.Add((0,5),0)
            if UseList[item]['Show']:
                scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
                scaleRef = wx.CheckBox(dataDisplay,label=' Scale factor: ')
                scaleRef.SetValue(UseList[item]['Scale'][1])
                Indx[scaleRef.GetId()] = item
                scaleRef.Bind(wx.EVT_CHECKBOX, OnScaleRef)
                scaleSizer.Add(scaleRef,0,wx.ALIGN_CENTER_VERTICAL)
                scaleVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,
                    '%.4f'%(UseList[item]['Scale'][0]),style=wx.TE_PROCESS_ENTER)
                Indx[scaleVal.GetId()] = item
                scaleVal.Bind(wx.EVT_TEXT_ENTER,OnScaleVal)
                scaleVal.Bind(wx.EVT_KILL_FOCUS,OnScaleVal)
                scaleSizer.Add(scaleVal,0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(scaleSizer)
                mainSizer.Add((0,5),0)
                
            if item[:4] == 'PWDR' and UseList[item]['Show']:
                cutoffSizer = wx.BoxSizer(wx.HORIZONTAL)
                cutoffSizer.Add(wx.StaticText(dataDisplay,label=' Peak cutoff ratio: '),0,wx.ALIGN_CENTER_VERTICAL)
                cutoffVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,'%.3f'%(UseList[item]['Cutoff']),
                    style=wx.TE_PROCESS_ENTER)                
                Indx[cutoffVal.GetId()] = item
                cutoffVal.Bind(wx.EVT_TEXT_ENTER,OnCutoffVal)
                cutoffVal.Bind(wx.EVT_KILL_FOCUS,OnCutoffVal)
                cutoffSizer.Add(cutoffVal,0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(cutoffSizer)
                mainSizer.Add((0,5),0)
                sizeSizer = wx.BoxSizer(wx.HORIZONTAL)
                choices = ['isotropic','uniaxial',]
                sizeType = wx.ComboBox(dataDisplay,wx.ID_ANY,value=UseList[item]['Size'][0],choices=choices,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                sizeType.Bind(wx.EVT_COMBOBOX, OnSizeType)
                Indx[sizeType.GetId()] = item
                sizeSizer.Add(sizeType)
                sizeSizer.Add((5,0),0)
                if UseList[item]['Size'][0] == 'isotropic':
                    sizeRef = wx.CheckBox(dataDisplay,label=' Cryst. size: ')
                    sizeRef.SetValue(UseList[item]['Size'][2][0])
                    Indx[sizeRef.GetId()] = [item,0]
                    sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
                    sizeSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
                    sizeVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,
                        '%.1f'%(UseList[item]['Size'][1][0]),style=wx.TE_PROCESS_ENTER)
                    Indx[sizeVal.GetId()] = [item,0]
                    sizeVal.Bind(wx.EVT_TEXT_ENTER,OnSizeVal)
                    sizeVal.Bind(wx.EVT_KILL_FOCUS,OnSizeVal)
                    sizeSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(sizeSizer)
                    mainSizer.Add((0,5),0)
                elif UseList[item]['Size'][0] == 'uniaxial':
                    sizeSizer.Add(wx.StaticText(dataDisplay,-1,' Unique axis: '),0,wx.ALIGN_CENTER_VERTICAL)
                    axes = zip(['H:','K:','L:'],UseList[item]['Size'][3],range(3))                    
                    for ax,H,i in axes:                            
                        Axis = wx.SpinCtrl(dataDisplay,wx.ID_ANY,ax,min=-3,max=3,size=wx.Size(40,20))
                        Axis.SetValue(H)
                        Indx[Axis.GetId()] = [item,i]
                        sizeSizer.Add(Axis)
                        Axis.Bind(wx.EVT_SPINCTRL, OnSizeAxis)
                    mainSizer.Add(sizeSizer)
                    mainSizer.Add((0,5),0)
                    sizeSizer = wx.BoxSizer(wx.HORIZONTAL)
                    parms = zip([' Equatorial size: ',' Axial size: '],UseList[item]['Size'][1],
                        UseList[item]['Size'][2],range(2))
                    for Pa,val,ref,id in parms:
                        sizeRef = wx.CheckBox(dataDisplay,label=Pa)
                        sizeRef.SetValue(ref)
                        Indx[sizeRef.GetId()] = [item,id]
                        sizeRef.Bind(wx.EVT_CHECKBOX, OnSizeRef)
                        sizeSizer.Add(sizeRef,0,wx.ALIGN_CENTER_VERTICAL)
                        sizeVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,'%.1f'%(val),style=wx.TE_PROCESS_ENTER)
                        Indx[sizeVal.GetId()] = [item,id]
                        sizeVal.Bind(wx.EVT_TEXT_ENTER,OnSizeVal)
                        sizeVal.Bind(wx.EVT_KILL_FOCUS,OnSizeVal)
                        sizeSizer.Add(sizeVal,0,wx.ALIGN_CENTER_VERTICAL)
                        sizeSizer.Add((5,0),0)
                    sizeSizer.Add((5,0),0)                    
                    mainSizer.Add(sizeSizer)
                
                strainSizer = wx.BoxSizer(wx.HORIZONTAL)
                choices = ['isotropic','uniaxial','generalized',]
                strainType = wx.ComboBox(dataDisplay,wx.ID_ANY,value=UseList[item]['Mustrain'][0],choices=choices,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                strainType.Bind(wx.EVT_COMBOBOX, OnStrainType)
                Indx[strainType.GetId()] = item
                strainSizer.Add(strainType)
                strainSizer.Add((5,0),0)
                if UseList[item]['Mustrain'][0] == 'isotropic':
                    strainRef = wx.CheckBox(dataDisplay,label=' microstrain: ')
                    strainRef.SetValue(UseList[item]['Mustrain'][2][0])
                    Indx[strainRef.GetId()] = [item,0]
                    strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
                    strainSizer.Add(strainRef,0,wx.ALIGN_CENTER_VERTICAL)
                    strainVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,
                        '%.4f'%(UseList[item]['Mustrain'][1][0]),style=wx.TE_PROCESS_ENTER)
                    Indx[strainVal.GetId()] = [item,0]
                    strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
                    strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
                    strainSizer.Add(strainVal,0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(strainSizer)
                    mainSizer.Add((0,5),0)
                elif UseList[item]['Mustrain'][0] == 'uniaxial':
                    strainSizer.Add(wx.StaticText(dataDisplay,-1,' Unique axis: '),0,wx.ALIGN_CENTER_VERTICAL)
                    axes = zip(['H:','K:','L:'],UseList[item]['Mustrain'][3],range(3))                    
                    for ax,H,i in axes:                            
                        Axis = wx.SpinCtrl(dataDisplay,wx.ID_ANY,ax,min=-3,max=3,size=wx.Size(40,20))
                        Axis.SetValue(H)
                        Indx[Axis.GetId()] = [item,i]
                        strainSizer.Add(Axis)
                        Axis.Bind(wx.EVT_SPINCTRL, OnStrainAxis)
                    mainSizer.Add(strainSizer)
                    mainSizer.Add((0,5),0)
                    strainSizer = wx.BoxSizer(wx.HORIZONTAL)
                    parms = zip([' Equatorial mustrain: ',' Axial mustrain: '],
                        UseList[item]['Mustrain'][1],UseList[item]['Mustrain'][2],range(2))
                    for Pa,val,ref,id in parms:
                        strainRef = wx.CheckBox(dataDisplay,label=Pa)
                        strainRef.SetValue(ref)
                        Indx[strainRef.GetId()] = [item,id]
                        strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
                        strainSizer.Add(strainRef,0,wx.ALIGN_CENTER_VERTICAL)
                        strainVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,'%.4f'%(val),style=wx.TE_PROCESS_ENTER)
                        Indx[strainVal.GetId()] = [item,id]
                        strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
                        strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
                        strainSizer.Add(strainVal,0,wx.ALIGN_CENTER_VERTICAL)
                        strainSizer.Add((5,0),0)
                    strainSizer.Add((5,0),0)                    
                    mainSizer.Add(strainSizer)
                elif UseList[item]['Mustrain'][0] == 'generalized':
                    strainSizer.Add(wx.StaticText(dataDisplay,-1,' Coefficients: '),0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(strainSizer)
                    mainSizer.Add((0,5),0)
                    Snames = G2spc.MustrainNames(SGData)
                    numb = len(Snames)
                    if len(UseList[item]['Mustrain'][1]) < numb:
                        UseList[item]['Mustrain'][1] = numb*[0.0,]
                        UseList[item]['Mustrain'][2] = numb*[False,]
                    parms = zip(Snames,UseList[item]['Mustrain'][1],UseList[item]['Mustrain'][2],range(numb))
                    strainSizer = wx.FlexGridSizer(numb%3+1,6,5,5)
                    for Pa,val,ref,id in parms:
                        strainRef = wx.CheckBox(dataDisplay,label=Pa)
                        strainRef.SetValue(ref)
                        Indx[strainRef.GetId()] = [item,id]
                        strainRef.Bind(wx.EVT_CHECKBOX, OnStrainRef)
                        strainSizer.Add(strainRef,0,wx.ALIGN_CENTER_VERTICAL)
                        strainVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,'%.4f'%(val),style=wx.TE_PROCESS_ENTER)
                        Indx[strainVal.GetId()] = [item,id]
                        strainVal.Bind(wx.EVT_TEXT_ENTER,OnStrainVal)
                        strainVal.Bind(wx.EVT_KILL_FOCUS,OnStrainVal)
                        strainSizer.Add(strainVal,0,wx.ALIGN_CENTER_VERTICAL)
                    mainSizer.Add(strainSizer)
                #MD texture  'MDtexture':[1.0,False,[0,0,1]]
                mdSizer = wx.BoxSizer(wx.HORIZONTAL)
                mdRef = wx.CheckBox(dataDisplay,label=' March-Dollase texture ratio: ')
                mdRef.SetValue(UseList[item]['MDtexture'][1])
                Indx[mdRef.GetId()] = item
                mdRef.Bind(wx.EVT_CHECKBOX, OnMDRef)
                mdSizer.Add(mdRef,0,wx.ALIGN_CENTER_VERTICAL)
                mdVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,
                    '%.3f'%(UseList[item]['MDtexture'][0]),style=wx.TE_PROCESS_ENTER)
                Indx[mdVal.GetId()] = item
                mdVal.Bind(wx.EVT_TEXT_ENTER,OnMDVal)
                mdVal.Bind(wx.EVT_KILL_FOCUS,OnMDVal)
                mdSizer.Add(mdVal,0,wx.ALIGN_CENTER_VERTICAL)
                mdSizer.Add(wx.StaticText(dataDisplay,-1,' Unique axis: '),0,wx.ALIGN_CENTER_VERTICAL)
                axes = zip(['H:','K:','L:'],UseList[item]['MDtexture'][2],range(3))                    
                for ax,H,i in axes:                            
                    Axis = wx.SpinCtrl(dataDisplay,wx.ID_ANY,ax,min=-3,max=3,size=wx.Size(40,20))
                    Axis.SetValue(H)
                    Indx[Axis.GetId()] = [item,i]
                    mdSizer.Add(Axis)
                    Axis.Bind(wx.EVT_SPINCTRL, OnMDAxis)
                mainSizer.Add(mdSizer)
                mainSizer.Add((0,5),0)
                
                #Extinction  'Extinction':[0.0,False]
                extSizer = wx.BoxSizer(wx.HORIZONTAL)
                extRef = wx.CheckBox(dataDisplay,label=' Extinction: ')
                extRef.SetValue(UseList[item]['Extinction'][1])
                Indx[extRef.GetId()] = item
                extRef.Bind(wx.EVT_CHECKBOX, OnExtRef)
                extSizer.Add(extRef,0,wx.ALIGN_CENTER_VERTICAL)
                extVal = wx.TextCtrl(dataDisplay,wx.ID_ANY,
                    '%.4f'%(UseList[item]['Extinction'][0]),style=wx.TE_PROCESS_ENTER)
                Indx[extVal.GetId()] = item
                extVal.Bind(wx.EVT_TEXT_ENTER,OnExtVal)
                extVal.Bind(wx.EVT_KILL_FOCUS,OnExtVal)
                extSizer.Add(extVal,0,wx.ALIGN_CENTER_VERTICAL)
                mainSizer.Add(extSizer)
                mainSizer.Add((0,5),0)
            elif item[:4] == 'HKLF' and UseList[item]['Show']:
                pass
        mainSizer.Add((5,5),0)

        dataDisplay.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.dataFrame)
        Size[0] = max(Size[0],300)+20
        Size[1] += 30                           #compensate for status bar
        DData.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-10)
        self.dataFrame.setSizePosLeft(Size)
        dataDisplay.SetSize(Size)
        
    def OnHklfAdd(event):
        UseList = data['Histograms']
        keyList = UseList.keys()
        TextList = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name not in keyList and 'HKLF' in name:
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                        
            dlg = wx.MultiChoiceDialog(self, 'Which new data to use?', 'Use data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result:
                        histoName = TextList[i]
                        UseList[histoName] = {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                            'Extinction':['Lorentzian','Secondary Type I',{'Eg':[0.0,False]},]}                        
                    data['Histograms'] = UseList
                    UpdateDData()
            finally:
                dlg.Destroy()
        
    def OnPwdrAdd(event):
        UseList = data['Histograms']
        keyList = UseList.keys()
        TextList = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name not in keyList and 'PWDR' in name:
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            dlg = wx.MultiChoiceDialog(self, 'Which new data to use?', 'Use data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: 
                        histoName = TextList[i]
                        UseList[histoName] = {'Histogram':histoName,'Show':False,
                            'Scale':[1.0,False],'MDtexture':[1.0,False,[0,0,1]],
                            'Size':['isotropic',[10000.,0,],[False,False],[0,0,1]],
                            'Mustrain':['isotropic',[0.0,0,],[False,False],[0,0,1]],                            
                            'Extinction':[0.0,False],'Cutoff':0.01}
                    data['Histograms'] = UseList
                    UpdateDData()
            finally:
                dlg.Destroy()
        
    def OnDataDelete(event):
        UseList = data['Histograms']
        keyList = UseList.keys()
        keyList.sort()
        DelList = []
        if UseList:
            DelList = []
            dlg = wx.MultiChoiceDialog(self, 
                'Which histogram to delete from this phase?', 'Delete histogram', 
                keyList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append(keyList[i])
                    for i in DelList:
                        del UseList[i]
                    UpdateDData()
            finally:
                dlg.Destroy()

    def FillPawleyReflectionsGrid():
        generalData = data['General']
        print 'Pawley reflections'

    def OnPageChanged(event):
        page = event.GetSelection()
        text = self.dataDisplay.GetPageText(page)
        if text == 'Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.AtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, OnAtomAdd, id=G2gd.wxID_ATOMSEDITADD)
            self.dataFrame.Bind(wx.EVT_MENU, OnAtomTestAdd, id=G2gd.wxID_ATOMSTESTADD)
            self.dataFrame.Bind(wx.EVT_MENU, OnAtomInsert, id=G2gd.wxID_ATOMSEDITINSERT)
            self.dataFrame.Bind(wx.EVT_MENU, OnAtomTestInsert, id=G2gd.wxID_ATONTESTINSERT)
            self.dataFrame.Bind(wx.EVT_MENU, AtomDelete, id=G2gd.wxID_ATOMSEDITDELETE)
            self.dataFrame.Bind(wx.EVT_MENU, AtomRefine, id=G2gd.wxID_ATOMSREFINE)
            self.dataFrame.Bind(wx.EVT_MENU, AtomModify, id=G2gd.wxID_ATOMSMODIFY)
            self.dataFrame.Bind(wx.EVT_MENU, AtomTransform, id=G2gd.wxID_ATOMSTRANSFORM)
            FillAtomsGrid()
        elif text == 'General':
            UpdateGeneral()
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        elif text == 'Data':
            self.dataFrame.SetMenuBar(self.dataFrame.DataMenu)
            self.dataFrame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2gd.wxID_PWDRADD)
            self.dataFrame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2gd.wxID_HKLFADD)
            self.dataFrame.Bind(wx.EVT_MENU, OnDataDelete, id=G2gd.wxID_DATADELETE)
            UpdateDData()            
        elif text == 'Draw Options':
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
            UpdateDrawOptions()
            G2plt.PlotStructure(self,data)
        elif text == 'Draw Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.DrawAtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomStyle, id=G2gd.wxID_DRAWATOMSTYLE)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomLabel, id=G2gd.wxID_DRAWATOMLABEL)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomColor, id=G2gd.wxID_DRAWATOMCOLOR)
            self.dataFrame.Bind(wx.EVT_MENU, ResetAtomColors, id=G2gd.wxID_DRAWATOMRESETCOLOR)
            self.dataFrame.Bind(wx.EVT_MENU, SetViewPoint, id=G2gd.wxID_DRAWVIEWPOINT)
            self.dataFrame.Bind(wx.EVT_MENU, AddSymEquiv, id=G2gd.wxID_DRAWADDEQUIV)
            self.dataFrame.Bind(wx.EVT_MENU, TransformSymEquiv, id=G2gd.wxID_DRAWTRANSFORM)
            self.dataFrame.Bind(wx.EVT_MENU, FillCoordSphere, id=G2gd.wxID_DRAWFILLCOORD)            
            self.dataFrame.Bind(wx.EVT_MENU, FillUnitCell, id=G2gd.wxID_DRAWFILLCELL)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomsDelete, id=G2gd.wxID_DRAWDELETE)
            UpdateDrawAtoms()
            G2plt.PlotStructure(self,data)
        else:
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        event.Skip()
        
    General = wx.Window(self.dataDisplay)
    self.dataDisplay.AddPage(General,'General')
    SetupGeneral()
    GeneralData = data['General']
    UpdateGeneral()

    if GeneralData['Type'] == 'Pawley':
        PawleyRefl = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(PawleyRefl,'Pawley reflections')
    else:
        DData = wx.ScrolledWindow(self.dataDisplay)
        self.dataDisplay.AddPage(DData,'Data')
        Atoms = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(Atoms,'Atoms')
        drawOptions = wx.Window(self.dataDisplay)
        self.dataDisplay.AddPage(drawOptions,'Draw Options')
        drawAtoms = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(drawAtoms,'Draw Atoms')

    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    
            
