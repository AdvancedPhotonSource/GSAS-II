#GSASII - phase data display routines
import wx
import wx.grid as wg
import matplotlib as mpl
import math
import copy
import time
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
            
def UpdatePhaseData(self,item,data,oldPage):
    
    class SymOpDialog(wx.Dialog):
        def __init__(self,parent,SGData):
            wx.Dialog.__init__(self,parent,-1,'Select symmetry operator', 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            panel = wx.Panel(self)
            self.SGData = SGData
            self.OpSelected = [0,0,0,0,0,0]
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
            self.OpSelected = [0,0,0,[],0]
            if self.SGData['SGInv']:
                self.OpSelected[0] = self.inv.GetSelection()
            if self.SGData['SGLatt'] != 'P':            
                self.OpSelected[1] = self.latt.GetSelection()
            self.OpSelected[2] = self.oprs.GetSelection()
            for i in range(3):
                self.OpSelected[3].append(float(self.cell[i].GetValue()))
            self.OpSelected[4] = self.new.GetSelection()
                
        def GetSelection(self):
            return self.OpSelected
                            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.SetReturnCode(wx.ID_OK)
            self.MakeModal(False)              
            self.Destroy()
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.SetReturnCode(wx.ID_CANCEL)
            self.MakeModal(False)              
            self.Destroy()            
        
    Atoms = []
    self.SelectedRow = 0
    
    def BookResize(event):
        w,h = self.GetSize()
        self.dataDisplay.SetSize(wx.Size(w,h))
        
    def FillGeneralGrid():
        def SetLatticeParametersStyle(SGData,table):
            clist = [1,2,3,4,5,6]
            if SGData['SGLaue'] in ['m3','m3m']:
                table[4][2] = table[4][3] = table[4][1]
                table[4][4] = table[4][5] = table[4][6] = 90.
                clist = [2,3,4,5,6]
            elif SGData['SGLaue'] in ['3R','3mR']:
                table[4][2] = table[4][3] = table[4][1]
                table[4][5] = table[4][6] = table[4][4]
                clist = [2,3,5,6]
            elif SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
                table[4][2] = table[4][1]
                table[4][4] = table[4][5] = 90.
                table[4][6] = 120.
                clist = [2,4,5,6]
            elif SGData['SGLaue'] in ['4/m','4/mmm']:
                table[4][2] = table[4][1]
                table[4][4] = table[4][5] = table[4][6] = 90.
                clist = [2,4,5,6]
            elif SGData['SGLaue'] in ['mmm']:
                table[4][4] = table[4][5] = table[4][6] = 90.
                clist = [4,5,6]
            elif SGData['SGLaue'] in ['2/m']:
                if SGData['SGUniq'] == 'a':
                    table[4][5]= table[4][6] = 90.
                    clist = [5,6]
                if SGData['SGUniq'] == 'b':
                    table[4][4]= table[4][6] = 90.
                    clist = [4,6]
                if SGData['SGUniq'] == 'c':
                    table[4][4]= table[4][5] = 90.
                    clist = [4,5]
            for c in clist:
                General.SetCellStyle(4,c,VERY_LIGHT_GREY,True)                
            
        def RefreshGeneralGrid(event):
                
            r,c =  event.GetRow(),event.GetCol()
            generalData['Name'] = table[0][0]
            self.PatternTree.SetItemText(item,generalData['Name'])
            generalData['Type'] = table[1][0]
            SpcGp = table[2][0]
            SGErr,SGData = G2spc.SpcGroup(SpcGp)
            if r == 2 and c == 0:
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    table[2][0] = generalData['SGData']['SpGrp']
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
            General.SetCellValue(4,0,str(generalData['Cell'][0]))
            for c in range(1,7):
                General.SetCellStyle(4,c,"white",False)
                generalData['Cell'][c] = float(General.GetCellValue(4,c))
            generalData['Cell'][7] = G2lat.calc_V(G2lat.cell2A(generalData['Cell'][1:7]))
            SetLatticeParametersStyle(SGData,table)
            generalData['Scale'][1] = float(General.GetCellValue(5,1))
            General.ForceRefresh()
                        
        rowLabels = ['Phase name','Phase type','Space group',
            'Lattice ',' parameters','Scale factor','Elements','No. per cell','Atom weight','','Bond radii','Angle radii']
        generalData = data['General']
        atomData = data['Atoms']
        generalData['AtomTypes'] = []
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['AtomMass'] = []
        colType = 1
        colSS = 7
        self.dataFrame.setSizePosLeft([600,350])
        if generalData['Type'] =='macromolecular':
            colType = 4
            colSS = 10
        for atom in atomData:
            if generalData['AtomTypes'].count(atom[colType]):
                generalData['NoAtoms'][atom[colType]] += atom[colSS-1]*atom[colSS+1]
            else:
                Info = G2elem.GetAtomInfo(atom[colType])
                generalData['AtomTypes'].append(Info['Symbol'])
                generalData['BondRadii'].append(Info['Drad'])
                generalData['AngleRadii'].append(Info['Arad'])
                generalData['AtomMass'].append(Info['Mass'])
                generalData['NoAtoms'][atom[colType]] = atom[colSS-1]*atom[colSS+1]
        colLabels = []
        colLabels += ['' for i in range(max(8,len(generalData['AtomTypes'])))]
        table = []
        table.append([generalData['Name'],'','','','','','','',''])      #phase name
        table.append([generalData['Type'],'','','','','','','',''])      #phase type
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        table.append([SGData['SpGrp'],'','','','','','','',''])     #space group symbol
        table.append(['refine','a    ','b    ','c    ','alpha ','beta ','gamma','volume  '])
        table.append(generalData['Cell'])                      #lattice parameters
        table.append([generalData['Scale'][0],generalData['Scale'][1],'','','','','',''])   #scale factor
        line = []
        if generalData['Type'] == 'Pawley':
            table.append(max(0.25,generalData['Pawley dmin']))
            rowLabels[6] = 'd min'
        else:
            table.append(generalData['AtomTypes']+['' for i in range(max(8,len(generalData['AtomTypes'])))]) #element list
            mass = 0.
            for i,elem in enumerate(generalData['AtomTypes']):
                mass += generalData['NoAtoms'][elem]*generalData['AtomMass'][i]
                line.append(generalData['NoAtoms'][elem])
            Volume = generalData['Cell'][7]
            table.append(line+['' for i in range(max(8,len(generalData['AtomTypes'])))]) #No. per cell
            table.append(generalData['AtomMass']+['' for i in range(max(8,len(generalData['AtomTypes'])))])  #At. wt.
            if generalData['Type'] == 'macromolecular' and mass > 0.0:
                table.append(['density',mass/(0.6022137*Volume),'Matthews coeff.',Volume/mass,'','','','',''])            
            else:
                table.append(['density',mass/(0.6022137*Volume),'','','','','','',''])
            table.append(generalData['BondRadii']+['' for i in range(max(8,len(generalData['AtomTypes'])))])
            table.append(generalData['AngleRadii']+['' for i in range(max(8,len(generalData['AtomTypes'])))])
        Types = [wg.GRID_VALUE_STRING for i in range(max(8,len(generalData['AtomTypes'])))]
        generalTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        General.SetTable(generalTable, True)
        General.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshGeneralGrid)
        General.SetMargins(0,0)
        General.SetColSize(0,100)
        General.SetColLabelSize(0)
        for c in range(max(8,len(generalData['AtomTypes']))):
            if c > 0:
                General.SetReadOnly(0,c,isReadOnly=True)
                General.SetReadOnly(1,c,isReadOnly=True)
                General.SetReadOnly(2,c,isReadOnly=True)
            General.SetReadOnly(3,c,isReadOnly=True)                         #unit cell labels
            General.SetCellAlignment(3,c,wx.ALIGN_RIGHT, wx.ALIGN_CENTRE)
            if c < 4:
                General.SetCellRenderer(4,c,wg.GridCellFloatRenderer(10,5))
                General.SetCellEditor(4,c,wg.GridCellFloatEditor(10,5))
                General.SetReadOnly(9,c,isReadOnly=True)
            else:
                General.SetCellRenderer(4,c,wg.GridCellFloatRenderer(10,3))
                General.SetCellEditor(4,c,wg.GridCellFloatEditor(10,3))
            for r in range(6,12):
                General.SetReadOnly(r,c,isReadOnly=True)
        General.SetReadOnly(4,7,isReadOnly=True)                            #cell volume - no edit
        General.SetCellEditor(1,0,wg.GridCellChoiceEditor(['nuclear','modulated',   #phase type 
            'magnetic','macromolecular','Pawley'],False))                           #- change only if no atoms
        if line:                                                    #no.of atoms not zero!
            General.SetReadOnly(1,0,isReadOnly=True)                #can't change phase type
        General.SetCellRenderer(4,0,wg.GridCellBoolRenderer())              #lattice parameters            
        General.SetCellEditor(4,0,wg.GridCellBoolEditor())
        SetLatticeParametersStyle(SGData,table)
        General.SetCellRenderer(5,1,wg.GridCellFloatRenderer(10,4))         #scale factor
        General.SetCellEditor(5,1,wg.GridCellFloatEditor(10,4))
        General.SetCellRenderer(5,0,wg.GridCellBoolRenderer())            
        General.SetCellEditor(5,0,wg.GridCellBoolEditor())
        General.SetCellRenderer(9,1,wg.GridCellFloatRenderer(8,3))
        General.SetCellRenderer(9,3,wg.GridCellFloatRenderer(8,3))
    
    def FillAtomsGrid():
        
        def RefreshAtomGrid(event):
            
            def chkUij(Uij,CSI):
                return Uij
                
            r,c =  event.GetRow(),event.GetCol()
            if r < 0 and c < 0:
                Atoms.ClearSelection()
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
            elif c < 0:                    #picked atom row
                self.SelectedRow = r
            elif Atoms.GetColLabelValue(c) in ['x','y','z']:
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
            elif Atoms.GetColLabelValue(c) == 'I/A':
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
                CSI = G2spc.GetCSuinel(atomData[r][colLabels.index('site sym')])
                value = atomData[r][c]
                iUij = CSI[0][c-colLabels.index('U11')]
                for i in range(6):
                    if iUij == CSI[0][i]:
                        atomData[r][i+colLabels.index('U11')] = value*CSI[1][i]                
            Atoms.ForceRefresh()
            
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
                PE.Destroy()
                Atoms.ForceRefresh()
            else:
                event.Skip()
        
        generalData = data['General']
        atomData = data['Atoms']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,X,XU,U,F,FX,FXU,FU",
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,4', #x,y,z,frac
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_CHOICE+":I,A",
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
            Types = [wg.GRID_VALUE_NUMBER,
                wg.GRID_VALUE_CHOICE+": ,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,MSE,HOH,UNK",
                wg.GRID_VALUE_STRING] + Types
        elif generalData['Type'] == 'modulated':
            Types += []
            colLabels += []        
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
                            
    def AtomAdd(event):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        Sytsym,Mult = G2spc.SytSym([0,0,0],SGData)
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','UNK',Sytsym,Mult,0,0,1,'',0,'I',0.10,0,0,0,0,0,0])
        elif generalData['Type'] == 'nuclear':
            atomData.append(['UNK','UNK','',0,0,0,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0])
        elif generalData['Type'] == 'magnetic':
            atomData.append(['UNK','UNK','',0,0,0,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0])
        FillAtomsGrid()            
        event.StopPropagation()
            
    def AtomInsert(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            indx = indx[0]
            atomData = data['Atoms']
            generalData = data['General']
            Ncol = Atoms.GetNumberCols()
            E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
            Sytsym,Mult = G2spc.SytSym([0,0,0],SGData)
            if generalData['Type'] == 'macromolecular':
                atomData.insert(indx,[0,'UNK','','UNK','UNK',Sytsym,Mult,0,0,1,'',0,'I',0.10,0,0,0,0,0,0])
            elif generalData['Type'] == 'nuclear':
                atomData.insert(indx,['UNK','UNK','',0,0,0,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0])
            elif generalData['Type'] == 'magnetic':
                atomData.insert(indx,['UNK','UNK','',0,0,0,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,0,0,0])
            FillAtomsGrid()            
        event.StopPropagation()
        
    def AtomDelete(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            indx.reverse()
            for ind in indx:
                del atomData[ind]                
            FillAtomsGrid()            
        event.StopPropagation()
        
    def AtomRefine(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            generalData = data['General']
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
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            c = colLabels.index('refine')
            for r in indx:
                atomData[r][c] = parms
            Atoms.ForceRefresh()                           
        
    def AtomModify(event):
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
            dlg = SymOpDialog(self,SGData)
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
                            U = np.inner(np.inner(M.T,U),M)
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
        
    def UpdateDrawing():
        print 'Drawing'
        
    def FillPawleyReflectionsGrid():
        generalData = data['General']
        
        print 'Pawley reflections'
        
    def OnPageChanged(event):
        page = event.GetSelection()
        text = self.dataDisplay.GetPageText(page)
        if text == 'Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.AtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, AtomAdd, id=G2gd.wxID_ATOMSEDITADD)
            self.dataFrame.Bind(wx.EVT_MENU, AtomInsert, id=G2gd.wxID_ATOMSEDITINSERT)
            self.dataFrame.Bind(wx.EVT_MENU, AtomDelete, id=G2gd.wxID_ATOMSEDITDELETE)
            self.dataFrame.Bind(wx.EVT_MENU, AtomRefine, id=G2gd.wxID_ATOMSREFINE)
            self.dataFrame.Bind(wx.EVT_MENU, AtomModify, id=G2gd.wxID_ATOMSMODIFY)
            self.dataFrame.Bind(wx.EVT_MENU, AtomTransform, id=G2gd.wxID_ATOMSTRANSFORM)
            FillAtomsGrid()            
        else:
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        event.Skip()
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()                    
    PhaseName = self.PatternTree.GetItemText(item)
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataFrame.SetLabel('Phase Data for '+PhaseName)
    self.dataFrame.CreateStatusBar()
    self.dataDisplay = G2gd.GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())
    
    General = G2gd.GSGrid(self.dataDisplay)
    FillGeneralGrid()
    self.dataDisplay.AddPage(General,'General')
     
    GeneralData = data['General']
    if GeneralData['Type'] == 'Pawley':
        PawleyRefl = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(PawleyRefl,'Pawley reflections')
    else:
        Atoms = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(Atoms,'Atoms')
        Drawing = wx.Window(self.dataDisplay)
        self.dataDisplay.AddPage(Drawing,'Drawing')
   
    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    