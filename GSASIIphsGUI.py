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

class SymOpDialog(wx.Dialog):
    def __init__(self,parent,SGData,New=True):
        wx.Dialog.__init__(self,parent,-1,'Select symmetry operator',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        self.SGData = SGData
        self.New = New
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
        self.OpSelected = [0,0,0,[]]
        if self.SGData['SGInv']:
            self.OpSelected[0] = self.inv.GetSelection()
        if self.SGData['SGLatt'] != 'P':
            self.OpSelected[1] = self.latt.GetSelection()
        self.OpSelected[2] = self.oprs.GetSelection()
        for i in range(3):
            self.OpSelected[3].append(float(self.cell[i].GetValue()))
        if self.New:
            self.OpSelected.append(self.new.GetSelection())

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
    self.SelectedRow = 0

#    def BookResize(event):
#        w,h = self.GetSize()
#        self.dataDisplay.SetSize(wx.Size(w,h))
#        
    def UpdateGeneral():
        generalData = data['General']
        atomData = data['Atoms']
        generalData['AtomTypes'] = []
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['vdWRadii'] = []
        generalData['AtomMass'] = []
        generalData['Color'] = []
        generalData['Myself'] = self
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

    def FillGeneralGrid():
        rowLabels = ['Phase name','Phase type','Space group',
            'Lattice ',' parameters','Scale factor','Density','Elements','No. per cell',
            'Atom weight','Bond radii','Angle radii','vdw radii','Color']
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
            if r == 0:
                self.G2plotNB.Rename(oldName,generalData['Name'])
            elif r == 2 and c == 0:
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

        UpdateGeneral()
        generalData = data['General']
        self.dataFrame.setSizePosLeft([750,340])
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
            mass = 0.
            for i,elem in enumerate(generalData['AtomTypes']):
                mass += generalData['NoAtoms'][elem]*generalData['AtomMass'][i]
            Volume = generalData['Cell'][7]
            if generalData['Type'] == 'macromolecular' and mass > 0.0:
                table.append([mass/(0.6022137*Volume),'Matthews coeff.',Volume/mass,'','','','','',''])
            else:
                table.append([mass/(0.6022137*Volume),'','','','','','','',''])
            for i,elem in enumerate(generalData['AtomTypes']):
                line.append(generalData['NoAtoms'][elem])
            table.append(generalData['AtomTypes']+['' for i in range(max(8,len(generalData['AtomTypes'])))]) #element list
            table.append(line+['' for i in range(max(8,len(generalData['AtomTypes'])))]) #No. per cell
            table.append(generalData['AtomMass']+['' for i in range(max(8,len(generalData['AtomTypes'])))])  #At. wt.
            table.append(generalData['BondRadii']+['' for i in range(max(8,len(generalData['AtomTypes'])))])
            table.append(generalData['AngleRadii']+['' for i in range(max(8,len(generalData['AtomTypes'])))])
            table.append(generalData['vdWRadii']+['' for i in range(max(8,len(generalData['AtomTypes'])))])
            table.append(['','','','','','','',''])                        #contains colors
        Types = [wg.GRID_VALUE_STRING for i in range(max(8,len(generalData['AtomTypes'])))]
        generalTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        General.SetTable(generalTable, True)
        General.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshGeneralGrid)
        General.SetMargins(0,0)
        General.SetColSize(0,100)
        General.SetColLabelSize(0)
        attr = wg.GridCellAttr()
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
            else:
                General.SetCellRenderer(4,c,wg.GridCellFloatRenderer(10,3))
                General.SetCellEditor(4,c,wg.GridCellFloatEditor(10,3))
            for r in range(6,13):
                General.SetReadOnly(r,c,isReadOnly=True)
        r = rowLabels.index('Color')
        for c in range(len(generalData['AtomTypes'])):
            General.SetCellBackgroundColour(r,c,generalData['Color'][c])

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
        General.SetCellRenderer(6,0,wg.GridCellFloatRenderer(8,3))
        General.SetCellRenderer(6,2,wg.GridCellFloatRenderer(8,3))

    def FillAtomsGrid():

        self.dataFrame.setSizePosLeft([700,300])
        generalData = data['General']
        atomData = data['Atoms']
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
                atomData[r][c] = Atoms.GetCellValue(r,c)
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
                UpdateGeneral()
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
                UpdateGeneral()
                FillAtomsGrid()
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
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','H','',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0])
        elif generalData['Type'] == 'nuclear':
            atomData.append(['UNK','H','',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0])
        elif generalData['Type'] == 'magnetic':
            atomData.append(['UNK','H','',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0])
        UpdateGeneral()

    def OnAtomInsert(event):
        AtomInsert(0,0,0)
        FillAtomsGrid()
        event.StopPropagation()
        
    def OnAtomTestInsert(event):
        try:
            drawData = data['Drawing']
            x,y,z = drawData['testPos']
            AtomAdd(x,y,z)
        except:
            AtomAdd(0,0,0)
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
            if generalData['Type'] == 'macromolecular':
                atomData.insert(indx,[0,'UNK','','UNK','UNK','',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0])
            elif generalData['Type'] == 'nuclear':
                atomData.insert(indx,['UNK','UNK','',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0])
            elif generalData['Type'] == 'magnetic':
                atomData.insert(indx,['UNK','UNK','',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,0,0,0])
            UpdateGeneral()

    def AtomDelete(event):
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            indx.reverse()
            for ind in indx:
                atom = atomData[ind]
                del atomData[ind]
            FillAtomsGrid()
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
                    print M,T
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
            'Zclip':50.0,'cameraPos':50.,'pickItem':'Atoms','showBadContacts':False,
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
        if not drawingData['Atoms']:
            for atom in atomData:
                if generalData['Type'] == 'nuclear':
                    drawingData['Atoms'].append(atom[:2]+atom[3:6]+['1',]+['lines',]+
                        ['',]+atom[9:17]+[[]])
                    drawingData['atomPtrs'] = [2,1,6]         #x, type & style
                elif generalData['Type'] == 'macromolecular':
                    try:
                        oneLetter = AA3letter.index(atom[1])
                    except ValueError:
                        oneLetter = -1
                    drawingData['Atoms'].append([atom[1].strip()+atom[0],]+
                        [AA1letter[oneLetter]+atom[0],]+atom[2:5]+
                        atom[6:9]+['1',]+['lines',]+['',]+atom[12:20]+[[]])
                    drawingData['atomPtrs'] = [5,4,9]         #x, type & style
                elif generalData['Type'] == 'magnetic':
                    drawingData['Atoms'].append(atom[:2]+atom[3:6]+['lines',]+['',]+atom[9:20]+[[]])
#            elif generalData['Type'] == 'modulated':
#                ?????   for future
            data['Drawing'] = drawingData

    def UpdateDrawAtoms():
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
        atomData = drawingData['Atoms']
        Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_STRING,
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',    #x,y,z
            wg.GRID_VALUE_STRING,wg.GRID_VALUE_CHOICE+": ,lines,vdW balls,sticks,balls & sticks,ellipsoids,polyhedra",
            wg.GRID_VALUE_CHOICE+": ,type,name,number",wg.GRID_VALUE_STRING,]
        styleChoice = [' ','lines','vdW balls','sticks','balls & sticks','ellipsoids','polyhedra']
        labelChoice = [' ','type','name','number']
        colLabels = ['Name','Type','x','y','z','Sym Op','Style','Label','I/A']
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
        drawAtoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        drawAtoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        indx = drawingData['selectedAtoms']
        if indx:
            for r in range(len(atomData)):
                if r in indx:
                    drawAtoms.SelectRow(r)
        for c in range(len(colLabels)):
            if colLabels[c] not in ['Style','Label']:
                attr = wg.GridCellAttr()                #needs to be here - gets lost if outside loop!
                attr.SetReadOnly(True)
                attr.SetBackgroundColour(VERY_LIGHT_GREY)
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
                    Inv,Cent,Opr,Cell = dlg.GetSelection()
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
                    Inv,Cent,Opr,Cell = dlg.GetSelection()
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
                        if 0 < dist <= 0.85*sumR:
                            if noDuplicate(xyz,atomData):
                                newAtom = atomB[:]
                                newAtom[cx:cx+3] = xyz
                                atomData.append(newAtom)
            data['Drawing']['Atoms'] = atomData
            UpdateDrawAtoms()
            FindBonds()
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
        Types = []
        Radii = []
        for atom in atomData:
            Atoms.append(np.array(atom[cx:cx+3]))
            Types.append(atom[cs])
            Radii.append(radii[atomTypes.index(atom[ct])])
        Atoms = np.array(Atoms)
        Radii = np.array(Radii)
        IATR = zip(Indx,Atoms,Types,Radii)
        for atomA in IATR:
            if atomA[2] in ['lines','sticks','ellipsoids','balls & sticks','polyhedra']:
                Dx = Atoms-atomA[1]
                dist = ma.masked_less(np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0)),0.5) #gets rid of self & disorder "bonds" < 0.5A
                sumR = atomA[3]+Radii
                IndB = ma.nonzero(ma.masked_greater(dist-0.85*sumR,0.))                 #get indices of bonded atoms
                i = atomA[0]
                for j in IndB[0]:
                    if j > i:
                        if Types[i] == 'polyhedra':
                            atomData[i][-1].append(np.inner(Amat,Dx[j]))
                        elif Types[j] != 'polyhedra':
                            atomData[i][-1].append(Dx[j]*Radii[i]/sumR[j])
                            atomData[j][-1].append(-Dx[j]*Radii[j]/sumR[j])

    def DrawAtomsDelete(event):   
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        if indx:
            atomData = data['Drawing']['Atoms']
            indx.reverse()
            for ind in indx:
                atom = atomData[ind]
                del atomData[ind]
            UpdateDrawAtoms()
            G2plt.PlotStructure(self,data)
        event.StopPropagation()

    def UpdateDrawOptions():
        import copy
        import wx.lib.colourselect as wcs
        self.dataFrame.setSizePosLeft([300,430])
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
            G2plt.PlotStructure(self,data)

        def OnSizeHatoms(event):
            try:
                value = max(0.1,min(1.2,float(sizeH.GetValue())))
            except ValueError:
                value = 0.5
            drawingData['sizeH'] = value
            sizeH.SetValue("%.2f"%(value))
            G2plt.PlotStructure(self,data)

        def OnPickItem(event):
            drawingData['pickItem'] = pickChoice[pickItem.GetSelection()]

        dataDisplay = wx.Panel(drawOptions)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(dataDisplay,-1,'Drawing controls:'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        
        slopSizer = wx.BoxSizer(wx.HORIZONTAL)
        slideSizer = wx.FlexGridSizer(6,2,5,0)
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
        slopSizer.SetMinSize(wx.Size(300,180))
        mainSizer.Add(slopSizer,1,wx.EXPAND)

        flexSizer = wx.FlexGridSizer(6,2,5,0)
        flexSizer.Add(wx.StaticText(dataDisplay,-1,'View Point:  '),0,wx.ALIGN_CENTER_VERTICAL)
        VP = drawingData['viewPoint'][0]
        viewPoint = wx.TextCtrl(dataDisplay,value='%.3f, %.3f, %.3f'%(VP[0],VP[1],VP[2]),
            style=wx.TE_READONLY,size=wx.Size(120,20),name='viewPoint')
        viewPoint.SetBackgroundColour(VERY_LIGHT_GREY)
        flexSizer.Add(viewPoint,0,wx.ALIGN_CENTER_VERTICAL)
        
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        lineSizer.Add(wx.StaticText(dataDisplay,-1,'Background color:'),0,wx.ALIGN_CENTER_VERTICAL)
        backColor = wcs.ColourSelect(dataDisplay, -1,colour=drawingData['backColor'],size=wx.Size(25,25))
        backColor.Bind(wcs.EVT_COLOURSELECT, OnBackColor)
        lineSizer.Add(backColor,0,wx.ALIGN_CENTER_VERTICAL)
        flexSizer.Add(lineSizer,0,)

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

        flexSizer.Add(wx.StaticText(dataDisplay,-1,'Hydrogen radius, A:  '),0,wx.ALIGN_CENTER_VERTICAL)
        sizeH = wx.TextCtrl(dataDisplay,-1,value='%.2f'%(drawingData['sizeH']),style=wx.TE_PROCESS_ENTER)
        sizeH.Bind(wx.EVT_TEXT_ENTER,OnSizeHatoms)
        sizeH.Bind(wx.EVT_KILL_FOCUS,OnSizeHatoms)
        flexSizer.Add(sizeH,0,wx.ALIGN_CENTER_VERTICAL)

        flexSizer.Add(wx.StaticText(dataDisplay,-1,'Pick items on drawing by:  '),0,wx.ALIGN_CENTER_VERTICAL)
        pickItem = wx.Choice(dataDisplay,-1,choices=pickChoice)
        pickItem.Bind(wx.EVT_CHOICE, OnPickItem)
        pickItem.SetSelection(pickChoice.index(drawingData['pickItem']))
        flexSizer.Add(pickItem,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(flexSizer,0,)
#        mainSizer.SetMinSize(wx.Size(300,340))          #to get sliders long enough

        dataDisplay.SetSizer(mainSizer)
        self.dataFrame.SetSize(dataDisplay.Fit())

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
            FillGeneralGrid()
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        elif text == 'Draw Options':
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
            UpdateDrawOptions()
            G2plt.PlotStructure(self,data)
        elif text == 'Draw Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.DrawAtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomStyle, id=G2gd.wxID_DRAWATOMSTYLE)
            self.dataFrame.Bind(wx.EVT_MENU, DrawAtomLabel, id=G2gd.wxID_DRAWATOMLABEL)
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
        drawOptions = wx.Window(self.dataDisplay)
        self.dataDisplay.AddPage(drawOptions,'Draw Options')
        drawAtoms = G2gd.GSGrid(self.dataDisplay)
        self.dataDisplay.AddPage(drawAtoms,'Draw Atoms')

    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    
            
