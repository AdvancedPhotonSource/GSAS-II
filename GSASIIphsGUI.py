#GSASII - phase data display routines
import wx
import wx.grid as wg
import matplotlib as mpl
import math
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
                    generalData[SGData] = SGData
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
            if r < 0:                          #on col label!
                sel = -1
                if Atoms.GetColLabelValue(c) == 'refine':
                    choice = ['F - site fraction','X - coordinates','U - thermal parameters']
                    dlg = wx.MultiChoiceDialog(self,'Select','Refinement controls',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelections()
                        parms = ''
                        for x in sel:
                            parms += choice[x][0]                            
                elif Atoms.GetColLabelValue(c) == 'I/A':
                    choice = ['Isotropic','Anisotropic']
                    dlg = wx.SingleChoiceDialog(self,'Select','Thermal Motion',choice)
                    if dlg.ShowModal() == wx.ID_OK:
                        sel = dlg.GetSelection()
                        parms = choice[sel][0]
                if sel >= 0:
                    for r in range(Atoms.GetNumberRows()):
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
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',    #Uij
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
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
        table = []
        rowLabels = []
        for i,atom in enumerate(atomData):
            table.append(atom)
            rowLabels.append(str(i+1))
        atomTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        Atoms.SetTable(atomTable, True)
        Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_SELECT_CELL, AtomTypeSelect)
        Atoms.SetMargins(0,0)
        Atoms.AutoSizeColumns(True)
        colType = colLabels.index('Type')
        colSS = colLabels.index('site sym')
        colIA = colLabels.index('I/A')
        for row in range(Atoms.GetNumberRows()):
            Atoms.SetReadOnly(row,colSS,True)                         #site sym
            Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
            if Atoms.GetCellValue(row,colIA) == 'I':
                for i in range(2,8):
                    Atoms.SetCellRenderer(row,colIA+i,wg.GridCellStringRenderer())
                    Atoms.SetCellValue(row,colIA+i,'')
                    Atoms.SetCellStyle(row,colIA+i,VERY_LIGHT_GREY,True)
            elif Atoms.GetCellValue(row,colIA) == 'A':
                CSI = G2spc.GetCSuinel(atomData[row][colLabels.index('site sym')])
                Atoms.SetCellRenderer(row,colIA+1,wg.GridCellStringRenderer())
                Atoms.SetCellStyle(row,colIA+1,VERY_LIGHT_GREY,True)
                Atoms.SetCellValue(row,colIA+1,'')
                for i in range(6):
                    ci = colIA+i+2
                    Atoms.SetCellStyle(row,ci,VERY_LIGHT_GREY,True)
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
            FillAtomsGrid()            
        else:
            self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
        event.Skip()
        
    if self.dataDisplay:
        self.dataDisplay.Destroy()                    
    PhaseName = self.PatternTree.GetItemText(item)
    self.dataFrame.SetMenuBar(self.dataFrame.BlankMenu)
    self.dataFrame.SetLabel('Phase Data for '+PhaseName)
    self.dataDisplay = G2gd.GSNoteBook(parent=self.dataFrame,size=self.dataFrame.GetClientSize())
    
    General = G2gd.GSGrid(parent=self.dataDisplay)
    FillGeneralGrid()
    self.dataDisplay.AddPage(General,'General')
     
    GeneralData = data['General']
    if GeneralData['Type'] == 'Pawley':
        PawleyRefl = G2gd.GSGrid(parent=self.dataDisplay)
        self.dataDisplay.AddPage(PawleyRefl,'Pawley reflections')
    else:
        Atoms = G2gd.GSGrid(parent=self.dataDisplay)
        self.dataDisplay.AddPage(Atoms,'Atoms')
        Drawing = wx.Window(parent=self.dataDisplay)
        self.dataDisplay.AddPage(Drawing,'Drawing')
   
    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    