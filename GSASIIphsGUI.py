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

VERY_LIGHT_GREY = wx.Colour(235,235,235)

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
            if SGData['SGLaue'] in ['m3','m3m']:
                table[4][2] = table[4][3] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                General.SetCellStyle(4,3,"light grey",True)
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['3R','3mR']:
                table[4][2] = table[4][3] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                General.SetCellStyle(4,3,"light grey",True)
                table[4][5] = table[4][6] = table[4][4]
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['3','3m1','31m','6/m','6/mmm']:
                table[4][2] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                table[4][4] = table[4][5] = 90.
                table[4][6] = 120.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['4/m','4/mmm']:
                table[4][2] = table[4][1]
                General.SetCellStyle(4,2,"light grey",True)
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['mmm']:
                table[4][4] = table[4][5] = table[4][6] = 90.
                General.SetCellStyle(4,4,"light grey",True)
                General.SetCellStyle(4,5,"light grey",True)
                General.SetCellStyle(4,6,"light grey",True)
            elif SGData['SGLaue'] in ['2/m']:
                if SGData['SGUniq'] == 'a':
                    table[4][5]= table[4][6] = 90.
                    General.SetCellStyle(4,5,"light grey",True)
                    General.SetCellStyle(4,6,"light grey",True)
                if SGData['SGUniq'] == 'b':
                    table[4][4]= table[4][6] = 90.
                    General.SetCellStyle(4,4,"light grey",True)
                    General.SetCellStyle(4,6,"light grey",True)
                if SGData['SGUniq'] == 'c':
                    table[4][4]= table[4][5] = 90.
                    General.SetCellStyle(4,4,"light grey",True)
                    General.SetCellStyle(4,5,"light grey",True)
            
        def RefreshGeneralGrid(event):
                
            r,c =  event.GetRow(),event.GetCol()
            generalData[0] = table[0][0]
            self.PatternTree.SetItemText(item,generalData[0])
            generalData[1] = table[1][0]
            SpcGp = table[2][0]
            SGErr,SGData = G2spc.SpcGroup(SpcGp)
            if r == 2 and c == 0:
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    table[2][0] = generalData[2]['SpGrp']
                    msg = 'Space Group Error'
                    Style = wx.ICON_EXCLAMATION
                else:
                    text = G2spc.SGPrint(SGData)
                    generalData[2] = SGData
                    msg = 'Space Group Information'
                    Style = wx.ICON_INFORMATION
                Text = ''
                for line in text:
                    Text += line+'\n'
                wx.MessageBox(Text,caption=msg,style=Style)
            General.SetCellValue(4,0,str(generalData[3][0]))
            for c in range(1,7):
                General.SetCellStyle(4,c,"white",False)
                generalData[3][c] = float(General.GetCellValue(4,c))
            generalData[3][7] = G2lat.calc_V(G2lat.cell2A(generalData[3][1:7]))
            SetLatticeParametersStyle(SGData,table)
            generalData[4][1] = float(General.GetCellValue(5,1))
            General.ForceRefresh()
                        
        rowLabels = ['Phase name','Phase type','Space group',
            'Lattice ',' parameters','Scale factor','Elements','No. per cell','Atom weight','','Bond radii','Angle radii']
        generalData = data['General']
        atomData = data['Atoms']
        AtomTypes = []
        NoAtoms = {}
        BondRadii = []
        AngleRadii = []
        AtomMass = []
        colType = 1
        colSS = 7
        self.dataFrame.setSizePosLeft([600,350])
        if generalData[1] =='macromolecular':
            colType = 4
            colSS = 10
        for atom in atomData:
            if AtomTypes.count(atom[colType]):
                NoAtoms[atom[colType]] += atom[colSS-1]*atom[colSS+1]
            else:
                Info = G2elem.GetAtomInfo(atom[colType])
                AtomTypes.append(Info['Symbol'])
                BondRadii.append(Info['Drad'])
                AngleRadii.append(Info['Arad'])
                AtomMass.append(Info['Mass'])
                NoAtoms[atom[colType]] = atom[colSS-1]*atom[colSS+1]
        generalData[5:9] = [AtomTypes,NoAtoms,AtomMass,BondRadii,AngleRadii]
        colLabels = []
        colLabels += ['' for i in range(max(8,len(generalData[5])))]
        table = []
        table.append([generalData[0],'','','','','','','',''])      #phase name
        table.append([generalData[1],'','','','','','','',''])      #phase type
        E,SGData = G2spc.SpcGroup(generalData[2]['SpGrp'])
        table.append([SGData['SpGrp'],'','','','','','','',''])     #space group symbol
        table.append(['refine','a    ','b    ','c    ','alpha ','beta ','gamma','volume  '])
        table.append(generalData[3])                      #lattice parameters
        table.append([generalData[4][0],generalData[4][1],'','','','','',''])   #scale factor
        table.append(generalData[5]+['' for i in range(max(8,len(generalData[5])))]) #element list
        line = []
        mass = 0.
        for i,elem in enumerate(generalData[5]):
            mass += generalData[6][elem]*generalData[7][i]
            line.append(generalData[6][elem])
        Volume = generalData[3][7]
        table.append(line+['' for i in range(max(8,len(generalData[5])))]) #No. per cell
        table.append(generalData[7]+['' for i in range(max(8,len(generalData[5])))])  #At. wt.
        if generalData[1] == 'macromolecular' and mass > 0.0:
            table.append(['density',mass/(0.6022137*Volume),'Matthews coeff.',Volume/mass,'','','','',''])            
        else:
            table.append(['density',mass/(0.6022137*Volume),'','','','','','',''])
        table.append(generalData[8]+['' for i in range(max(8,len(generalData[5])))])
        table.append(generalData[9]+['' for i in range(max(8,len(generalData[5])))])
        Types = [wg.GRID_VALUE_STRING for i in range(max(8,len(generalData[5])))]
        generalTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
        General.SetTable(generalTable, True)
        General.Bind(wg.EVT_GRID_CELL_CHANGE, RefreshGeneralGrid)
        General.SetMargins(0,0)
        General.SetColSize(0,100)
        General.SetColLabelSize(0)
        for c in range(max(8,len(generalData[5]))):
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
                colLabel = Atoms.GetColLabelValue(c)
                if colLabel == 'x':
                    XYZ = [atomData[r][c],atomData[r][c+1],atomData[r][c+2]]
                elif colLabel == 'y':
                    XYZ = [atomData[r][c-1],atomData[r][c],atomData[r][c+1]]
                elif colLabel == 'z':
                    XYZ = [atomData[r][c-2],atomData[r][c-1],atomData[r][c]]
                if None in XYZ:
                    XYZ = [0,0,0]
                SScol = colLabels.index('site sym')
                Mulcol = colLabels.index('mult')
                E,SGData = G2spc.SpcGroup(generalData[2]['SpGrp'])
                Sytsym,Mult = G2spc.SytSym(XYZ,SGData)
                atomData[r][SScol] = Sytsym
                atomData[r][Mulcol] = Mult
                Atoms.ForceRefresh()
                    
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
            wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',wg.GRID_VALUE_FLOAT+':10,5',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_STRING,wg.GRID_VALUE_NUMBER,wg.GRID_VALUE_CHOICE+":I,A",
            wg.GRID_VALUE_FLOAT+':10,4',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',
            wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
        colLabels = ['Name','Type','refine','x','y','z','frac','site sym','mult','I/A','Uiso','U11','U22','U33','U12','U13','U23']
        if generalData[1] == 'magnetic':
            colLabels += ['Mx','My','Mz']
            Types[2] = wg.GRID_VALUE_CHOICE+": ,X,XU,U,M,MX,MXU,MU,F,FX,FXU,FU,FM,FMX,FMU,"
            Types += [
                wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_FLOAT+':10,4']
        elif generalData[1] == 'macromolecular':
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
                    Atoms.SetReadOnly(row,colIA+i,isReadOnly=True)
                    Atoms.SetCellValue(row,colIA+i,'')
            elif Atoms.GetCellValue(row,colIA) == 'A':
                Atoms.SetCellRenderer(row,colIA+1,wg.GridCellStringRenderer())
                Atoms.SetReadOnly(row,colIA+1,isReadOnly=True)
                Atoms.SetCellValue(row,colIA+1,'')
        
    def AtomAdd(event):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        if generalData[1] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','UNK','',0,0,0,0,'',0,'I',0.10,0,0,0,0,0,0])
        elif generalData[1] == 'nuclear':
            atomData.append(['UNK','UNK','',0,0,0,0,'',0,'I',0.01,0,0,0,0,0,0])
        event.StopPropagation()
        FillAtomsGrid()
            
    def AtomInsert(event):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        if generalData[1][0] == 'macromolecular':
            atomData.append([0,'UNK','','UNK','UNK','',0,0,0,0,'',0,'I',0.10,0,0,0,0,0,0])
        elif generalData[1][0] == 'nuclear':
            atomData.append(['UNK','UNK','',0,0,0,0,'',0,'I',0.01,0,0,0,0,0,0])
        event.StopPropagation()
        FillAtomsGrid()
        
    def UpdateDrawing():
        print 'Drawing'
        
    def FillPawleyReflectionsGrid():
        
        print 'Pawley reflections'
        
    def OnPageChanged(event):
        page = event.GetSelection()
        text = self.dataDisplay.GetPageText(page)
        if text == 'Atoms':
            self.dataFrame.SetMenuBar(self.dataFrame.AtomsMenu)
            self.dataFrame.Bind(wx.EVT_MENU, AtomAdd, id=G2gd.wxID_ATOMSEDITADD)
            self.dataFrame.Bind(wx.EVT_MENU, AtomInsert, id=G2gd.wxID_ATOMSEDITINSERT)
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
    if GeneralData[3] == 'Pawley':
        PawleyRefl = G2gd.GSGrid(parent=self.dataDisplay)
        self.dataDisplay.AddPage(PawleyRefl,'Pawley reflections')
        FillPawleyReflectionsGrid()
    else:
        Atoms = G2gd.GSGrid(parent=self.dataDisplay)
        FillAtomsGrid()
        self.dataDisplay.AddPage(Atoms,'Atoms')

    Drawing = wx.Window(parent=self.dataDisplay)
    self.dataDisplay.AddPage(Drawing,'Drawing')
   
    self.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
    self.dataDisplay.SetSelection(oldPage)
    