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

'''
import os.path
import wx
import wx.grid as wg
import wx.lib.gridmovers as wgmove
import wx.wizard as wz
import wx.lib.scrolledpanel as wxscroll
import matplotlib as mpl
import math
import copy
import time
import sys
import random as ran
import cPickle
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIddataGUI as G2ddG
import GSASIIplot as G2plt
import GSASIIgrid as G2gd
import GSASIIIO as G2IO
import GSASIIstrMain as G2stMn
import GSASIImath as G2mth
import GSASIIpwd as G2pwd
import GSASIIpy3 as G2py3
import GSASIIobj as G2obj
import numpy as np
import numpy.linalg as nl
import numpy.ma as ma

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

def SetPhaseWindow(mainFrame,phasePage,mainSizer):
    phasePage.SetSizer(mainSizer)
    Size = mainSizer.GetMinSize()
    Size[0] += 40
#    Size[1] = 500
    Size[1] = min(Size[1]+ 150,500) 
    phasePage.SetSize(Size)
    phasePage.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
    mainFrame.setSizePosLeft(Size)
    
def UpdatePhaseData(G2frame,Item,data,oldPage):
    '''Create the data display window contents when a phase is clicked on
    in the main (data tree) window.
    Called only from :meth:`GSASIIgrid.MovePatternTreeToGrid`,
    which in turn is called from :meth:`GSASII.GSASII.OnPatternTreeSelChanged`
    when a tree item is selected.

    :param wx.frame G2frame: the main GSAS-II frame object
    :param wx.TreeItemId Item: the tree item that was selected
    :param dict data: all the information on the phase in a dictionary
    :param int oldPage: This sets a tab to select when moving
      from one phase to another, in which case the same tab is selected
      to display first. This is set only when the previous data tree
      selection is a phase, if not the value is None. The default action
      is to bring up the General tab.

    '''

    # UpdatePhaseData execution continues below
    
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
            generalData['Map'] = {}
            generalData['Map'].update(mapDefault)
        if 'Flip' not in generalData:
            generalData['Flip'] = {'RefList':'','Resolution':0.5,'Norm element':'None',
                'k-factor':0.1,'k-Max':20.}
        if 'doPawley' not in generalData:
            generalData['doPawley'] = False
        if 'Pawley dmin' not in generalData:
            generalData['Pawley dmin'] = 1.0
        if 'Pawley neg wt' not in generalData:
            generalData['Pawley neg wt'] = 0.0
        if 'Algolrithm' in generalData.get('MCSA controls',{}) or \
            'MCSA controls' not in generalData:
            generalData['MCSA controls'] = {'Data source':'','Annealing':[50.,0.001,50],
            'dmin':2.0,'Algorithm':'log','Jump coeff':[0.95,0.5],'boltzmann':1.0,
            'fast parms':[1.0,1.0,1.0],'log slope':0.9,'Cycles':1,'Results':[],'newDmin':True}
        if 'AtomPtrs' not in generalData:
            generalData['AtomPtrs'] = [3,1,7,9]
            if generalData['Type'] =='macromolecular':
                generalData['AtomPtrs'] = [6,4,10,12]
        if generalData['Type'] in ['modulated','magnetic',]: 
            if 'Super' not in generalData:
                generalData['Super'] = 1
                generalData['SuperVec'] = [[0,0,.1],False,4]
                generalData['SSGData'] = {}
            if '4DmapData' not in generalData:
                generalData['4DmapData'] = {}
                generalData['4DmapData'].update(mapDefault)
                generalData['4DmapData'].update({'MapType':'Fobs'})
# end of patches
        cx,ct,cs,cia = generalData['AtomPtrs']
        generalData['NoAtoms'] = {}
        generalData['BondRadii'] = []
        generalData['AngleRadii'] = []
        generalData['vdWRadii'] = []
        generalData['AtomMass'] = []
        generalData['Color'] = []
        generalData['Mydir'] = G2frame.dirname
        badList = {}
        for atom in atomData:
            atom[ct] = atom[ct].lower().capitalize()              #force to standard form
            if generalData['AtomTypes'].count(atom[ct]):
                generalData['NoAtoms'][atom[ct]] += atom[cs-1]*float(atom[cs+1])
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
                    generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]]['Mass'])
                else:
                    generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                    generalData['AtomMass'].append(Info['Mass'])
                generalData['NoAtoms'][atom[ct]] = atom[cs-1]*float(atom[cs+1])
                generalData['Color'].append(Info['Color'])
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
       

################################################################################
##### General phase routines
################################################################################

    def UpdateGeneral():
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
        # UpdateGeneral execution starts here
        phaseTypes = ['nuclear','modulated','magnetic','macromolecular']
        SetupGeneral()
        generalData = data['General']
        Map = generalData['Map']
        Flip = generalData['Flip']
        MCSAdata = generalData['MCSA controls']  
        PWDR = any(['PWDR' in item for item in data['Histograms'].keys()])
        # UpdateGeneral execution continues below
        
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
                oldName = generalData['Name']
                phaseRIdList,usedHistograms = G2frame.GetPhaseInfofromTree()
                phaseNameList = usedHistograms.keys() # phase names in use
                newName = NameTxt.GetValue().strip()
                if newName and newName != oldName:
                    newName = G2obj.MakeUniqueLabel(newName,phaseNameList)             
                    generalData['Name'] = newName
                    G2frame.G2plotNB.Rename(oldName,generalData['Name'])
                    G2frame.dataFrame.SetLabel('Phase Data for '+generalData['Name'])
                    G2frame.PatternTree.SetItemText(Item,generalData['Name'])
                    # change phase name key in Reflection Lists for each histogram
                    for hist in data['Histograms']:
                        ht = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,hist)
                        rt = G2gd.GetPatternTreeItemId(G2frame,ht,'Reflection Lists')
                        if not rt: continue
                        RfList = G2frame.PatternTree.GetItemPyData(rt)
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
                    if generalData['Type'] in ['modulated',]:
                        if 'SuperSg' not in generalData:
                            generalData['SuperSg'] = SetDefaultSSsymbol()
                            generalData['SSGData'] = G2spc.SSpcGroup(generalData['SGData'],generalData['SuperSg'])[1]
                    wx.CallAfter(UpdateGeneral)
                else:
                    G2frame.ErrorDialog('Phase type change error','Can change phase type only if there are no atoms')
                    TypeTxt.SetValue(generalData['Type'])                
                
            def OnSpaceGroup(event):
                Flds = SGTxt.GetValue().split()
                #get rid of extra spaces between fields first
                for fld in Flds: fld = fld.strip()
                SpcGp = ' '.join(Flds)
                # try a lookup on the user-supplied name
                SpGrpNorm = G2spc.StandardizeSpcName(SpcGp)
                if SpGrpNorm:
                    SGErr,SGData = G2spc.SpcGroup(SpGrpNorm)
                else:
                    SGErr,SGData = G2spc.SpcGroup(SpcGp)
                if SGErr:
                    text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                    SGTxt.SetValue(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Error'
                    Style = wx.ICON_EXCLAMATION
                    Text = '\n'.join(text)
                    wx.MessageBox(Text,caption=msg,style=Style)
                else:
                    text,table = G2spc.SGPrint(SGData)
                    generalData['SGData'] = SGData
                    SGTxt.SetValue(generalData['SGData']['SpGrp'])
                    msg = 'Space Group Information'
                    G2gd.SGMessageBox(General,msg,text,table).Show()
                if generalData['Type'] in ['modulated',]:
                    generalData['SuperSg'] = SetDefaultSSsymbol()
                    generalData['SSGData'] = G2spc.SSpcGroup(generalData['SGData'],generalData['SuperSg'])[1]
                wx.CallAfter(UpdateGeneral)
                
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
            SGTxt = wx.TextCtrl(General,-1,value=generalData['SGData']['SpGrp'],style=wx.TE_PROCESS_ENTER)
            SGTxt.Bind(wx.EVT_TEXT_ENTER,OnSpaceGroup)
            nameSizer.Add(SGTxt,0,WACV)
            return nameSizer
            
        def CellSizer():
            
            cellGUIlist = [[['m3','m3m'],4,zip([" Unit cell: a = "," Vol = "],["%.5f","%.3f"],[True,False],[0,0])],
            [['3R','3mR'],6,zip([" a = "," alpha = "," Vol = "],["%.5f","%.3f","%.3f"],[True,True,False],[0,3,0])],
            [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],6,zip([" a = "," c = "," Vol = "],["%.5f","%.5f","%.3f"],[True,True,False],[0,2,0])],
            [['mmm'],8,zip([" a = "," b = "," c = "," Vol = "],["%.5f","%.5f","%.5f","%.3f"],
                [True,True,True,False],[0,1,2,0])],
            [['2/m'+'a'],10,zip([" a = "," b = "," c = "," alpha = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,5,0])],
            [['2/m'+'b'],10,zip([" a = "," b = "," c = "," beta = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,4,0])],
            [['2/m'+'c'],10,zip([" a = "," b = "," c = "," gamma = "," Vol = "],
                ["%.5f","%.5f","%.5f","%.3f","%.3f"],[True,True,True,True,False],[0,1,2,3,0])],
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
            
            def OnIsotope(event):   #how can I update Atom weight on isotope change?
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                isotope = Obj.GetValue()
                generalData['Isotope'][item] = isotope
                indx = generalData['AtomTypes'].index(item)
                data['General']['AtomMass'][indx] = generalData['Isotopes'][item][isotope]['Mass']
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
                choices = generalData['Isotopes'][elem].keys()
                isoSel = wx.ComboBox(General,-1,value=generalData['Isotope'][elem],choices=choices,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
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
            return elemSizer
        
        def DenSizer():
            
            mass = G2mth.getMass(generalData)
            density,mattCoeff = G2mth.getDensity(generalData)
            denSizer = wx.BoxSizer(wx.HORIZONTAL)
            denSizer.Add(wx.StaticText(General,-1,' Density: '),0,WACV)
            denTxt = wx.TextCtrl(General,-1,'%.3f'%(density),style=wx.TE_READONLY)
            denTxt.SetBackgroundColour(VERY_LIGHT_GREY)
            denSizer.Add(denTxt,0,WACV)
            mattTxt = None        
            if generalData['Type'] == 'macromolecular' and mass > 0.0:
                denSizer.Add(wx.StaticText(General,-1,' Matthews coeff.: '),
                    0,WACV)
                mattTxt = wx.TextCtrl(General,-1,'%.3f'%(mattCoeff),style=wx.TE_READONLY)
                mattTxt.SetBackgroundColour(VERY_LIGHT_GREY)
                denSizer.Add(mattTxt,0,WACV)
            return denSizer,denTxt,mattTxt
            
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
                pawlVal.SetValue("%.5f"%(generalData['Pawley dmin']))          #reset in case of error                
            
            def OnPawleyNegWt(event):
                try:
                    wt = float(pawlNegWt.GetValue())
                    if 0. <= wt <= 1.:
                        generalData['Pawley neg wt'] = wt
                except ValueError:
                    pass
                pawlNegWt.SetValue("%.4f"%(generalData['Pawley neg wt']))          #reset in case of error                

            pawleySizer = wx.BoxSizer(wx.HORIZONTAL)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley controls: '),0,WACV)
            pawlRef = wx.CheckBox(General,-1,label=' Do Pawley refinement?')
            pawlRef.SetValue(generalData['doPawley'])
            pawlRef.Bind(wx.EVT_CHECKBOX,OnPawleyRef)
            pawleySizer.Add(pawlRef,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley dmin: '),0,WACV)
            pawlVal = wx.TextCtrl(General,value='%.5f'%(generalData['Pawley dmin']),style=wx.TE_PROCESS_ENTER)
            pawlVal.Bind(wx.EVT_TEXT_ENTER,OnPawleyVal)        
            pawlVal.Bind(wx.EVT_KILL_FOCUS,OnPawleyVal)
            pawleySizer.Add(pawlVal,0,WACV)
            pawleySizer.Add(wx.StaticText(General,label=' Pawley neg. wt.: '),0,WACV)
            pawlNegWt = wx.TextCtrl(General,value='%.4f'%(generalData['Pawley neg wt']),style=wx.TE_PROCESS_ENTER)
            pawlNegWt.Bind(wx.EVT_TEXT_ENTER,OnPawleyNegWt)        
            pawlNegWt.Bind(wx.EVT_KILL_FOCUS,OnPawleyNegWt)
            pawleySizer.Add(pawlNegWt,0,WACV)
            return pawleySizer
            
        def ModulatedSizer(name):
            
            def OnSuperGp(event):
                SSymbol = superGp.GetValue()
                E,SSGData = G2spc.SSpcGroup(generalData['SGData'],SSymbol)
                if SSGData:
                    Vec = generalData['SuperVec'][0]     #(3+1) only
                    modSymb = SSGData['modSymb']
                    generalData['SuperVec'][0] = G2spc.SSGModCheck(Vec,modSymb)[0]
                    text,table = G2spc.SSGPrint(generalData['SGData'],SSGData)
                    generalData['SSGData'] = SSGData
                    generalData['SuperSg'] = SSymbol
                    msg = 'Superspace Group Information'
                    G2gd.SGMessageBox(General,msg,text,table).Show()
                else:
                    text = [E+'\nSuperspace Group set to previous']
                    superGp.SetValue(generalData['SuperSg'])
                    msg = 'Superspace Group Error'
                    Style = wx.ICON_EXCLAMATION
                    Text = '\n'.join(text)
                    wx.MessageBox(Text,caption=msg,style=Style)
                wx.CallAfter(UpdateGeneral)                
            
            def OnDim(event):
                generalData['Super'] = dim.GetValue()
                wx.CallAfter(UpdateGeneral)
                
            def OnVec(event):
                Obj = event.GetEventObject()
                ind = Indx[Obj.GetId()]
                val = Obj.GetValue()
                try:
                    val = min(1.0,max(0.0,float(val)))
                except ValueError:
                    val = generalData['SuperVec'][0][ind]
                generalData['SuperVec'][0][ind] = val
                Obj.SetValue('%.4f'%(generalData['SuperVec'][0][ind])) 
                
            def OnVecRef(event):
                generalData['SuperVec'][1] = Ref.GetValue()
                
            def OnMax(event):
                generalData['SuperVec'][2] = int(Max.GetValue())
            
            Indx = {}
            ssSizer = wx.BoxSizer(wx.VERTICAL)
            modSizer = wx.BoxSizer(wx.HORIZONTAL)
            modSizer.Add(wx.StaticText(General,label=' '+name.capitalize()+' structure controls: '),0,WACV)
            modSizer.Add(wx.StaticText(General,label=' Superspace group: '+generalData['SGData']['SpGrp']),0,WACV)
            SSChoice = G2spc.ssdict.get(generalData['SGData']['SpGrp'],[])
            if SSChoice:
                superGp = wx.ComboBox(General,value=generalData['SuperSg'],choices=SSChoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                superGp.Bind(wx.EVT_COMBOBOX,OnSuperGp)
            else:   #nonstandard space group symbol not in my dictionary
                superGp = wx.TextCtrl(General,value=generalData['SuperSg'],style=wx.TE_PROCESS_ENTER)
                superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)                        
            modSizer.Add(superGp,0,WACV)
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
            vec = generalData['SuperVec'][0]
            for i,[val,show] in enumerate(zip(generalData['SuperVec'][0],ifShow)):
                if show:
                    modVal = wx.TextCtrl(General,value=('%.4f'%(val)),
                        size=wx.Size(50,20),style=wx.TE_PROCESS_ENTER)
                    modVal.Bind(wx.EVT_TEXT_ENTER,OnVec)        
                    modVal.Bind(wx.EVT_KILL_FOCUS,OnVec)
                    vecSizer.Add(modVal,0,WACV)
                    Indx[modVal.GetId()] = i
                else:
                    modVal = wx.TextCtrl(General,value=('%.3f'%(val)),
                        size=wx.Size(50,20),style=wx.TE_READONLY)
                    modVal.SetBackgroundColour(VERY_LIGHT_GREY)
                    vecSizer.Add(modVal,0,WACV)
            Ref = wx.CheckBox(General,label='Refine?')
            Ref.SetValue(generalData['SuperVec'][1])
            Ref.Bind(wx.EVT_CHECKBOX, OnVecRef)
            vecSizer.Add(Ref,0,WACV)
            ssSizer.Add(vecSizer)
            return ssSizer
            
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
                    if 10.0 <= res <= 100.:
                        Map['cutOff'] = res
                except ValueError:
                    pass
                cutOff.SetValue("%.1f"%(Map['cutOff']))          #reset in case of error
            
            #patch
            if 'cutOff' not in Map:
                Map['cutOff'] = 100.0
            mapTypes = ['Fobs','Fcalc','delt-F','2*Fo-Fc','Omit','2Fo-Fc Omit','Patterson']
            refList = data['Histograms'].keys()
            if not generalData['AtomTypes']:
                 mapTypes = ['Patterson',]
                 Map['MapType'] = 'Patterson'
            mapSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Fourier map controls: Map type: '),0,WACV)
            mapType = wx.ComboBox(General,-1,value=Map['MapType'],choices=mapTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            mapType.Bind(wx.EVT_COMBOBOX,OnMapType)
            lineSizer.Add(mapType,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' Reflection set from: '),0,WACV)
            refList = wx.ComboBox(General,-1,value=Map['RefList'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,WACV)
            mapSizer.Add(lineSizer,0,WACV)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Resolution: '),0,WACV)
            mapRes =  wx.TextCtrl(General,value='%.2f'%(Map['Resolution']),style=wx.TE_PROCESS_ENTER)
            mapRes.Bind(wx.EVT_TEXT_ENTER,OnResVal)        
            mapRes.Bind(wx.EVT_KILL_FOCUS,OnResVal)
            line2Sizer.Add(mapRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' Peak cutoff %: '),0,WACV)
            cutOff =  wx.TextCtrl(General,value='%.1f'%(Map['cutOff']),style=wx.TE_PROCESS_ENTER)
            cutOff.Bind(wx.EVT_TEXT_ENTER,OnCutOff)        
            cutOff.Bind(wx.EVT_KILL_FOCUS,OnCutOff)
            line2Sizer.Add(cutOff,0,WACV)
            mapSizer.Add(line2Sizer,0,WACV)
            return mapSizer
                
        def FlipSizer():
            if 'k-Max' not in Flip: Flip['k-Max'] = 20.
            
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
            
            def OnkMax(event):
                try:
                    res = float(kMax.GetValue())
                    if res >= 10.:
                        Flip['k-Max'] = res
                except ValueError:
                    pass
                kMax.SetValue("%.1f"%(Flip['k-Max']))          #reset in case of error

            refList = data['Histograms'].keys()
            flipSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Charge flip controls: Reflection set from: '),0,WACV)
            refList = wx.ComboBox(General,-1,value=Flip['RefList'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' Normalizing element: '),0,WACV)
            normElem = wx.Button(General,label=Flip['Norm element'],style=wx.TE_READONLY)
            normElem.Bind(wx.EVT_BUTTON,OnNormElem)
            lineSizer.Add(normElem,0,WACV)
            flipSizer.Add(lineSizer,0,WACV)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' Resolution: '),0,WACV)
            flipRes =  wx.TextCtrl(General,value='%.2f'%(Flip['Resolution']),style=wx.TE_PROCESS_ENTER)
            flipRes.Bind(wx.EVT_TEXT_ENTER,OnResVal)        
            flipRes.Bind(wx.EVT_KILL_FOCUS,OnResVal)
            line2Sizer.Add(flipRes,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Factor (0.1-1.2): '),0,WACV)
            kFactor =  wx.TextCtrl(General,value='%.3f'%(Flip['k-factor']),style=wx.TE_PROCESS_ENTER)
            kFactor.Bind(wx.EVT_TEXT_ENTER,OnkFactor)        
            kFactor.Bind(wx.EVT_KILL_FOCUS,OnkFactor)
            line2Sizer.Add(kFactor,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label=' k-Max (>=10.0): '),0,WACV)
            kMax = wx.TextCtrl(General,value='%.1f'%(Flip['k-Max']),style=wx.TE_PROCESS_ENTER)
            kMax.Bind(wx.EVT_TEXT_ENTER,OnkMax)        
            kMax.Bind(wx.EVT_KILL_FOCUS,OnkMax)
            line2Sizer.Add(kMax,0,WACV)
            flipSizer.Add(line2Sizer,0,WACV)
            return flipSizer
            
        def MCSASizer():
            Ind = {}
            
            def OnRefList(event):
                MCSAdata['Data source'] = refList.GetValue()
            
            def OnDmin(event):
                try:
                    val = float(dmin.GetValue())
                    if 1.0 <= val < 5.0:
                        MCSAdata['dmin'] = val
                except ValueError:
                    pass
                dmin.SetValue("%.3f"%(MCSAdata['dmin']))          #reset in case of error
                MCSAdata['newDmin'] = True

            def OnCycles(event):
                MCSAdata['Cycles'] = int(cycles.GetValue())
                               
            def OnAlist(event):
                MCSAdata['Algorithm'] = Alist.GetValue()
                wx.CallAfter(UpdateGeneral)
                
            def OnSlope(event):
                try:
                    val = float(slope.GetValue())
                    if .25 <= val < 1.0:
                        MCSAdata['log slope'] = val
                except ValueError:
                    pass
                slope.SetValue("%.3f"%(MCSAdata['log slope']))          #reset in case of error                
            
            def OnAjump(event):
                Obj = event.GetEventObject()
                name,ind = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    if .0 <= val <= 1.0:
                        MCSAdata[name][ind] = val
                except ValueError:
                    pass
                Obj.SetValue("%.3f"%(MCSAdata[name][ind]))
                
            def OnRanStart(event):
                MCSAdata['ranStart'] = ranStart.GetValue()
                
            def OnAutoRan(event):
                MCSAdata['autoRan'] = autoRan.GetValue()
                
            def OnRanRange(event):
                try:
                    val = float(ranRange.GetValue())/100
                    if 0.01 <= val <= 0.99:
                        MCSAdata['ranRange'] = val
                except ValueError:
                    pass
                ranRange.SetValue('%.1f'%(MCSAdata['ranRange']*100.))
            
            def OnAnneal(event):
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
                       
            refList = []
            if len(data['Pawley ref']):
                refList = ['Pawley reflections']
            for item in data['Histograms'].keys():
                if 'HKLF' in item or 'PWDR' in item:
                    refList.append(item)
            mcsaSizer = wx.BoxSizer(wx.VERTICAL)
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(General,label=' Monte Carlo/Simulated Annealing controls: Reflection set from: '),0,WACV)
            refList = wx.ComboBox(General,-1,value=MCSAdata['Data source'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            lineSizer.Add(refList,0,WACV)
            lineSizer.Add(wx.StaticText(General,label=' d-min: '),0,WACV)
            dmin = wx.TextCtrl(General,-1,value='%.3f'%(MCSAdata['dmin']),style=wx.TE_PROCESS_ENTER)
            dmin.Bind(wx.EVT_TEXT_ENTER,OnDmin)        
            dmin.Bind(wx.EVT_KILL_FOCUS,OnDmin)
            lineSizer.Add(dmin,0,WACV)
            mcsaSizer.Add(lineSizer)
            mcsaSizer.Add((5,5),)
            line2Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line2Sizer.Add(wx.StaticText(General,label=' MC/SA runs: '),0,WACV)
            Cchoice = ['1','2','4','8','16','32','64','128','256']
            cycles = wx.ComboBox(General,-1,value=str(MCSAdata.get('Cycles',1)),choices=Cchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            cycles.Bind(wx.EVT_COMBOBOX,OnCycles)        
            line2Sizer.Add(cycles,0,WACV)
            line2Sizer.Add((5,0),)
            ranStart = wx.CheckBox(General,-1,label=' MC/SA Refine at ')
            ranStart.Bind(wx.EVT_CHECKBOX, OnRanStart)
            ranStart.SetValue(MCSAdata.get('ranStart',False))
            line2Sizer.Add(ranStart,0,WACV)
            ranRange = wx.TextCtrl(General,-1,value='%.1f'%(MCSAdata.get('ranRange',0.10)*100),style=wx.TE_PROCESS_ENTER)
            ranRange.Bind(wx.EVT_TEXT_ENTER,OnRanRange)        
            ranRange.Bind(wx.EVT_KILL_FOCUS,OnRanRange)
            line2Sizer.Add(ranRange,0,WACV)
            line2Sizer.Add(wx.StaticText(General,label='% of ranges. '),0,WACV)
#            autoRan = wx.CheckBox(General,-1,label=' Do auto range reduction? ')
#            autoRan.Bind(wx.EVT_CHECKBOX, OnAutoRan)
#            autoRan.SetValue(MCSAdata.get('autoRan',False))
#            line2Sizer.Add(autoRan,0,WACV)
            mcsaSizer.Add(line2Sizer)
            mcsaSizer.Add((5,5),)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            Achoice = ['log','fast']                #these work
#            Achoice = ['log','fast','cauchy','boltzmann']
            line3Sizer.Add(wx.StaticText(General,label=' MC/SA schedule: '),0,WACV)
            Alist = wx.ComboBox(General,-1,value=MCSAdata['Algorithm'],choices=Achoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Alist.Bind(wx.EVT_COMBOBOX,OnAlist)
            line3Sizer.Add(Alist,0,WACV)
            if MCSAdata['Algorithm'] in ['fast','boltzmann','cauchy']:
                Names = [' A-jump: ',' B-jump: ']
                parms = 'Jump coeff'
                if MCSAdata['Algorithm'] in ['boltzmann','cauchy']:
                    Names = [' A-jump: ']
                elif 'fast' in MCSAdata['Algorithm']:
                    Names = [' quench: ',' m-factor: ',' n-factor: ']
                    parms = 'fast parms'
                for i,name in enumerate(Names):
                    line3Sizer.Add(wx.StaticText(General,label=name),0,WACV)
                    Ajump =  wx.TextCtrl(General,-1,value='%.3f'%(MCSAdata[parms][i]),style=wx.TE_PROCESS_ENTER)
                    Ajump.Bind(wx.EVT_TEXT_ENTER,OnAjump)        
                    Ajump.Bind(wx.EVT_KILL_FOCUS,OnAjump)
                    Indx[Ajump.GetId()] = [parms,i]
                    line3Sizer.Add(Ajump,0,WACV)
            elif 'log' in MCSAdata['Algorithm']:
                line3Sizer.Add(wx.StaticText(General,label=' slope: '),0,WACV)
                slope =  wx.TextCtrl(General,-1,value='%.3f'%(MCSAdata['log slope']),style=wx.TE_PROCESS_ENTER)
                slope.Bind(wx.EVT_TEXT_ENTER,OnSlope)        
                slope.Bind(wx.EVT_KILL_FOCUS,OnSlope)
                line3Sizer.Add(slope,0,WACV)
            mcsaSizer.Add(line3Sizer)
            mcsaSizer.Add((5,5),)
            line3Sizer = wx.BoxSizer(wx.HORIZONTAL)
            line3Sizer.Add(wx.StaticText(General,label=' Annealing schedule: '),0,WACV)
            names = [' Start temp: ',' Final temp: ',' No. trials: ']
            fmts = ['%.1f','%.5f','%d']
            for i,[name,fmt] in enumerate(zip(names,fmts)):
                if MCSAdata['Annealing'][i]:
                    text = fmt%(MCSAdata['Annealing'][i])
                else:
                    text = 'None'
                line3Sizer.Add(wx.StaticText(General,label=name),0,WACV)
                anneal =  wx.TextCtrl(General,-1,value=text,style=wx.TE_PROCESS_ENTER)
                anneal.Bind(wx.EVT_TEXT_ENTER,OnAnneal)        
                anneal.Bind(wx.EVT_KILL_FOCUS,OnAnneal)
                Indx[anneal.GetId()] = [i,fmt]
                line3Sizer.Add(anneal,0,WACV)
            mcsaSizer.Add(line3Sizer)            
            return mcsaSizer

        # UpdateGeneral execution continues here
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
        G2gd.HorizontalLine(mainSizer,General)
        
        if generalData['Type'] in ['modulated','magnetic',]:
            mainSizer.Add(ModulatedSizer(generalData['Type']))
            G2gd.HorizontalLine(mainSizer,General)

        mainSizer.Add(PawleySizer())
        G2gd.HorizontalLine(mainSizer,General)
        
        mainSizer.Add(MapSizer())
        G2gd.HorizontalLine(mainSizer,General)
        
        mainSizer.Add(FlipSizer())
        if generalData['Type'] in ['nuclear','macromolecular']:
            G2gd.HorizontalLine(mainSizer,General)
            mainSizer.Add(MCSASizer())
        SetPhaseWindow(G2frame.dataFrame,General,mainSizer)
        G2frame.dataFrame.SetStatusText('')

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
            if r < 0:                          #double click on col label! Change all atoms!
                sel = -1
                noSkip = True
                if Atoms.GetColLabelValue(c) == 'refine':
                    Type = generalData['Type']
                    if Type in ['nuclear','macromolecular','modulated']:
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
                if sel >= 0 and noSkip:
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
                            DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                    wx.CallAfter(Paint)
                    
        def ChangeAtomCell(event):

            def chkUij(Uij,CSI): #needs to do something!!!
                return Uij

            r,c =  event.GetRow(),event.GetCol()
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
                elif Atoms.GetColLabelValue(c) == 'refine':
                    ci = colLabels.index('I/A')
                    atomData[r][c] = atomData[r][c].replace(rbAtmDict.get(atomData[r][ci+8],''),'')
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                wx.CallAfter(Paint)

        def AtomTypeSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if Atoms.GetColLabelValue(c) == 'Type':
                PE = G2elemGUI.PickElement(G2frame)
                if PE.ShowModal() == wx.ID_OK:
                    if PE.Elem != 'None':                        
                        atomData[r][c] = PE.Elem.strip()
                        name = atomData[r][c]
                        if len(name) in [2,4]:
                            atomData[r][c-1] = name[:2]+'(%d)'%(r+1)
                        else:
                            atomData[r][c-1] = name[:1]+'(%d)'%(r+1)
                PE.Destroy()
                SetupGeneral()
                wx.CallAfter(Paint)
                value = Atoms.GetCellValue(r,c)
                atomData[r][c] = value
                ci = colLabels.index('I/A')
                ID = atomData[r][ci+8]
                if 'Atoms' in data['Drawing']:
                    DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                SetupGeneral()
            else:
                event.Skip()

        def RowSelect(event):
            r,c =  event.GetRow(),event.GetCol()
            if not (event.AltDown() or (event.ShiftDown() and event.ControlDown())):
                Atoms.frm = -1
                G2frame.dataFrame.SetStatusText('')                    
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
                            G2frame.dataFrame.SetStatusText('Atom '+atomData[r][n]+' is to be moved')
                        else:                       #move it
                            item = atomData.pop(Atoms.frm)
                            atomData.insert(r,item)
                            Atoms.frm = -1
                            G2frame.dataFrame.SetStatusText('')
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
            atomTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Atoms.SetTable(atomTable, True)
            Atoms.frm = -1            
            colType = colLabels.index('Type')
            colR = colLabels.index('refine')
            colSS = colLabels.index('site sym')
            colX = colLabels.index('x')
            colIA = colLabels.index('I/A')
            colU11 = colLabels.index('U11')
            colUiso = colLabels.index('Uiso')
            attr = wx.grid.GridCellAttr()
            attr.IncRef()               #fix from Jim Hester
            attr.SetEditor(G2gd.GridFractionEditor(Atoms))
            for c in range(colX,colX+3):
                attr = wx.grid.GridCellAttr()
                attr.IncRef()               #fix from Jim Hester
                attr.SetEditor(G2gd.GridFractionEditor(Atoms))
                Atoms.SetColAttr(c, attr)
            for i in range(colU11-1,colU11+6):
                Atoms.SetColSize(i,50)            
            for row in range(Atoms.GetNumberRows()):
                atId = atomData[row][colIA+8]
                rbExcl = rbAtmDict.get(atId,'')
                Atoms.SetReadOnly(row,colSS,True)                         #site sym
                Atoms.SetReadOnly(row,colSS+1,True)                       #Mult
                if Atoms.GetCellValue(row,colIA) == 'A':
                    try:    #patch for sytsym name changes
                        CSI = G2spc.GetCSuinel(atomData[row][colSS])
                    except KeyError:
                        Sytsym = G2spc.SytSym(atomData[row][colX:colX+3],SGData)[0]
                        atomData[row][colSS] = Sytsym
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
                if 'X' in rbExcl:
                    for c in range(0,colX+3):
                        if c != colR:
                            Atoms.SetCellStyle(row,c,VERY_LIGHT_GREY,True)
            Atoms.AutoSizeColumns(False)

        # FillAtomsGrid executable code starts here
        generalData = data['General']
        atomData = data['Atoms']
        DData = data['Drawing']
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
        Items = [G2gd.wxID_ATOMSEDITINSERT,G2gd.wxID_ATOMSEDITDELETE,G2gd.wxID_ATOMSREFINE, 
            G2gd.wxID_ATOMSMODIFY,G2gd.wxID_ATOMSTRANSFORM,G2gd.wxID_ATOMVIEWINSERT,G2gd.wxID_ATOMMOVE]
        if atomData:
            for item in Items:    
                G2frame.dataFrame.AtomsMenu.Enable(item,True)
        else:
            for item in Items:
                G2frame.dataFrame.AtomsMenu.Enable(item,False)
        Items = [G2gd.wxID_ATOMVIEWINSERT, G2gd.wxID_ATOMSVIEWADD,G2gd.wxID_ATOMMOVE]
        if 'showABC' in data['Drawing']:
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
        if generalData['Type'] == 'macromolecular':
            colLabels = ['res no','residue','chain'] + colLabels
            Types = [wg.GRID_VALUE_STRING,
                wg.GRID_VALUE_CHOICE+AAchoice,
                wg.GRID_VALUE_STRING] + Types
        SGData = data['General']['SGData']
        G2frame.dataFrame.SetStatusText('')
        if SGData['SGPolax']:
            G2frame.dataFrame.SetStatusText('Warning: The location of the origin is arbitrary in '+SGData['SGPolax'])
        Atoms.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeAtomCell)
        Atoms.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, AtomTypeSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, RefreshAtomGrid)
        Atoms.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
        Atoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, ChangeSelection)
        Atoms.SetMargins(0,0)
        
        G2frame.dataFrame.setSizePosLeft([700,300])
        Paint()

    def OnAtomAdd(event):
        AtomAdd(0,0,0)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        if data['Drawing']:
            G2plt.PlotStructure(G2frame,data)
        
    def OnAtomViewAdd(event):
        try:
            drawData = data['Drawing']
            x,y,z = drawData['viewPoint'][0]
            AtomAdd(x,y,z)
        except:
            AtomAdd(0,0,0)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        G2plt.PlotStructure(G2frame,data)
                
    def AtomAdd(x,y,z,El='H',Name='UNK'):
        atomData = data['Atoms']
        generalData = data['General']
        Ncol = Atoms.GetNumberCols()
        atId = ran.randint(0,sys.maxint)
        E,SGData = G2spc.SpcGroup(generalData['SGData']['SpGrp'])
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)
        if generalData['Type'] == 'macromolecular':
            atomData.append([0,Name,'',Name,El,'',x,y,z,1,Sytsym,Mult,'I',0.10,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'nuclear':
            atomData.append([Name,El,'',x,y,z,1,Sytsym,Mult,'I',0.01,0,0,0,0,0,0,atId])
        elif generalData['Type'] in ['modulated','magnetic']:
            atomData.append([Name,El,'',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,atId,[],[],[],[]])
        SetupGeneral()
        if 'Atoms' in data['Drawing']:            
            DrawAtomAdd(data['Drawing'],atomData[-1])

    def OnAtomInsert(event):
        AtomInsert(0,0,0)
        FillAtomsGrid(Atoms)
        event.StopPropagation()
        G2plt.PlotStructure(G2frame,data)
        
    def OnAtomViewInsert(event):
        if 'Drawing' in data:
            drawData = data['Drawing']
            x,y,z = drawData['viewPoint'][0]
            AtomAdd(x,y,z)
            FillAtomsGrid(Atoms)
        event.StopPropagation()
        
    def OnAtomMove(event):
        drawData = data['Drawing']
        atomData = data['Atoms']
        x,y,z = drawData['viewPoint'][0]
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        cx = colLabels.index('x')
        ci = colLabels.index('I/A')
        indx = Atoms.GetSelectedRows()
        if len(indx) != 1:
            G2frame.ErrorDialog('Atom move error','Only one atom can be moved')
        elif atomData[indx[0]][ci+8] in rbAtmDict:
            G2frame.ErrorDialog('Atom move error','Atoms in rigid bodies can not be moved')
        else:
            atomData[indx[0]][cx:cx+3] = [x,y,z]
            SetupGeneral()
            FillAtomsGrid(Atoms)
            ID = atomData[indx[0]][ci+8]
            DrawAtomsReplaceByID(data['Drawing'],atomData[indx[0]],ID)
            G2plt.PlotStructure(G2frame,data)
        event.StopPropagation()
            
    def DrawAtomsReplaceByID(drawingData,atom,ID):
        IDs = [ID,]
        atomData = drawingData['Atoms']
        indx = G2mth.FindAtomIndexByIDs(atomData,IDs)
        for ind in indx:
            atomData[ind] = MakeDrawAtom(atom,atomData[ind])
                
    def MakeDrawAtom(atom,oldatom=None):
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        generalData = data['General']
        SGData = generalData['SGData']
        if generalData['Type'] in ['nuclear','modulated',]:
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
            elif generalData['Type'] in ['modulated','magnetic']:
                atomData.insert(indx,['UNK','UNK','',x,y,z,1,Sytsym,Mult,0,'I',0.01,0,0,0,0,0,0,atId,[],[],[],[]])
            SetupGeneral()

    def AtomDelete(event):
        colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
        ci = colLabels.index('I/A')
        indx = Atoms.GetSelectedRows()
        IDs = []
        if indx:
            atomData = data['Atoms']
            indx.reverse()
            for ind in indx:
                atom = atomData[ind]
                if atom[ci+8] in rbAtmDict:
                    G2frame.dataFrame.SetStatusText('**** ERROR - atom is in a rigid body and can not be deleted ****')
                else:
                    IDs.append(atom[ci+8])
                    del atomData[ind]
            if 'Atoms' in data['Drawing']:
                DrawAtomsDeleteByIDs(IDs)
                wx.CallAfter(FillAtomsGrid,Atoms)
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
                choice = ['F - site fraction','X - coordinates','U - thermal parameters']
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
        indx = Atoms.GetSelectedRows()
        if indx:
            atomData = data['Atoms']
            generalData = data['General']
            colLabels = [Atoms.GetColLabelValue(c) for c in range(Atoms.GetNumberCols())]
            ci = colLabels.index('I/A')
            choices = ['Type','Name','x','y','z','frac','I/A','Uiso']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Atom parameter',choices)
            parm = ''
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                parm = choices[sel]
                cid = colLabels.index(parm)
            dlg.Destroy()
            if parm in ['Type']:
                dlg = G2elemGUI.PickElement(G2frame)
                if dlg.ShowModal() == wx.ID_OK:
                    if dlg.Elem not in ['None']:
                        El = dlg.Elem.strip()
                        for r in indx:                        
                            if not Atoms.IsReadOnly(r,cid):
                                atomData[r][cid] = El
                                if len(El) in [2,4]:
                                    atomData[r][cid-1] = El[:2]+'(%d)'%(r+1)
                                else:
                                    atomData[r][cid-1] = El[:1]+'(%d)'%(r+1)
                        SetupGeneral()
                        if 'Atoms' in data['Drawing']:
                            for r in indx:
                                ID = atomData[r][ci+8]
                                DrawAtomsReplaceByID(data['Drawing'],atomData[r],ID)
                    FillAtomsGrid(Atoms)
                dlg.Destroy()
            elif parm in ['Name',]:
                dlg = wx.MessageDialog(G2frame,'Do you really want to rename the selected atoms?','Rename', 
                    wx.YES_NO | wx.ICON_QUESTION)
                try:
                    result = dlg.ShowModal()
                    if result == wx.ID_YES:
                        for r in indx:
                            if not Atoms.IsReadOnly(r,cid+1):
                                El = atomData[r][cid+1]
                                if len(El) in [2,4]:
                                    atomData[r][cid] = El[:2]+'(%d)'%(r+1)
                                else:
                                    atomData[r][cid] = El[:1]+'(%d)'%(r+1)
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
                        if not Atoms.IsReadOnly(r,cid):
                            atomData[r][cid] = parm
                    FillAtomsGrid(Atoms)
                dlg.Destroy()
            elif parm in ['frac','Uiso']:
                limits = [0.,1.]
                val = 1.0
                if  parm in ['Uiso']:
                    limits = [0.,0.25]
                    val = 0.01
                dlg = G2gd.SingleFloatDialog(G2frame,'New value','Enter new value for '+parm,val,limits)
                if dlg.ShowModal() == wx.ID_OK:
                    parm = dlg.GetValue()
                    for r in indx:                        
                        if not Atoms.IsReadOnly(r,cid):
                            atomData[r][cid] = parm
                    SetupGeneral()
                    FillAtomsGrid(Atoms)
                dlg.Destroy()
            elif parm in ['x','y','z']:
                limits = [-1.,1.]
                val = 0.
                dlg = G2gd.SingleFloatDialog(G2frame,'Atom shift','Enter shift for '+parm,val,limits)
                if dlg.ShowModal() == wx.ID_OK:
                    parm = dlg.GetValue()
                    for r in indx:                        
                        if not Atoms.IsReadOnly(r,cid):
                            atomData[r][cid] += parm
                    SetupGeneral()
                    FillAtomsGrid(Atoms)
                dlg.Destroy()

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
            dlg = G2gd.SymOpDialog(G2frame,SGData,True,True)
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
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        if Force:
                            XYZ = G2spc.MoveToUnitCell(XYZ)
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
                FillAtomsGrid(Atoms)
            else:
                Atoms.ForceRefresh()

    def OnDistAnglePrt(event):
        'save distances and angles to a file'    
        fp = file(os.path.abspath(os.path.splitext(G2frame.GSASprojectfile
                                                   )[0]+'.disagl'),'w')
        OnDistAngle(event,fp=fp)
        fp.close()
    
    def OnDistAngle(event,fp=None):
        'Compute distances and angles'    
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
            dlg = G2gd.DisAglDialog(G2frame,DisAglCtls,generalData)
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
                DisAglData['covData'] = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Covariance'))
            try:
                if fp:
                    G2stMn.PrintDistAngle(DisAglCtls,DisAglData,fp)
                else:    
                    G2stMn.PrintDistAngle(DisAglCtls,DisAglData)
            except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
                G2frame.ErrorDialog('Distance/Angle calculation','try again but do "Reset" to fill in missing atom types')
        else:
            print "select one or more rows of atoms"
            G2frame.ErrorDialog('Select atom',"select one or more rows of atoms then redo")
                        
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
            raise Exception,"Should not happen: 'ISODISTORT' not in data"
        if len(data.get('Histograms',[])) == 0:
            G2frame.ErrorDialog(
                'No data',
                'Sorry, this computation requires that a histogram first be added to the phase'
                )
            return
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree() # init for constraint
        # make a lookup table for constraints
        sub = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Constraints') 
        Constraints = G2frame.PatternTree.GetItemPyData(sub)
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
                    ch = G2gd.HelpButton(panel2,fmtHelp(constDict[G2var],var))
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
                    ch = G2gd.HelpButton(panel2,fmtHelp(constDict[G2var],var))
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
        reqrdr = G2frame.dataFrame.ReImportMenuId.get(event.GetId())
        rdlist = G2frame.OnImportGeneric(reqrdr,
            G2frame.ImportPhaseReaderlist,'phase')
        if len(rdlist) == 0: return
        # rdlist is only expected to have one element
        rd = rdlist[0]
        G2frame.OnFileSave(event)
        # rd contains all info for a phase
        PhaseName = rd.Phase['General']['Name']
        print 'Read phase '+str(PhaseName)+' from file '+str(G2frame.lastimport)
        atomData = data['Atoms']
        atomNames = []
        for atom in atomData:
            atomNames.append(atom[:ct+1])
        for atom in rd.Phase['Atoms']:
            try:
                idx = atomNames.index(atom[:ct+1])
                atId = atom[cia+8]
                atomData[idx][:-1] = atom[:-1]
                atomData[idx][cia+8] = atId
            except ValueError:
                print atom[:ct+1], 'not in Atom array; not updated'
        wx.CallAfter(FillAtomsGrid,Atoms)
        
################################################################################
#### Wave Data page
################################################################################

    def UpdateWavesData():
        
        def AtomSizer(SS,atom):
            
            def OnWaveType(event):
                atom[-1][SS]['waveType']=waveType.GetValue()
                
            def OnShowWave(event):
                Obj = event.GetEventObject()
                atom = Indx[Obj.GetId()]               
                Ax = Obj.GetValue()
                G2plt.ModulationPlot(G2frame,data,atom,Ax)
                
            atomSizer = wx.BoxSizer(wx.HORIZONTAL)
            atomSizer.Add(wx.StaticText(waveData,label=' Modulation data for atom:    '+atom[0]+'    WaveType: '),0,WACV)            
            waveType = wx.ComboBox(waveData,value=atom[-1][SS]['waveType'],choices=waveTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            waveType.Bind(wx.EVT_COMBOBOX,OnWaveType)
            atomSizer.Add(waveType,0,WACV)
            axchoice = ['x','y','z']
            atomSizer.Add(wx.StaticText(waveData,label=' Show contour map for axis:'),0,WACV)
            mapSel = wx.ComboBox(waveData,value=' ',choices=axchoice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            mapSel.Bind(wx.EVT_COMBOBOX,OnShowWave)
            Indx[mapSel.GetId()] = atom
            atomSizer.Add(mapSel,0,WACV)
            return atomSizer
            
        def WaveSizer(waveBlk,Stype,typeName,Names):
            
            def OnAddWave(event):
                Obj = event.GetEventObject()
                iatm,item = Indx[Obj.GetId()]
                atomData[iatm][-1][SS][item].append([[0.0 for i in range(numVals[Stype])],False])
                UpdateWavesData()
                
            def OnWaveVal(event):
                Obj = event.GetEventObject()
                iatm,item,iwave,ival = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                except ValueError:
                    val = atomData[iatm][-1][SS][item][iwave][0][ival]
                Obj.SetValue('%.4f'%val)
                atomData[iatm][-1][SS][item][iwave][0][ival] = val
                
            def OnRefWave(event):
                Obj = event.GetEventObject()
                iatm,item,iwave = Indx[Obj.GetId()]
                atomData[iatm][-1][SS][item][iwave][1] = not atomData[iatm][-1][SS][item][iwave][1]
                
            def OnDelWave(event):
                Obj = event.GetEventObject()
                iatm,item,iwave = Indx[Obj.GetId()]
                del atomData[iatm][-1][SS][item][iwave]
                UpdateWavesData()                
                
            waveSizer = wx.BoxSizer(wx.VERTICAL)
            waveHead = wx.BoxSizer(wx.HORIZONTAL)
            waveHead.Add(wx.StaticText(waveData,label=typeName+' modulation parameters: '),0,WACV)
            waveAdd = wx.CheckBox(waveData,label='Add wave?')
            waveAdd.Bind(wx.EVT_CHECKBOX, OnAddWave)
            Indx[waveAdd.GetId()] = [iatm,Stype]
            waveHead.Add(waveAdd,0,WACV)
            waveSizer.Add(waveHead)
            if len(waveBlk):
                waveSizer.Add(wx.StaticText(waveData,label=' Parameters: '+str(Names).rstrip(']').lstrip('[').replace("'",'')),0,WACV)
                if Stype == 'Sfrac':
                    Waves = wx.FlexGridSizer(1,4,5,5)
                else:
                    Waves = wx.FlexGridSizer(1,8,5,5)
                for iwave,wave in enumerate(waveBlk):
                    for ival,val in enumerate(wave[0]):
                        waveVal = wx.TextCtrl(waveData,value='%.4f'%(val),style=wx.TE_PROCESS_ENTER)
                        waveVal.Bind(wx.EVT_TEXT_ENTER,OnWaveVal)
                        waveVal.Bind(wx.EVT_KILL_FOCUS,OnWaveVal)
                        Indx[waveVal.GetId()] = [iatm,Stype,iwave,ival]
                        Waves.Add(waveVal,0,WACV)
                        if len(wave[0]) > 6 and ival == 5:
                            Waves.Add((5,5),0)
                            Waves.Add((5,5),0)
                    waveRef = wx.CheckBox(waveData,label='Refine?')
                    waveRef.SetValue(wave[1])
                    Indx[waveRef.GetId()] = [iatm,Stype,iwave]
                    waveRef.Bind(wx.EVT_CHECKBOX, OnRefWave)
                    Waves.Add(waveRef,0,WACV)
                    waveDel = wx.CheckBox(waveData,label='Delete?')
                    Indx[waveDel.GetId()] = [iatm,Stype,iwave]
                    waveDel.Bind(wx.EVT_CHECKBOX, OnDelWave)
                    Waves.Add(waveDel,0,WACV)                
                waveSizer.Add(Waves)                    
            return waveSizer
            
        def MapSizer():

            def OnRefList(event):
                Map['RefList'] = refList.GetValue()
                
            def OnMapType(event):
                Map['MapType'] = mapType.GetValue()
                
            Map = generalData['4DmapData']
            Map['Resolution'] = 0.25
            refList = data['Histograms'].keys()
            mapSizer = wx.BoxSizer(wx.HORIZONTAL)
            mapSizer.Add(wx.StaticText(waveData,label=' 4D map data: Reflection set from: '),0,WACV)
            refList = wx.ComboBox(waveData,-1,value=Map['RefList'],choices=refList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            refList.Bind(wx.EVT_COMBOBOX,OnRefList)
            mapSizer.Add(refList,0,WACV)
            mapTypes = ['Fobs','delt-F']
            mapSizer.Add(wx.StaticText(waveData,label=' Map type: '),0,WACV)
            mapType = wx.ComboBox(waveData,-1,value=Map['MapType'],choices=mapTypes,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            mapType.Bind(wx.EVT_COMBOBOX,OnMapType)
            mapSizer.Add(mapType,0,WACV)
            return mapSizer
            
        Indx = {}
        G2frame.dataFrame.SetStatusText('')
        generalData = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        atomData = data['Atoms']
        if waveData.GetSizer():
            waveData.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        typeNames = {'Sfrac':' Site fraction','Spos':' Position','Sadp':' Thermal motion','Smag':' Magnetic moment'}
        numVals = {'Sfrac':2,'Spos':6,'Sadp':12,'Smag':6}
        posNames = ['Xsin','Ysin','Zsin','Xcos','Ycos','Zcos']
        adpNames = ['U11sin','U22sin','U33sin','U12sin','U13sin','U23sin',
            'U11cos','U22cos','U33cos','U12cos','U13cos','U23cos']
        magNames = ['MXsin','MYsin','MZsin','MXcos','MYcos','MZcos']
        fracNames = ['Flen','Fcent','Fsin','Fcos']
        waveTypes = ['Fourier','Sawtooth','ZigZag','Crenel/Fourier']
        Labels = {'Spos':posNames,'Sfrac':fracNames,'Sadp':adpNames,'Smag':magNames}
        mainSizer.Add(wx.StaticText(waveData,label=' Incommensurate propagation wave data:'),0,WACV)
        if generalData['Type'] in ['modulated','magnetic']:
            mainSizer.Add(MapSizer(),0,WACV)            
            for iatm,atom in enumerate(atomData):
                for SS in ['SS1',]:  #future SS2 & SS3 - I doubt it!
                    G2gd.HorizontalLine(mainSizer,waveData)
                    mainSizer.Add(AtomSizer(SS,atom))
                    for Stype in ['Sfrac','Spos','Sadp','Smag']:
                        if generalData['Type'] == 'modulated' and Stype == 'Smag':
                            break
                        mainSizer.Add(WaveSizer(atom[-1][SS][Stype],Stype,typeNames[Stype],Labels[Stype]))
                        
        SetPhaseWindow(G2frame.dataFrame,waveData,mainSizer)
                       
    def On4DMapCompute(event):
        generalData = data['General']
        mapData = generalData['4DmapData']
        reflName = mapData['RefList']
        if not reflName:
            G2frame.ErrorDialog('Fourier map','No reflections defined for Fourier map')
            return
        phaseName = generalData['Name']
        if 'PWDR' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflData = reflSets[phaseName]
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflData = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        mapData.update(G2mth.Fourier4DMap(data,reflData))
        mapSig = np.std(mapData['rho'])
        print mapData['MapType']+' computed: rhomax = %.3f rhomin = %.3f sigma = %.3f'%(np.max(mapData['rho']),np.min(mapData['rho']),mapSig)
            
################################################################################
#Structure drawing GUI stuff                
################################################################################

    def SetupDrawingData():
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        atomData = data['Atoms']
        AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH','WAT','UNK']
        AA1letter = ['A','R','N','D','C','Q','E','G','H','I',
            'L','K','M','F','P','S','T','W','Y','V','M',' ',' ',' ']
        defaultDrawing = {'Atoms':[],'viewPoint':[[0.5,0.5,0.5],[]],'showHydrogen':True,
            'backColor':[0,0,0],'depthFog':False,'Zclip':50.0,'cameraPos':50.,'Zstep':0.5,
            'radiusFactor':0.85,'contourLevel':1.,'bondRadius':0.1,'ballScale':0.33,
            'vdwScale':0.67,'ellipseProb':50,'sizeH':0.50,'unitCellBox':True,
            'showABC':True,'selectedAtoms':[],'Atoms':[],'oldxy':[],
            'bondList':{},'viewDir':[1,0,0]}
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
            drawingData.update(defaultDrawing)
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
        cx,ct,cs,ci = [0,0,0,0]
        if generalData['Type'] in ['nuclear','modulated']:
            cx,ct,cs,ci = [2,1,6,17]         #x, type, style & index
        elif generalData['Type'] == 'macromolecular':
            cx,ct,cs,ci = [5,4,9,20]         #x, type, style & index
        elif generalData['Type'] == 'magnetic':
            cx,ct,cs,ci = [2,1,6,20]         #x, type, style & index
#        elif generalData['Type'] == 'modulated':
#           ?????   for future
        drawingData['atomPtrs'] = [cx,ct,cs,ci]
        if not drawingData.get('Atoms'):
            for atom in atomData:
                DrawAtomAdd(drawingData,atom)
            data['Drawing'] = drawingData
            
    def DrawAtomAdd(drawingData,atom):
        drawingData['Atoms'].append(MakeDrawAtom(atom))
        
    def OnRestraint(event):        
        indx = drawAtoms.GetSelectedRows()
        restData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Restraints'))
        drawingData = data['Drawing']
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])            
        cx,ct,cs,ci = drawingData['atomPtrs']
        atomData = drawingData['Atoms']
        atNames = []
        atXYZ = []
        atSymOp = []
        atIndx = []
        for item in indx:
            atXYZ.append(np.array(atomData[item][cx:cx+3]))
            atSymOp.append(atomData[item][cs-1])
            atIndx.append(atomData[item][ci])
        if event.GetId() == G2gd.wxID_DRAWRESTRBOND and len(indx) == 2:
            try:
                bondData = restData[PhaseName]['Bond']
            except KeyError:
                bondData = {'wtFactor':1.0,'Bonds':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Bond'] = bondData
            dist = G2mth.getRestDist(atXYZ,Amat)
            bondData['Bonds'].append([atIndx,atSymOp,1.54,0.01])
        elif event.GetId() == G2gd.wxID_DRAWRESTRANGLE and len(indx) == 3:
            try:
                angleData = restData[PhaseName]['Angle']
            except KeyError:
                angleData = {'wtFactor':1.0,'Angles':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Angle'] = angleData
            angle = G2mth.getRestAngle(atXYZ,Amat)
            angleData['Angles'].append([atIndx,atSymOp,109.5,1.0])            
        elif event.GetId() == G2gd.wxID_DRAWRESTRPLANE and len(indx) > 3:
            try:
                planeData = restData[PhaseName]['Plane']
            except KeyError:
                planeData = {'wtFactor':1.0,'Planes':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Plane'] = planeData
            plane = G2mth.getRestPlane(atXYZ,Amat)
            planeData['Planes'].append([atIndx,atSymOp,0.0,0.01])            
        elif event.GetId() == G2gd.wxID_DRAWRESTRCHIRAL and len(indx) == 4:
            try:
                chiralData = restData[PhaseName]['Chiral']
            except KeyError:
                chiralData = {'wtFactor':1.0,'Volumes':[],'Use':True}
                restData[PhaseName] = {}
                restData[PhaseName]['Chiral'] = chiralData
            volume = G2mth.getRestChiral(atXYZ,Amat)
            chiralData['Volumes'].append([atIndx,atSymOp,2.5,0.1])            
        else:
            print '**** ERROR wrong number of atoms selected for this restraint'
            return
        G2frame.PatternTree.SetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Restraints'),restData)

    def OnDefineRB(event):
        indx = drawAtoms.GetSelectedRows()
        indx.sort()
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
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
        rbId = ran.randint(0,sys.maxint)
        rbName = 'UNKRB'
        dlg = wx.TextEntryDialog(G2frame,'Enter the name for the new rigid body',
            'Edit rigid body name',rbName ,style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            rbName = dlg.GetValue()
        dlg.Destroy()
        RBData['Residue'][rbId] = {'RBname':rbName,'rbXYZ':rbXYZ,'rbTypes':rbType,
            'atNames':atNames,'rbRef':[0,1,2,False],'rbSeq':[],'SelSeq':[0,0],'useCount':0}
        RBData['RBIds']['Residue'].append(rbId)
        G2frame.dataFrame.SetStatusText('New rigid body UNKRB added to set of Residue rigid bodies')

################################################################################
##### Atom draw routines
################################################################################
            
    def UpdateDrawAtoms(atomStyle=''):
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
        G2frame.dataFrame.SetStatusText('')
        generalData = data['General']
        SetupDrawingData()
        drawingData = data['Drawing']
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
#        elif generalData['Type'] == 'modulated':
#            Types += []
#            colLabels += []
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
            FindBondsDraw()
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
                G2frame.dataFrame.SetStatusText('Select Custom Color, change color, Add to Custom Colors, then OK')
            else:
                G2frame.dataFrame.SetStatusText('Change color, Add to Custom Colors, then OK')
            generalData = data['General']
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
            G2frame.dataFrame.SetStatusText('')
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
            dlg = G2gd.SymOpDialog(G2frame,SGData,False,True)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Inv,Cent,Opr,Cell,New,Force = dlg.GetSelection()
                    Cell = np.array(Cell)
                    cent = SGData['SGCen'][Cent]
                    M,T = SGData['SGOps'][Opr]
                    for ind in indx:
                        XYZ = np.array(atomData[ind][cx:cx+3])
                        XYZ = np.inner(M,XYZ)+T
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        if Force:
                            XYZ = G2spc.MoveToUnitCell(XYZ)
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
            dlg = G2gd.SymOpDialog(G2frame,SGData,False,True)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    Inv,Cent,Opr,Cell,New,Force = dlg.GetSelection()
                    Cell = np.array(Cell)
                    cent = SGData['SGCen'][Cent]
                    M,T = SGData['SGOps'][Opr]
                    for ind in indx:
                        XYZ = np.array(atomData[ind][cx:cx+3])
                        XYZ = np.inner(M,XYZ)+T
                        if Inv:
                            XYZ = -XYZ
                        XYZ = XYZ+cent+Cell
                        if Force:
                            XYZ = G2spc.MoveToUnitCell(XYZ)
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
            cx,ct,cs,ci = data['Drawing']['atomPtrs']
            generalData = data['General']
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
                                    oprB = atomB[cx+3]
                                    C = xyz-xyzB
                                    newOp = '1+'+str(int(round(C[0])))+','+str(int(round(C[1])))+','+str(int(round(C[2])))
                                    newAtom = atomB[:]
                                    newAtom[cx:cx+3] = xyz
                                    newAtom[cx+3] = G2spc.StringOpsProd(oprB,newOp,SGData)
                                    atomData.append(newAtom)
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
            cx = colLabels.index('x')
            cuia = colLabels.index('I/A')
            cuij = cuia+2
            generalData = data['General']
            SGData = generalData['SGData']
            wx.BeginBusyCursor()
            try:
                for ind in indx:
                    atom = atomData[ind]
                    XYZ = np.array(atom[cx:cx+3])
                    if atom[cuia] == 'A':
                        Uij = atom[cuij:cuij+6]
                        result = G2spc.GenAtom(XYZ,SGData,False,Uij,True)
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
                        result = G2spc.GenAtom(XYZ,SGData,False,Move=True)
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
            finally:
                wx.EndBusyCursor()
            UpdateDrawAtoms()
            drawAtoms.ClearSelection()
            G2plt.PlotStructure(G2frame,data)
            
    def FindBondsToo():                         #works but slow for large structures - keep as reference
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
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
        cx,ct,cs,ci = data['Drawing']['atomPtrs']
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
        
    def DrawAtomsDeleteByIDs(IDs):
        atomData = data['Drawing']['Atoms']
        indx = G2mth.FindAtomIndexByIDs(atomData,IDs)
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
        indx = G2mth.FindAtomIndexByIDs(atomData,IDs)
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
        G2stMn.BestPlane(PlaneData)
        
    def OnDrawDistVP(event):
        # distance to view point
        indx = drawAtoms.GetSelectedRows()
        if not indx:
            print '***** ERROR - no atoms selected'
            return
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])            
        drawingData = data['Drawing']
        viewPt = np.array(drawingData['viewPoint'][0])
        print ' Distance from view point at %.3f %.3f %.3f to:'%(viewPt[0],viewPt[1],viewPt[2])
        atomDData = drawingData['Atoms']
        colLabels = [drawAtoms.GetColLabelValue(c) for c in range(drawAtoms.GetNumberCols())]
        cx = colLabels.index('x')
        cn = colLabels.index('Name')
        for i in indx:
            atom = atomDData[i]
            Dx = np.array(atom[cx:cx+3])-viewPt
            dist = np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0))
            print 'Atom: %8s (%12s) distance = %.3f'%(atom[cn],atom[cx+3],dist)
    
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
            id = G2mth.FindAtomIndexByIDs(atomData,[atom[cid],],False)[0]
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
        G2stMn.DisAglTor(DATData)
                        
################################################################################
#### Draw Options page
################################################################################

    def UpdateDrawOptions():
        import copy
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
                
            def OnZstep(event):
                try:
                    step = float(Zstep.GetValue())
                    if not (0.01 <= step <= 1.0):
                        raise ValueError
                except ValueError:
                    step = drawingData['Zstep']
                drawingData['Zstep'] = step
                Zstep.SetValue('%.2fA'%(drawingData['Zstep']))
                
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
                
            def OnContourLevel(event):
                drawingData['contourLevel'] = contourLevel.GetValue()/100.
                contourLevelTxt.SetLabel(' Contour level: '+'%.2f'%(drawingData['contourLevel']*generalData['Map']['rhoMax']))
                G2plt.PlotStructure(G2frame,data)

            def OnMapSize(event):
                drawingData['mapSize'] = mapSize.GetValue()/10.
                mapSizeTxt.SetLabel(' Map radius, A: '+'%.1f'%(drawingData['mapSize']))
                G2plt.PlotStructure(G2frame,data)

            
            slopSizer = wx.BoxSizer(wx.HORIZONTAL)
            slideSizer = wx.FlexGridSizer(0,2)
            slideSizer.AddGrowableCol(1,1)
    
            cameraPosTxt = wx.StaticText(drawOptions,-1,
                ' Camera Distance: '+'%.2f'%(drawingData['cameraPos']),name='cameraPos')
            G2frame.dataDisplay.cameraPosTxt = cameraPosTxt
            slideSizer.Add(cameraPosTxt,0,WACV)
            cameraPos = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=drawingData['cameraPos'],name='cameraSlider')
            cameraPos.SetRange(10,500)
            cameraPos.Bind(wx.EVT_SLIDER, OnCameraPos)
            G2frame.dataDisplay.cameraSlider = cameraPos
            slideSizer.Add(cameraPos,1,wx.EXPAND|wx.RIGHT)
            
            ZclipTxt = wx.StaticText(drawOptions,-1,' Z clipping: '+'%.2fA'%(drawingData['Zclip']*drawingData['cameraPos']/100.))
            slideSizer.Add(ZclipTxt,0,WACV)
            Zclip = wx.Slider(drawOptions,style=wx.SL_HORIZONTAL,value=drawingData['Zclip'])
            Zclip.SetRange(1,99)
            Zclip.Bind(wx.EVT_SLIDER, OnZclip)
            slideSizer.Add(Zclip,1,wx.EXPAND|wx.RIGHT)
            
            ZstepSizer = wx.BoxSizer(wx.HORIZONTAL)
            ZstepSizer.Add(wx.StaticText(drawOptions,-1,' Z step:'),0,WACV)
            Zstep = wx.TextCtrl(drawOptions,value='%.2f'%(drawingData['Zstep']),
                style=wx.TE_PROCESS_ENTER)
            Zstep.Bind(wx.EVT_TEXT_ENTER,OnZstep)
            Zstep.Bind(wx.EVT_KILL_FOCUS,OnZstep)
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
                
            def OnShowRB(event):
                drawingData['showRigidBodies'] = showRB.GetValue()
                FindBondsDraw()
                G2plt.PlotStructure(G2frame,data)
                
            def OnViewPoint(event):
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
            G2frame.dataDisplay.viewDir = viewDir
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
            G2frame.dataDisplay.viewPoint = viewPoint
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
            radSizer.Add(wx.StaticText(drawOptions,-1,' Hydrogen radius, A:  '),0,WACV)
            sizeH = wx.TextCtrl(drawOptions,-1,value='%.2f'%(drawingData['sizeH']),size=wx.Size(60,20),style=wx.TE_PROCESS_ENTER)
            sizeH.Bind(wx.EVT_TEXT_ENTER,OnSizeHatoms)
            sizeH.Bind(wx.EVT_KILL_FOCUS,OnSizeHatoms)
            radSizer.Add(sizeH,0,WACV)
    
            radSizer.Add(wx.StaticText(drawOptions,-1,' Bond search factor:  '),0,WACV)
            radFactor = wx.TextCtrl(drawOptions,value='%.2f'%(drawingData['radiusFactor']),size=wx.Size(60,20),style=wx.TE_PROCESS_ENTER)
            radFactor.Bind(wx.EVT_TEXT_ENTER,OnRadFactor)
            radFactor.Bind(wx.EVT_KILL_FOCUS,OnRadFactor)
            radSizer.Add(radFactor,0,WACV)
            return radSizer

        # UpdateDrawOptions exectable code starts here
        generalData = data['General']
        Amat,Bmat = G2lat.cell2AB(generalData['Cell'][1:7])
        SetupDrawingData()
        drawingData = data['Drawing']
        if generalData['Type'] == 'nuclear':
            pickChoice = ['Atoms','Bonds','Torsions','Planes']
        elif generalData['Type'] == 'macromolecular':
            pickChoice = ['Atoms','Residues','Chains','Bonds','Torsions','Planes','phi/psi']

        G2frame.dataFrame.SetStatusText('')
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
        SetPhaseWindow(G2frame.dataFrame,drawOptions,mainSizer)

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

        # UpdateTexture executable starts here
        G2frame.dataFrame.SetStatusText('')
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
        shAngles = ['omega','chi','phi']
        if Texture.GetSizer():
            Texture.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        titleSizer = wx.BoxSizer(wx.HORIZONTAL)
        titleSizer.Add(wx.StaticText(Texture,-1,'Spherical harmonics texture data for '+PhaseName+':'),0,WACV)
        titleSizer.Add(wx.StaticText(Texture,-1,
            ' Texture Index J = %7.3f'%(G2lat.textureIndex(textureData['SH Coeff'][1]))),
            0,WACV)
        mainSizer.Add(titleSizer,0)
        mainSizer.Add((0,5),0)
        shSizer = wx.FlexGridSizer(0,6,5,5)
        shSizer.Add(wx.StaticText(Texture,-1,'Texture model: '),0,WACV)
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
        PTSizer = wx.FlexGridSizer(0,4,5,5)
        PTSizer.Add(wx.StaticText(Texture,-1,' Texture plot type: '),0,WACV)
        choices = ['Axial pole distribution','Pole figure','Inverse pole figure']            
        pfType = wx.ComboBox(Texture,-1,value=str(textureData['PlotType']),choices=choices,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        pfType.Bind(wx.EVT_COMBOBOX,OnPfType)
        PTSizer.Add(pfType,0,WACV)
        if 'Axial' not in textureData['PlotType']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Projection type: '),0,WACV)
            projSel = wx.ComboBox(Texture,-1,value=G2frame.Projection,choices=['equal area','stereographic','3D display'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            projSel.Bind(wx.EVT_COMBOBOX,OnProjSel)
            PTSizer.Add(projSel,0,WACV)
        if textureData['PlotType'] in ['Pole figure','Axial pole distribution']:
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
        if 'Axial' not in textureData['PlotType']:
            PTSizer.Add(wx.StaticText(Texture,-1,' Color scheme'),0,WACV)
            choice = [m for m in mpl.cm.datad.keys() if not m.endswith("_r")]
            choice.sort()
            colorSel = wx.ComboBox(Texture,-1,value=G2frame.ContourColor,choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            colorSel.Bind(wx.EVT_COMBOBOX,OnColorSel)
            PTSizer.Add(colorSel,0,WACV)        
        mainSizer.Add(PTSizer,0,WACV)
        mainSizer.Add((0,5),0)
        if textureData['SHShow']:
            mainSizer.Add(wx.StaticText(Texture,-1,'Spherical harmonic coefficients: '),0,WACV)
            mainSizer.Add((0,5),0)
            ODFSizer = wx.FlexGridSizer(0,8,2,2)
            ODFIndx = {}
            ODFkeys = textureData['SH Coeff'][1].keys()
            ODFkeys.sort()
            for item in ODFkeys:
                ODFSizer.Add(wx.StaticText(Texture,-1,item),0,WACV)
                ODFval = wx.TextCtrl(Texture,wx.ID_ANY,'%8.3f'%(textureData['SH Coeff'][1][item]),style=wx.TE_PROCESS_ENTER)
                ODFIndx[ODFval.GetId()] = item
                ODFval.Bind(wx.EVT_TEXT_ENTER,OnODFValue)
                ODFval.Bind(wx.EVT_KILL_FOCUS,OnODFValue)
                ODFSizer.Add(ODFval,0,WACV)
            mainSizer.Add(ODFSizer,0,WACV)
            mainSizer.Add((0,5),0)
        mainSizer.Add((0,5),0)
        mainSizer.Add(wx.StaticText(Texture,-1,'Sample orientation angles: '),0,WACV)
        mainSizer.Add((0,5),0)
        angSizer = wx.BoxSizer(wx.HORIZONTAL)
        angIndx = {}
        valIndx = {}
        for item in ['Sample omega','Sample chi','Sample phi']:
            angRef = wx.CheckBox(Texture,-1,label=item+': ')
            angRef.SetValue(textureData[item][0])
            angIndx[angRef.GetId()] = item
            angRef.Bind(wx.EVT_CHECKBOX, OnAngRef)
            angSizer.Add(angRef,0,WACV)
            angVal = wx.TextCtrl(Texture,wx.ID_ANY,'%8.2f'%(textureData[item][1]),style=wx.TE_PROCESS_ENTER)
            valIndx[angVal.GetId()] = item
            angVal.Bind(wx.EVT_TEXT_ENTER,OnAngValue)
            angVal.Bind(wx.EVT_KILL_FOCUS,OnAngValue)
            angSizer.Add(angVal,0,WACV)
            angSizer.Add((5,0),0)
        mainSizer.Add(angSizer,0,WACV)
        SetPhaseWindow(G2frame.dataFrame,Texture,mainSizer)

################################################################################
##### DData routines - GUI stuff in GSASIIddataGUI.py
################################################################################
        
    def OnHklfAdd(event):
        UseList = data['Histograms']
        keyList = UseList.keys()
        TextList = []
        if not G2frame.PatternTree.GetCount():
            return
        
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
            UseList[histoName] = {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                                  'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
                                  'Extinction':['Lorentzian','None',
                                                {'Tbar':0.1,'Cos2TM':0.955,'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]},]}                        
            UpdateHKLFdata(histoName)
            data['Histograms'] = UseList
        wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
        wx.EndBusyCursor()
                
    def UpdateHKLFdata(histoName):
        generalData = data['General']
        Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,histoName)
        refDict,reflData = G2frame.PatternTree.GetItemPyData(Id)
        SGData = generalData['SGData']
        Cell = generalData['Cell'][1:7]
        G,g = G2lat.cell2Gmat(Cell)
        for iref,ref in enumerate(reflData['RefList']):
            H = list(ref[:3])
            ref[4] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
            iabsnt,ref[3],Uniq,phi = G2spc.GenHKLf(H,SGData)
        #G2frame.PatternTree.SetItemPyData(Id,[refDict,reflData]) #removed by BHT -- not needed!
        
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
                        Id = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,histoName)
                        UseList[histoName] = {'Histogram':histoName,'Show':False,
                            'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{}],
                            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                                [1.,1.,1.,0.,0.,0.],6*[False,]],
                            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                            'HStrain':[NDij*[0.0,],NDij*[False,]],                          
                            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}
                        refList = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,Id,'Reflection Lists'))
                        refList[generalData['Name']] = []                       
                    data['Histograms'] = UseList
                    wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
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
                    wx.CallAfter(G2ddG.UpdateDData,G2frame,DData,data)
            finally:
                dlg.Destroy()
                
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
            
            def OnThermval(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    RBObj['ThermalMotion'][1][item] = val
                except ValueError:
                    pass
                Obj.SetValue('%8.4f'%(RBObj['ThermalMotion'][1][item]))
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
                thermVal = wx.TextCtrl(RigidBodies,-1,value='%8.4f'%(model[1][i]),
                    style=wx.TE_PROCESS_ENTER)
                thermVal.Bind(wx.EVT_TEXT_ENTER,OnThermval)
                thermVal.Bind(wx.EVT_KILL_FOCUS,OnThermval)
                Indx[thermVal.GetId()] = i
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
                
            def OnOrigX(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    RBObj['Orig'][0][item] = val
                    Obj.SetValue('%8.5f'%(val))
                    newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,rbType)[0]
                    for i,id in enumerate(RBObj['Ids']):
                        data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
                    data['Drawing']['Atoms'] = []
                    UpdateDrawAtoms(atomStyle)
                    G2plt.PlotStructure(G2frame,data)
                except ValueError:
                    pass
                
            def OnOrien(event):
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
                origX = wx.TextCtrl(RigidBodies,-1,value='%8.5f'%(x),style=wx.TE_PROCESS_ENTER)
                origX.Bind(wx.EVT_TEXT_ENTER,OnOrigX)
                origX.Bind(wx.EVT_KILL_FOCUS,OnOrigX)
                Indx[origX.GetId()] = ix
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
            G2frame.dataFrame.SetStatusText('NB: Rotation vector is in crystallographic space')
             
            def OnTorsionRef(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                RBObj['Torsions'][item][1] = Obj.GetValue()                
                
            def OnTorsion(event):
                Obj = event.GetEventObject()
                item = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    RBObj['Torsions'][item][0] = val
                    newXYZ = G2mth.UpdateRBXYZ(Bmat,RBObj,RBData,'Residue')[0]
                    for i,id in enumerate(RBObj['Ids']):
                        data['Atoms'][AtLookUp[id]][cx:cx+3] = newXYZ[i]
                except ValueError:
                    pass
                Obj.SetValue("%10.3f"%(RBObj['Torsions'][item][0]))                
                data['Drawing']['Atoms'] = []
                UpdateDrawAtoms(atomStyle)
                drawAtoms.ClearSelection()
                G2plt.PlotStructure(G2frame,data)
                
            def OnDelResRB(event):
                Obj = event.GetEventObject()
                RBId = Indx[Obj.GetId()]
                RBData['Residue'][RBId]['useCount'] -= 1
                RBObjs = data['RBModels']['Residue']
                for rbObj in RBObjs:
                    if RBId == rbObj['RBId']:
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
                torsTxt = wx.TextCtrl(RigidBodies,-1,value='%.3f'%(tors[0]),style=wx.TE_PROCESS_ENTER)
                torsTxt.Bind(wx.EVT_TEXT_ENTER,OnTorsion)
                torsTxt.Bind(wx.EVT_KILL_FOCUS,OnTorsion)
                Indx[torsTxt.GetId()] = itors
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
            G2frame.dataFrame.SetStatusText('NB: Rotation vector is in crystallographic space')
                   
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
        
        # FillRigidBodyGrid executable code starts here
        if refresh:
            RigidBodies.DestroyChildren()
        general = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        AtLookUp = G2mth.FillAtomLookUp(data['Atoms'],cia+8)
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        Indx = {}
        atomStyle = 'balls & sticks'
        if 'macro' in general['Type']:
            atomStyle = 'sticks'
        G2frame.dataFrame.SetStatusText('')
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if not data['RBModels']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'No rigid body models:'),0,WACV)
            mainSizer.Add((5,5),0)
        if 'Residue' in data['RBModels']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'Residue rigid bodies:'),0,WACV)
            mainSizer.Add((5,5),0)
            for RBObj in data['RBModels']['Residue']:
                mainSizer.Add(ResrbSizer(RBObj))
        if 'Vector' in data['RBModels']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(RigidBodies,-1,'Vector rigid bodies:'),0,WACV)
            mainSizer.Add((5,5),0)
            for RBObj in data['RBModels']['Vector']:
                mainSizer.Add(VecrbSizer(RBObj))

        SetPhaseWindow(G2frame.dataFrame,RigidBodies,mainSizer)

    def OnRBCopyParms(event):
        RBObjs = []
        for rbType in ['Vector','Residue']:            
            RBObjs += data['RBModels'].get(rbType,[])
        if not len(RBObjs):
            print '**** ERROR - no rigid bodies defined ****'
            return
        if len(RBObjs) == 1:
            print '**** INFO - only one rigid body defined; nothing to copy to ****'
            return
        Source = []
        sourceRB = {}
        for RBObj in RBObjs:
            Source.append(RBObj['RBname'])
        dlg = wx.SingleChoiceDialog(G2frame,'Select source','Copy rigid body parameters',Source)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            name = Source[sel]
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
        wx.CallAfter(FillRigidBodyGrid(True))
                
    def OnRBAssign(event):
        
        G2frame.dataFrame.SetStatusText('')
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbNames = {}
        for rbVec in RBData['Vector']:
            if rbVec != 'AtInfo':
                rbNames[RBData['Vector'][rbVec]['RBname']] =['Vector',rbVec]
        for rbRes in RBData['Residue']:
            if rbRes != 'AtInfo':
                rbNames[RBData['Residue'][rbRes]['RBname']] = ['Residue',rbRes]
        if not rbNames:
            print '**** ERROR - no rigid bodies defined ****'
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
                    Ids.append(atomData[id][-1])
                    atomData[id][cx:cx+3] = xyz
                if dmax > 1.0:
                    print '**** WARNING - some atoms not found or misidentified ****'
                    print '****           check torsion angles & try again      ****'
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
                    'Orient':[[0.,0.,0.,1.],' '],'Ids':[],'RBId':rbId,'Torsions':[],
                    'numChain':'','RBname':RBData[rbType][rbId]['RBname']}
                data['testRBObj']['torAtms'] = []                
                for item in RBData[rbType][rbId].get('rbSeq',[]):
                    data['testRBObj']['rbObj']['Torsions'].append([item[2],False])
                    data['testRBObj']['torAtms'].append([-1,-1,-1])
                Draw()
                
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
                VAC = np.inner(Amat,np.array(atomData[atInd[1]][cx:cx+3])-Orig)
                VBC = np.inner(Amat,np.array(atomData[atInd[2]][cx:cx+3])-Orig)
                VCC = np.cross(VAR,VAC)
                QuatA = G2mth.makeQuat(VAR,VAC,VCC)[0]
                VAR = G2mth.prodQVQ(QuatA,VAR)
                VBR = G2mth.prodQVQ(QuatA,VBR)
                QuatB = G2mth.makeQuat(VBR,VBC,VAR)[0]
                QuatC = G2mth.prodQQ(QuatB,QuatA)
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
                
            def OnTorAngle(event):
                OkBtn.SetLabel('OK')
                OkBtn.Enable(True)
                Obj = event.GetEventObject()
                [tor,torSlide] = Indx[Obj.GetId()]
                Tors = data['testRBObj']['rbObj']['Torsions'][tor]
                try:
                    value = float(Obj.GetValue())
                except ValueError:
                    value = Tors[0]
                Tors[0] = value
                Obj.SetValue('%8.3f'%(value))
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
                ang.SetValue('%8.3f'%(val))
                G2plt.PlotStructure(G2frame,data)

            if len(data['testRBObj']):
                G2plt.PlotStructure(G2frame,data)
                    
            RigidBodies.DestroyChildren()
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
                    choice = atNames[0].keys()
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
                        choice = ['',]+atNames[i].keys()
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
                    AtNames = data['testRBObj']['AtNames']
                    rbAtTypes = data['testRBObj']['rbAtTypes']
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
                        ang = wx.TextCtrl(RigidBodies,-1,value='%8.3f'%(torsion[0]),style=wx.TE_PROCESS_ENTER)
                        ang.Bind(wx.EVT_TEXT_ENTER,OnTorAngle)
                        ang.Bind(wx.EVT_KILL_FOCUS,OnTorAngle)
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
                rbSel = wx.ComboBox(RigidBodies,-1,value='',choices=rbNames.keys(),
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
            SetPhaseWindow(G2frame.dataFrame,RigidBodies,mainSizer)
        Draw()
        
    def OnAutoFindResRB(event):
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbKeys = RBData['Residue'].keys()
        rbKeys.remove('AtInfo')
        if not len(rbKeys):
            print '**** ERROR - no residue rigid bodies are defined ****'
            return
        RBNames = [RBData['Residue'][k]['RBname'] for k in rbKeys]
        RBIds = dict(zip(RBNames,rbKeys))
        general = data['General']
        cx,ct,cs,cia = generalData['AtomPtrs']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        Atoms = data['Atoms']
        AtLookUp = G2mth.FillAtomLookUp(Atoms,cia+8)
        if 'macro' not in general['Type']:
            print '**** ERROR - this phase is not a macromolecule ****'
            return
        if not len(Atoms):
            print '**** ERROR - this phase has no atoms ****'
            return
        RBObjs = []
        cx,ct = general['AtomPtrs'][:2]
        iatm = 0
        wx.BeginBusyCursor()
        try:
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
                rbObj = {'RBname':rbRes['RBname']+':'+str(rbRes['useCount']),'numChain':numChain}
                rbAtoms = []
                rbIds = []
                for iratm in range(len(rbRes['atNames'])):
                    rbAtoms.append(np.array(Atoms[iatm][cx:cx+3]))
                    rbIds.append(Atoms[iatm][20])
                    iatm += 1    #puts this at beginning of next residue?
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
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
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
                    elif 'T' in rbObj['ThermalMotion'][0]:
                        if 'Tii' in parms:
                            rbObj['ThermalMotion'][2][0:2] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][0:2] = [False,False,False]
                        if 'Tij' in parms:
                            rbObj['ThermalMotion'][2][3:6] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][3:6] = [False,False,False]
                    elif 'L' in rbObj['ThermalMotion'][0]:
                        if 'Lii' in parms:
                            rbObj['ThermalMotion'][2][6:9] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][6:9] = [False,False,False]
                        if 'Lij' in parms:
                            rbObj['ThermalMotion'][2][9:12] = [True,True,True]
                        else:
                            rbObj['ThermalMotion'][2][9:12] = [False,False,False]
                    elif 'S' in rbObj['ThermalMotion'][0]:
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
            
        def OnPosVal(event):
            Obj = event.GetEventObject()
            model,item,ix = Indx[Obj.GetId()]
            try:
                model[item][0][ix] = float(Obj.GetValue())
            except ValueError:
                pass
            Obj.SetValue("%.4f"%(model[item][0][ix]))
            G2plt.PlotStructure(G2frame,data)
            
        def OnPosRange(event):
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
            atomsizer.Add(wx.StaticText(MCSA,-1,' Atom: '+model['name']+': '),0,WACV)
            for ix,item in enumerate(['x','y','z']):
                posRef = wx.CheckBox(MCSA,-1,label=item+': ')
                posRef.SetValue(model['Pos'][1][ix])
                posRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                Indx[posRef.GetId()] = [model,'Pos',ix]
                atomsizer.Add(posRef,0,WACV)
                posVal = wx.TextCtrl(MCSA,-1,'%.4f'%(model['Pos'][0][ix]),style=wx.TE_PROCESS_ENTER)
                posVal.Bind(wx.EVT_TEXT_ENTER,OnPosVal)
                posVal.Bind(wx.EVT_KILL_FOCUS,OnPosVal)
                Indx[posVal.GetId()] = [model,'Pos',ix]
                atomsizer.Add(posVal,0,WACV)
            atomsizer.Add((5,5),0)
            for ix,item in enumerate(['x','y','z']):
                atomsizer.Add(wx.StaticText(MCSA,-1,' Range: '),0,WACV)
                rmin,rmax = model['Pos'][2][ix]
                posRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
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
#                UpdateMCSA()

            def OnMolCent(event):
                Obj = event.GetEventObject()
                model = Indx[Obj.GetId()]
                model['MolCent'][1] = Obj.GetValue()
                if model['MolCent'][1]:
                    G2mth.SetMolCent(model,RBData)                
                G2plt.PlotStructure(G2frame,data)
            
            rbsizer = wx.BoxSizer(wx.VERTICAL)
            rbsizer1 = wx.FlexGridSizer(0,7,5,5)
            rbsizer1.Add(wx.StaticText(MCSA,-1,model['Type']+': '+model['name']+': '),0,WACV)
            for ix,item in enumerate(['x','y','z']):
                posRef = wx.CheckBox(MCSA,-1,label=item+': ')
                posRef.SetValue(model['Pos'][1][ix])
                posRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                Indx[posRef.GetId()] = [model,'Pos',ix]
                rbsizer1.Add(posRef,0,WACV)
                posVal = wx.TextCtrl(MCSA,-1,'%.4f'%(model['Pos'][0][ix]),style=wx.TE_PROCESS_ENTER)
                posVal.Bind(wx.EVT_TEXT_ENTER,OnPosVal)
                posVal.Bind(wx.EVT_KILL_FOCUS,OnPosVal)
                Indx[posVal.GetId()] = [model,'Pos',ix]
                rbsizer1.Add(posVal,0,WACV)
            molcent = wx.CheckBox(MCSA,-1,label=' Use mol. center? ')
            molcent.SetValue(model['MolCent'][1])
            molcent.Bind(wx.EVT_CHECKBOX,OnMolCent)
            Indx[molcent.GetId()] = model
            rbsizer1.Add(molcent,0,WACV)
            for ix,item in enumerate(['x','y','z']):
                rbsizer1.Add(wx.StaticText(MCSA,-1,' Range: '),0,WACV)
                rmin,rmax = model['Pos'][2][ix]
                posRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[posRange.GetId()] = [model,'Pos',ix]
                posRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                posRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                rbsizer1.Add(posRange,0,WACV)
                
            rbsizer2 = wx.FlexGridSizer(0,6,5,5)
            Ori = model['Ori'][0]
            rbsizer2.Add(wx.StaticText(MCSA,-1,'Oa: '),0,WACV)
            angVal = wx.TextCtrl(MCSA,-1,'%.5f'%(Ori[0]),style=wx.TE_PROCESS_ENTER)
            angVal.Bind(wx.EVT_TEXT_ENTER,OnOriVal)
            angVal.Bind(wx.EVT_KILL_FOCUS,OnOriVal)
            rbsizer2.Add(angVal,0,WACV)
            rbsizer2.Add(wx.StaticText(MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
            vecVal = wx.TextCtrl(MCSA,-1,'%.3f %.3f %.3f'%(Ori[1],Ori[2],Ori[3]),style=wx.TE_PROCESS_ENTER)
            vecVal.Bind(wx.EVT_TEXT_ENTER,OnOriVal)
            vecVal.Bind(wx.EVT_KILL_FOCUS,OnOriVal)
            Indx[angVal.GetId()] = [model,0,angVal,vecVal]
            Indx[vecVal.GetId()] = [model,1,angVal,vecVal]
            rbsizer2.Add(vecVal,0,WACV)
            rbsizer2.Add(wx.StaticText(MCSA,-1,' Vary? '),0,WACV)
            choice = [' ','A','AV']
            orvar = wx.ComboBox(MCSA,-1,value=model['Ovar'],choices=choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            orvar.Bind(wx.EVT_COMBOBOX, OnOrVar)
            Indx[orvar.GetId()] = model
            rbsizer2.Add(orvar,0,WACV)
            rbsizer2.Add(wx.StaticText(MCSA,-1,' Range: Oa: '),0,WACV)
            Rge = model['Ori'][2]
            angRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f'%(Rge[0][0],Rge[0][1]),style=wx.TE_PROCESS_ENTER)
            Indx[angRange.GetId()] = [model,'Ori',0]
            angRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
            angRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
            rbsizer2.Add(angRange,0,WACV)
            rbsizer2.Add(wx.StaticText(MCSA,-1,'Oi,Oj,Ok: '),0,WACV)
            for io,item in enumerate(['Oi','Oj','Ok']):
                rmin,rmax = Rge[io+1]
                vecRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f '%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                Indx[vecRange.GetId()] = [model,'Ori',io+1]
                vecRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                vecRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                rbsizer2.Add(vecRange,0,WACV)
            rbsizer.Add(rbsizer1)    
            rbsizer.Add(rbsizer2)    
            if model['Type'] == 'Residue':
                atNames = RBData['Residue'][model['RBId']]['atNames']
                rbsizer.Add(wx.StaticText(MCSA,-1,'Torsions:'),0,WACV)
                rbsizer3 = wx.FlexGridSizer(0,8,5,5)
                for it,tor in enumerate(model['Tor'][0]):
                    iBeg,iFin = RBData['Residue'][model['RBId']]['rbSeq'][it][:2]
                    name = atNames[iBeg]+'-'+atNames[iFin]
                    torRef = wx.CheckBox(MCSA,-1,label=' %s: '%(name))
                    torRef.SetValue(model['Tor'][1][it])
                    torRef.Bind(wx.EVT_CHECKBOX,OnPosRef)
                    Indx[torRef.GetId()] = [model,'Tor',it]
                    rbsizer3.Add(torRef,0,WACV)
                    torVal = wx.TextCtrl(MCSA,-1,'%.4f'%(tor),style=wx.TE_PROCESS_ENTER)
                    torVal.Bind(wx.EVT_TEXT_ENTER,OnPosVal)
                    torVal.Bind(wx.EVT_KILL_FOCUS,OnPosVal)
                    Indx[torVal.GetId()] = [model,'Tor',it]
                    rbsizer3.Add(torVal,0,WACV)
                    rbsizer3.Add(wx.StaticText(MCSA,-1,' Range: '),0,WACV)
                    rmin,rmax = model['Tor'][2][it]
                    torRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
                    Indx[torRange.GetId()] = [model,'Tor',it]
                    torRange.Bind(wx.EVT_TEXT_ENTER,OnPosRange)
                    torRange.Bind(wx.EVT_KILL_FOCUS,OnPosRange)
                    rbsizer3.Add(torRange,0,WACV)
                rbsizer.Add(rbsizer3)
                
            return rbsizer
            
        def MDSizer(POData):
            
            def OnPORef(event):
                POData['Coef'][1] = poRef.GetValue()
                
            def OnPOVal(event):
                try:
                    mdVal = float(poVal.GetValue())
                    if mdVal > 0:
                        POData['Coef'][0] = mdVal
                except ValueError:
                    pass
                poVal.SetValue("%.3f"%(POData['Coef'][0]))
                
            def OnPORange(event):
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
            poRef = wx.CheckBox(MCSA,-1,label=' March-Dollase ratio: ')
            poRef.SetValue(POData['Coef'][1])
            poRef.Bind(wx.EVT_CHECKBOX,OnPORef)
            poSizer.Add(poRef,0,WACV)
            poVal = wx.TextCtrl(MCSA,-1,'%.3f'%(POData['Coef'][0]),style=wx.TE_PROCESS_ENTER)
            poVal.Bind(wx.EVT_TEXT_ENTER,OnPOVal)
            poVal.Bind(wx.EVT_KILL_FOCUS,OnPOVal)
            poSizer.Add(poVal,0,WACV)
            poSizer.Add(wx.StaticText(MCSA,-1,' Range: '),0,WACV)
            rmin,rmax = POData['Coef'][2]
            poRange = wx.TextCtrl(MCSA,-1,'%.3f %.3f'%(rmin,rmax),style=wx.TE_PROCESS_ENTER)
            poRange.Bind(wx.EVT_TEXT_ENTER,OnPORange)
            poRange.Bind(wx.EVT_KILL_FOCUS,OnPORange)
            poSizer.Add(poRange,0,WACV)                       
            poSizer.Add(wx.StaticText(MCSA,-1,' Unique axis, H K L: '),0,WACV)
            h,k,l = POData['axis']
            poAxis = wx.TextCtrl(MCSA,-1,'%3d %3d %3d'%(h,k,l),style=wx.TE_PROCESS_ENTER)
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
                    G2plt.PlotStructure(G2frame,data)
                    wx.CallAfter(UpdateMCSA,MCSA.GetScrollPos(wx.VERTICAL))
                    resultsGrid.ForceRefresh()
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
            resultsTable = G2gd.Table(resultVals,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            resultsGrid = G2gd.GSGrid(MCSA)
            resultsGrid.SetTable(resultsTable, True)
            resultsGrid.Bind(wg.EVT_GRID_CELL_LEFT_CLICK, OnCellChange)
            resultsGrid.AutoSizeColumns(True)
            for r in range(resultsGrid.GetNumberRows()):
                for c in range(resultsGrid.GetNumberCols()):
                    if c in [0,1]:
                        resultsGrid.SetReadOnly(r,c,isReadOnly=False)
                    else:
                        resultsGrid.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            resultsSizer.Add(resultsGrid)
            return resultsSizer
        
        # UpdateMCSA executable code starts here
        MCSA.DestroyChildren()
        if not data['Drawing']:                 #if new drawing - no drawing data!
            SetupDrawingData()
        general = data['General']
        Amat,Bmat = G2lat.cell2AB(general['Cell'][1:7])
        RBData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        Indx = {}
        atomStyle = 'balls & sticks'
        if 'macro' in general['Type']:
            atomStyle = 'sticks'
        G2frame.dataFrame.SetStatusText('')
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if not data['MCSA']['Models']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(MCSA,-1,'No MC/SA models:'),0,WACV)
            mainSizer.Add((5,5),0)
        else:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(MCSA,-1,'MC/SA models:'),0,WACV)
            mainSizer.Add((5,5),0)
            for model in data['MCSA']['Models']:
                Xsize = 500
                if model['Type'] == 'MD':
                    mainSizer.Add(MDSizer(model))
                elif model['Type'] == 'Atom':
                    Asizer = atomSizer(model)
                    mainSizer.Add(Asizer)
                    Xsize = max(Asizer.GetMinSize()[0],Xsize)
                else:
                    Rsizer = rbSizer(model)
                    mainSizer.Add(Rsizer)
                    Xsize = max(Rsizer.GetMinSize()[0],Xsize)
                G2gd.HorizontalLine(mainSizer,MCSA)
                
        if not data['MCSA']['Results']:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(MCSA,-1,'No MC/SA results:'),0,WACV)
            mainSizer.Add((5,5),0)
        else:
            mainSizer.Add((5,5),0)
            mainSizer.Add(wx.StaticText(MCSA,-1,'MC/SA results:'),0,WACV)
            mainSizer.Add((5,5),0)
            Results = data['MCSA']['Results']
            mainSizer.Add(ResultsSizer(Results))
            
        SetPhaseWindow(G2frame.dataFrame,MCSA,mainSizer)
        Size = MCSA.GetSize()
        Size[0] = Xsize+40
        G2frame.dataFrame.SetSize(Size)
        MCSA.Scroll(0,Scroll)
        
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
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            try:        #patch for old reflection data
                reflData = reflSets[phaseName]['RefList']
            except TypeError:
                reflData = reflSets[phaseName]
            reflType = 'PWDR'
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            try:
                reflData = G2frame.PatternTree.GetItemPyData(PatternId)[1]['RefList']
            except TypeError:
                reflData = G2frame.PatternTree.GetItemPyData(PatternId)[1]
            reflType = 'HKLF'
        elif reflName == 'Pawley reflections':
            reflData = data['Pawley ref']
            covData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,G2frame.root, 'Covariance'))
            reflType = 'Pawley'
        else:
            print '**** ERROR - No data defined for MC/SA run'
            return
        print 'MC/SA run:'
        print 'Reflection type:',reflType,' Total No. reflections: ',len(reflData)
        RBdata = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        MCSAmodels = MCSAdata['Models']
        if not len(MCSAmodels):
            print '**** ERROR - no models defined for MC/SA run****'
            return
        time1 = time.time()
        if process == 'single':
            pgbar = wx.ProgressDialog('MC/SA','Residual Rcf =',101.0, 
                style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
            screenSize = wx.ClientDisplayRect()
            Size = pgbar.GetSize()
            Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
            pgbar.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
            pgbar.SetSize(Size)
        else:
            pgbar = None
        try:
            tsf = 0.
            nCyc = mcsaControls['Cycles']
            if process == 'single':
                for i in range(nCyc):
                    pgbar.SetTitle('MC/SA run '+str(i+1)+' of '+str(nCyc))
                    Result,tsum = G2mth.mcsaSearch(data,RBdata,reflType,reflData,covData,pgbar)
                    MCSAdata['Results'].append(Result)
                    print ' MC/SA run completed: %d residual: %.3f%% SFcalc time: %.2fs'%(i,100*Result[2],tsum)
                    tsf += tsum
                print ' Structure factor time: %.2f'%(tsf)
            else:
                MCSAdata['Results'] = G2mth.MPmcsaSearch(nCyc,data,RBdata,reflType,reflData,covData)
            print ' MC/SA run time: %.2f'%(time.time()-time1)
        finally:
            if process == 'single':
                pgbar.Destroy()
        MCSAdata['Results'] = G2mth.sortArray(MCSAdata['Results'],2,reverse=False)
        MCSAdata['Results'][0][0] = True
        SetSolution(MCSAdata['Results'][0],data['MCSA']['Models'])
        G2frame.dataDisplay.SetFocus()
        Page = G2frame.dataDisplay.FindPage('MC/SA')
        G2frame.dataDisplay.SetSelection(Page)
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
        rbData = G2frame.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbNames = {}
        for rbVec in rbData['Vector']:
            if rbVec != 'AtInfo':
                rbNames[rbData['Vector'][rbVec]['RBname']] = ['Vector',rbVec]
        for rbRes in rbData['Residue']:
            if rbRes != 'AtInfo':
                rbNames[rbData['Residue'][rbRes]['RBname']] = ['Residue',rbRes]
        if not rbNames:
            print '**** ERROR - no rigid bodies defined ****'
            return
        dlg = wx.SingleChoiceDialog(G2frame.dataFrame,'Select','Rigid body',rbNames.keys())
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            rbname = rbNames.keys()[sel]
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
        G2plt.PlotStructure(G2frame,data)
        UpdateMCSA()
        
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
        G2frame.dataFrame.SetStatusText('To delete a Pawley reflection: select row & press Delete')                        
        generalData = data['General']
        if 'Pawley ref' in data:
            PawleyPeaks = data['Pawley ref']                        
            rowLabels = []
            for i in range(len(PawleyPeaks)): rowLabels.append(str(i))
            if generalData['Type'] in ['modulated','magnetic',]:
                colLabels = ['h','k','l','m','mul','d','refine','Fsq(hkl)','sig(Fsq)']
                Types = 5*[wg.GRID_VALUE_LONG,]+[wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,]+ \
                    2*[wg.GRID_VALUE_FLOAT+':10,2',]
                pos = [6,7]
            else:    
                colLabels = ['h','k','l','mul','d','refine','Fsq(hkl)','sig(Fsq)']
                Types = 4*[wg.GRID_VALUE_LONG,]+[wg.GRID_VALUE_FLOAT+':10,4',wg.GRID_VALUE_BOOL,]+ \
                    2*[wg.GRID_VALUE_FLOAT+':10,2',]
                pos = [5,6]
            PawleyTable = G2gd.Table(PawleyPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
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
            G2frame.dataFrame.setSizePosLeft([450,300])
                    
    def OnPawleyLoad(event):
        generalData = data['General']
        cell = generalData['Cell'][1:7]
        A = G2lat.cell2A(cell)
        SGData = generalData['SGData']
        dmin = generalData['Pawley dmin']
        PawleyPeaks = []
        HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
        if generalData['Type'] in ['modulated','magnetic',]:
            Vec,x,maxH = generalData['SuperVec']
            SSGData = G2spc.SSpcGroup(SGData,generalData['SuperSg'])[1]
            wx.BeginBusyCursor()
            try:
                HKLd = G2lat.GenSSHLaue(dmin,SGData,SSGData,Vec,maxH,A)
                for h,k,l,m,d in HKLd:
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
        if generalData['Type'] in ['modulated','magnetic',]:
            im = 1
        HistoNames = filter(lambda a:Histograms[a]['Use']==True,Histograms.keys())
        PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,HistoNames[0])
        xdata = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        Inst = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Instrument Parameters'))[0]
        Sample = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Sample Parameters'))
        wave = G2mth.getWave(Inst)
        posCorr = Inst['Zero'][1]
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
                    FWHM = max(0.001,G2pwd.getFWHM(pos,Inst))/100.0
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
        Patterson squared, leBail extraction, etc.
        '''
        try:
            Refs = data['Pawley ref']
            Histograms = data['Histograms']
        except KeyError:
            G2frame.ErrorDialog('Pawley update','No histograms defined for this phase')
            return
        HistoNames = Histograms.keys()
        PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root,HistoNames[0])
        refData = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,  \
            PatternId,'Reflection Lists'))[PhaseName]['RefList']
        im = 0
        if data['General']['Type'] in ['modulated','magnetic',]:
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
                    print 'skipped',ref
                    pass
        finally:
            wx.EndBusyCursor()
        wx.CallAfter(FillPawleyReflectionsGrid)
                            
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
                if colLabels[c] == 'mag':
                    mapPeaks = G2mth.sortArray(mapPeaks,c,reverse=True)
                elif colLabels[c] in ['x','y','z','dzero']:
                    mapPeaks = G2mth.sortArray(mapPeaks,c)
                else:
                    return
                data['Map Peaks'] = mapPeaks
                wx.CallAfter(FillMapPeaksGrid)
            G2plt.PlotStructure(G2frame,data)                    
            
        G2frame.dataFrame.setSizePosLeft([450,300])
        G2frame.dataFrame.SetStatusText('')
        if 'Map Peaks' in data:
            G2frame.dataFrame.SetStatusText('Select mag or dzero columns to sort')
            mapPeaks = data['Map Peaks']                        
            rowLabels = []
            for i in range(len(mapPeaks)): rowLabels.append(str(i))
            colLabels = ['mag','x','y','z','dzero']
            Types = 5*[wg.GRID_VALUE_FLOAT+':10,4',]
            G2frame.MapPeaksTable = G2gd.Table(mapPeaks,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            MapPeaks.SetTable(G2frame.MapPeaksTable, True)
            MapPeaks.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
            for r in range(MapPeaks.GetNumberRows()):
                for c in range(MapPeaks.GetNumberCols()):
                    MapPeaks.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            MapPeaks.SetMargins(0,0)
            MapPeaks.AutoSizeColumns(False)
                    
    def OnPeaksMove(event):
        if 'Map Peaks' in data:
            mapPeaks = np.array(data['Map Peaks'])
            peakMax = np.max(mapPeaks.T[0])
            Ind = MapPeaks.GetSelectedRows()
            for ind in Ind:
                mag,x,y,z,d = mapPeaks[ind]
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
            mapPeaks = data['Map Peaks']
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
            G2frame.dataFrame.MapPeaksEdit.SetLabel(G2gd.wxID_SHOWBONDS,'Show bonds')
        else:
            generalData['Map']['Show bonds'] = True
            G2frame.dataFrame.MapPeaksEdit.SetLabel(G2gd.wxID_SHOWBONDS,'Hide bonds')
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
        print ' Distance from view point at %.3f %.3f %.3f to:'%(viewPt[0],viewPt[1],viewPt[2])
        colLabels = [MapPeaks.GetColLabelValue(c) for c in range(MapPeaks.GetNumberCols())]
        cx = colLabels.index('x')
        cm = colLabels.index('mag')
        for i in indx:
            peak = mapPeaks[i]
            Dx = np.array(peak[cx:cx+3])-viewPt
            dist = np.sqrt(np.sum(np.inner(Amat,Dx)**2,axis=0))
            print 'Peak: %5d mag= %8.2f distance = %.3f'%(i,peak[cm],dist)

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
            print ' distance for atoms %s = %.3f'%(str(indx),G2mth.getRestDist(xyz,Amat))
        else:
            print ' angle for atoms %s = %.2f'%(str(indx),G2mth.getRestAngle(xyz,Amat))
                                    
    def OnFourierMaps(event):
        generalData = data['General']
        mapData = generalData['Map']
        reflName = mapData['RefList']
        if not reflName:
            G2frame.ErrorDialog('Fourier map','No reflections defined for Fourier map')
            return
        phaseName = generalData['Name']
        if 'PWDR' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflData = reflSets[phaseName]
            if 'list' in str(type(reflData)):       #patch for old reflection data
            #if isinstance(reflData,list):       #patch for old reflection data
                RefData = {'RefList':[],'FF':[]}
                for ref in reflData:
                    RefData['RefList'].append(ref[:11]+[ref[13],])
                    RefData['FF'].append(ref[14])
                RefData['RefList'] = np.array(RefData['RefList'])
                reflData = RefData
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflData = G2frame.PatternTree.GetItemPyData(PatternId)[1]
        if 'Omit' in mapData['MapType']:
            pgbar = wx.ProgressDialog('Omit map','Blocks done',65, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
            mapData.update(G2mth.OmitMap(data,reflData,pgbar))
            pgbar.Destroy()
        else:
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
        
    def OnFourClear(event):
        generalData = data['General']
        generalData['Map'] = mapDefault
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
        generalData = data['General']
        mapData = generalData['Map']
        if len(mapData['rho']):
            wx.BeginBusyCursor()
            try:
                peaks,mags,dzeros = G2mth.SearchMap(data)
            finally:
                wx.EndBusyCursor()
            if len(peaks):
                mapPeaks = np.concatenate((mags,peaks,dzeros),axis=1)
                data['Map Peaks'] = G2mth.sortArray(mapPeaks,0,reverse=True)            
            print ' Map search finished, time = %.2fs'%(time.time()-time0)
            print ' No.peaks found:',len(peaks)    
            Page = G2frame.dataDisplay.FindPage('Map peaks')
            G2frame.dataDisplay.SetSelection(Page)
            wx.CallAfter(FillMapPeaksGrid)
            UpdateDrawAtoms()
        else:
            print 'No map available'
        
    def OnChargeFlip(event):
        generalData = data['General']
        mapData = generalData['Map']
        flipData = generalData['Flip']
        reflName = flipData['RefList']
        phaseName = generalData['Name']
        if 'PWDR' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflSets = G2frame.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(G2frame,PatternId,'Reflection Lists'))
            reflDict = reflSets[phaseName]
            if 'list' in str(type(reflDict)):       #patch for old reflection data
            #if isinstance(reflDict,list):       #patch for old reflection data
                RefData = {'RefList':[],'FF':[]}
                for ref in reflDict:
                    RefData['RefList'].append(ref[:11]+[ref[13],])
                    RefData['FF'].append(ref[14])
                RefData['RefList'] = np.array(RefData['RefList'])
                reflDict = RefData
        elif 'HKLF' in reflName:
            PatternId = G2gd.GetPatternTreeItemId(G2frame,G2frame.root, reflName)
            reflDict = G2frame.PatternTree.GetItemPyData(PatternId)[1]
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
            mapData.update(G2mth.ChargeFlip(data,reflDict,pgbar))
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

    def FillSelectPageMenu(TabSelectionIdDict, menuBar):
        '''Fill "Select tab" menu with menu items for each tab and assign
        bindings to the menu ietm to switch between phase tabs
        '''
        def OnSelectPage(event):
            'Called when an item is selected from the Select page menu'
            # lookup the menu item that called us and get its text
            tabname = TabSelectionIdDict.get(event.GetId())
            if not tabname:
                print 'Warning: menu item not in dict! id=',event.GetId()
                return                
            # find the matching tab
            for PageNum in range(G2frame.dataDisplay.GetPageCount()):
                if tabname == G2frame.dataDisplay.GetPageText(PageNum):
                    G2frame.dataDisplay.SetSelection(PageNum)
                    return
            else:
                print "Warning: tab "+tabname+" was not found"
        mid = menuBar.FindMenu('Select tab')
        menu = menuBar.GetMenu(mid)
        for ipage,page in enumerate(Pages):
            if menu.FindItem(page) < 0: # is tab already in menu?
                Id = wx.NewId()
                TabSelectionIdDict[Id] = page
                menu.Append(id=Id,kind=wx.ITEM_NORMAL,text=page)
                G2frame.Bind(wx.EVT_MENU, OnSelectPage, id=Id)
        
    def OnPageChanged(event):
        '''This is called every time that a Notebook tab button is pressed
        on a Phase data item window
        '''
        for page in G2frame.dataDisplay.gridList: # clear out all grids, forcing edits in progress to complete
            page.ClearGrid()
        wx.Frame.Unbind(G2frame.dataFrame,wx.EVT_SIZE) # ignore size events during this routine
        page = event.GetSelection()
        ChangePage(page)
        
    def ChangePage(page):
        text = G2frame.dataDisplay.GetPageText(page)
        if text == 'General':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.DataGeneral)
            UpdateGeneral()
        elif text == 'Data':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.DataMenu)
            G2ddG.UpdateDData(G2frame,DData,data)
            G2plt.PlotSizeStrainPO(G2frame,data,Start=True)            
        elif text == 'Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.AtomsMenu)
            FillAtomsGrid(Atoms)
        elif text == 'Wave Data' and data['General']['Type'] in ['modulated','magnetic']:
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.WavesData)
            UpdateWavesData()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Draw Options':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.DataDrawOptions)
            UpdateDrawOptions()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Draw Atoms':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.DrawAtomsMenu)
            UpdateDrawAtoms()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'RB Models':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RigidBodiesMenu)
            FillRigidBodyGrid()
        elif text == 'Map peaks':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.MapPeaksMenu)
            FillMapPeaksGrid()
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'MC/SA':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.MCSAMenu)
            UpdateMCSA()                        
            wx.CallAfter(G2plt.PlotStructure,G2frame,data,firstCall=True)
        elif text == 'Texture':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.TextureMenu)
            UpdateTexture()                        
            G2plt.PlotTexture(G2frame,data,Start=True)            
        elif text == 'Pawley reflections':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.PawleyMenu)
            FillPawleyReflectionsGrid()
        else:
            G2gd.SetDataMenuBar(G2frame)
    def FillMenus():
        '''Create the Select tab menus and bind to all menu items
        '''
        # General
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.DataGeneral)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnFourierMaps, id=G2gd.wxID_FOURCALC)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnSearchMaps, id=G2gd.wxID_FOURSEARCH)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnChargeFlip, id=G2gd.wxID_CHARGEFLIP)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnFourClear, id=G2gd.wxID_FOURCLEAR)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRunSingleMCSA, id=G2gd.wxID_SINGLEMCSA)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRunMultiMCSA, id=G2gd.wxID_MULTIMCSA)
        # Data
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.DataMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPwdrAdd, id=G2gd.wxID_PWDRADD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnHklfAdd, id=G2gd.wxID_HKLFADD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDataDelete, id=G2gd.wxID_DATADELETE)
        # Atoms
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.AtomsMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomAdd, id=G2gd.wxID_ATOMSEDITADD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomViewAdd, id=G2gd.wxID_ATOMSVIEWADD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomInsert, id=G2gd.wxID_ATOMSEDITINSERT)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomViewInsert, id=G2gd.wxID_ATOMVIEWINSERT)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAtomMove, id=G2gd.wxID_ATOMMOVE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, AtomDelete, id=G2gd.wxID_ATOMSEDITDELETE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, AtomRefine, id=G2gd.wxID_ATOMSREFINE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, AtomModify, id=G2gd.wxID_ATOMSMODIFY)
        G2frame.dataFrame.Bind(wx.EVT_MENU, AtomTransform, id=G2gd.wxID_ATOMSTRANSFORM)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnReloadDrawAtoms, id=G2gd.wxID_RELOADDRAWATOMS)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDistAngle, id=G2gd.wxID_ATOMSDISAGL)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDistAnglePrt, id=G2gd.wxID_ATOMSPDISAGL)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnIsoDistortCalc, id=G2gd.wxID_ISODISP)
        for id in G2frame.dataFrame.ReImportMenuId:     #loop over submenu items
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnReImport, id=id)                
        # Wave Data
        if data['General']['Type'] in ['modulated','magnetic']:
            FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.WavesData)
            G2frame.dataFrame.Bind(wx.EVT_MENU, On4DMapCompute, id=G2gd.wxID_4DMAPCOMPUTE)
        # Draw Options
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.DataDrawOptions)
        # Draw Atoms
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.DrawAtomsMenu)
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
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDrawDistVP, id=G2gd.wxID_DRAWDISTVP)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDrawDAT, id=G2gd.wxID_DRAWDISAGLTOR)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDrawPlane, id=G2gd.wxID_DRAWPLANE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRestraint, id=G2gd.wxID_DRAWRESTRBOND)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRestraint, id=G2gd.wxID_DRAWRESTRANGLE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRestraint, id=G2gd.wxID_DRAWRESTRPLANE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRestraint, id=G2gd.wxID_DRAWRESTRCHIRAL)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnDefineRB, id=G2gd.wxID_DRAWDEFINERB)
        # RB Models
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.RigidBodiesMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAutoFindResRB, id=G2gd.wxID_AUTOFINDRESRB)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRBAssign, id=G2gd.wxID_ASSIGNATMS2RB)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRBCopyParms, id=G2gd.wxID_COPYRBPARMS)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnGlobalResRBTherm, id=G2gd.wxID_GLOBALTHERM)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnGlobalResRBRef, id=G2gd.wxID_GLOBALRESREFINE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnRBRemoveAll, id=G2gd.wxID_RBREMOVEALL)
        # Map peaks
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.MapPeaksMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksMove, id=G2gd.wxID_PEAKSMOVE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksViewPoint, id=G2gd.wxID_PEAKSVIEWPT)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksDistVP, id=G2gd.wxID_PEAKSDISTVP)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksDA, id=G2gd.wxID_PEAKSDA)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnShowBonds, id=G2gd.wxID_SHOWBONDS)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksEquiv, id=G2gd.wxID_FINDEQVPEAKS)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksUnique, id=G2gd.wxID_PEAKSUNIQUE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksDelete, id=G2gd.wxID_PEAKSDELETE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPeaksClear, id=G2gd.wxID_PEAKSCLEAR)
        # MC/SA
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.MCSAMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnMCSAaddAtom, id=G2gd.wxID_ADDMCSAATOM)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnMCSAaddRB, id=G2gd.wxID_ADDMCSARB)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnMCSAclear, id=G2gd.wxID_CLEARMCSARB)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnMCSAmove, id=G2gd.wxID_MOVEMCSA)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnClearResults, id=G2gd.wxID_MCSACLEARRESULTS)
        # Texture
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.TextureMenu)
        #G2frame.dataFrame.Bind(wx.EVT_MENU, OnTextureRefine, id=G2gd.wxID_REFINETEXTURE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnTextureClear, id=G2gd.wxID_CLEARTEXTURE)
        # Pawley reflections
        FillSelectPageMenu(TabSelectionIdDict, G2frame.dataFrame.PawleyMenu)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyLoad, id=G2gd.wxID_PAWLEYLOAD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyEstimate, id=G2gd.wxID_PAWLEYESTIMATE)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPawleyUpdate, id=G2gd.wxID_PAWLEYUPDATE)
        
    # UpdatePhaseData execution starts here
#patch
    if 'RBModels' not in data:
        data['RBModels'] = {}
    if 'MCSA' not in data:
        data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
    #if isinstance(data['MCSA']['Results'],dict):
    if 'dict' in str(type(data['MCSA']['Results'])):
        data['MCSA']['Results'] = []
#end patch    

    global rbAtmDict   
    rbAtmDict = {}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    PhaseName = G2frame.PatternTree.GetItemText(Item)
    G2gd.SetDataMenuBar(G2frame)
    G2frame.dataFrame.SetLabel('Phase Data for '+PhaseName)
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataDisplay = G2gd.GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    G2frame.dataDisplay.gridList = [] # list of all grids in notebook
    Pages = []    
    wx.Frame.Unbind(G2frame.dataFrame,wx.EVT_SIZE) # ignore size events during this routine
    G2frame.dataDisplay.gridList = []
    General = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(General,'General')
    Pages.append('General')
    DData = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(DData,'Data')
    Pages.append('Data')
    Atoms = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.gridList.append(Atoms)
    G2frame.dataDisplay.AddPage(Atoms,'Atoms')
    Pages.append('Atoms')
    if data['General']['Type'] in ['modulated','magnetic']:
        waveData = wx.ScrolledWindow(G2frame.dataDisplay)
        G2frame.dataDisplay.AddPage(waveData,'Wave Data')
        Pages.append('Wave Data')        
    drawOptions = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(drawOptions,'Draw Options')
    Pages.append('Draw Options')
    drawAtoms = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.gridList.append(drawAtoms)
    G2frame.dataDisplay.AddPage(drawAtoms,'Draw Atoms')
    Pages.append('Draw Atoms')
    RigidBodies = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(RigidBodies,'RB Models')
    Pages.append('RB Models')
    MapPeaks = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.gridList.append(MapPeaks)    
    G2frame.dataDisplay.AddPage(MapPeaks,'Map peaks')
    Pages.append('Map peaks')
    MCSA = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(MCSA,'MC/SA')
    Pages.append('MC/SA')
    Texture = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(Texture,'Texture')
    Pages.append('Texture')
    G2frame.PawleyRefl = G2gd.GSGrid(G2frame.dataDisplay)
    G2frame.dataDisplay.gridList.append(G2frame.PawleyRefl)
    G2frame.dataDisplay.AddPage(G2frame.PawleyRefl,'Pawley reflections')
    Pages.append('Pawley reflections')
    G2frame.dataFrame.AtomCompute.ISOcalc.Enable('ISODISTORT' in data)
    G2frame.dataDisplay.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    FillMenus()
    if oldPage is None or oldPage == 0:
        ChangePage(0)
    elif oldPage:
        SetupGeneral()    # not sure why one might need this when moving from phase to phase; but does not hurt
        G2frame.dataDisplay.SetSelection(oldPage)
