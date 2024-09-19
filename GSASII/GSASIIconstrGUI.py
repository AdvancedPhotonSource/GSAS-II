# -*- coding: utf-8 -*-
#GSASIIconstrGUI - constraint GUI routines
########### SVN repository information ###################
# $Date: 2024-06-13 07:33:46 -0500 (Thu, 13 Jun 2024) $
# $Author: toby $
# $Revision: 5790 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIconstrGUI.py $
# $Id: GSASIIconstrGUI.py 5790 2024-06-13 12:33:46Z toby $
########### SVN repository information ###################
'''
Constraints and rigid bodies GUI routines follow.

'''
from __future__ import division, print_function
import platform
import sys
import copy
import os.path
import wx
import wx.grid as wg
import wx.lib.scrolledpanel as wxscroll
import wx.lib.gridmovers as gridmovers
import random as ran
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 5790 $")
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIstrIO as G2stIO
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIlattice as G2lat
import GSASIIdataGUI as G2gd
import GSASIIctrlGUI as G2G
import GSASIIfiles as G2fl
import GSASIIplot as G2plt
import GSASIIobj as G2obj
import GSASIIspc as G2spc
import GSASIIphsGUI as G2phG
import GSASIIIO as G2IO
import GSASIIscriptable as G2sc
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WACV = wx.ALIGN_CENTER_VERTICAL

class G2BoolEditor(wg.GridCellBoolEditor):
    '''Substitute for wx.grid.GridCellBoolEditor except toggles 
    grid items immediately when opened, updates grid & table contents after every
    item change
    '''
    def __init__(self):
        self.saveVals = None
        wx.grid.GridCellBoolEditor.__init__(self)


    def Create(self, parent, id, evtHandler):
        '''Create the editing control (wx.CheckBox) when cell is opened 
        for edit
        '''
        self._tc = wx.CheckBox(parent, -1, "")
        self._tc.Bind(wx.EVT_CHECKBOX, self.onCheckSet)
        self.SetControl(self._tc)
        if evtHandler:
            self._tc.PushEventHandler(evtHandler)

    def onCheckSet(self, event):
        '''Callback used when checkbox is toggled.
        Makes change to table immediately (creating event)
        '''
        if self.saveVals:
            self.ApplyEdit(*self.saveVals)

        
    def SetSize(self, rect):
        '''Set position/size the edit control within the cell's rectangle.
        '''
#        self._tc.SetDimensions(rect.x, rect.y, rect.width+2, rect.height+2, # older
        self._tc.SetSize(rect.x, rect.y, rect.width+2, rect.height+2,
                               wx.SIZE_ALLOW_MINUS_ONE)

    def BeginEdit(self, row, col, grid):
        '''Prepares the edit control by loading the initial 
        value from the table (toggles it since you would not 
        click on it if you were not planning to change it), 
        buts saves the original, pre-change value.
        Makes change to table immediately.
        Saves the info needed to make updates in self.saveVals.
        Sets the focus.
        '''
        if grid.GetTable().GetValue(row,col) not in [True,False]:
            return
        self.startValue = int(grid.GetTable().GetValue(row, col))
        self.saveVals = row, col, grid
        # invert state and set in editor
        if self.startValue:
            grid.GetTable().SetValue(row, col, 0)
            self._tc.SetValue(0)
        else:
            grid.GetTable().SetValue(row, col, 1)
            self._tc.SetValue(1)
        self._tc.SetFocus()
        self.ApplyEdit(*self.saveVals)

    def EndEdit(self, row, col, grid, oldVal=None):
        '''End editing the cell.  This is supposed to
        return None if the value has not changed, but I am not 
        sure that actually works. 
        '''
        val = int(self._tc.GetValue())
        if val != oldVal:   #self.startValue:?
            return val
        else:
            return None

    def ApplyEdit(self, row, col, grid):
        '''Save the value into the table, and create event. 
        Called after EndEdit(), BeginEdit and onCheckSet.
        '''
        val = int(self._tc.GetValue())
        grid.GetTable().SetValue(row, col, val) # update the table

    def Reset(self):
        '''Reset the value in the control back to its starting value.
        '''
        self._tc.SetValue(self.startValue)

    def StartingClick(self):
        '''This seems to be needed for BeginEdit to work properly'''
        pass

    def Destroy(self):
        "final cleanup"
        super(G2BoolEditor, self).Destroy()

    def Clone(self):
        'required'
        return G2BoolEditor()
    
class DragableRBGrid(wg.Grid):
    '''Simple grid implentation for display of rigid body positions.

    :param parent: frame or panel where grid will be placed
    :param dict rb: dict with atom labels, types and positions
    :param function onChange: a callback used every time a value in
      rb is changed. 
    '''
    def __init__(self, parent, rb, onChange=None):
        #wg.Grid.__init__(self, parent, wx.ID_ANY,size=(-1,200))
        wg.Grid.__init__(self, parent, wx.ID_ANY)
        self.SetTable(RBDataTable(rb,onChange), True)
        # Enable Row moving
        gridmovers.GridRowMover(self)
        self.Bind(gridmovers.EVT_GRID_ROW_MOVE, self.OnRowMove, self)
        self.SetColSize(0, 60)
        self.SetColSize(1, 40)
        self.SetColSize(2, 35)
        for r in range(len(rb['RBlbls'])):
            self.SetReadOnly(r,0,isReadOnly=True)
            self.SetCellEditor(r, 1, G2BoolEditor())            
            self.SetCellRenderer(r, 1, wg.GridCellBoolRenderer())
            self.SetReadOnly(r,2,isReadOnly=True)
            self.SetCellEditor(r,3, wg.GridCellFloatEditor())
            self.SetCellEditor(r,4, wg.GridCellFloatEditor())
            self.SetCellEditor(r,6, wg.GridCellFloatEditor())

    def OnRowMove(self,evt):
        'called when a row move needs to take place'
        frm = evt.GetMoveRow()          # Row being moved
        to = evt.GetBeforeRow()         # Before which row to insert
        self.GetTable().MoveRow(frm,to)
        
    def completeEdits(self):
        'complete any outstanding edits'
        if self.IsCellEditControlEnabled(): # complete any grid edits in progress
            #if GSASIIpath.GetConfigValue('debug'): print ('Completing grid edit')
            self.SaveEditControlValue()
            self.HideCellEditControl()
            self.DisableCellEditControl()
            
class RBDataTable(wg.GridTableBase):
    '''A Table to support :class:`DragableRBGrid`
    '''
    def __init__(self,rb,onChange):
        wg.GridTableBase.__init__(self)
        self.colLabels = ['Label','Select','Type','x','y','z']
        self.coords = rb['RBcoords']
        self.labels = rb['RBlbls']
        self.types = rb['RBtypes']
        self.index = rb['RBindex']
        self.select = rb['RBselection']
        self.onChange = onChange

    # required methods
    def GetNumberRows(self):
        return len(self.labels)
    def GetNumberCols(self):
        return len(self.colLabels)
    def IsEmptyCell(self, row, col):
        return False
    def GetValue(self, row, col):
        row = self.index[row]
        if col == 0:
            return self.labels[row]
        elif col == 1:
            if self.select[row]:
                return '1'
            else:
                return ''
        elif col == 2:
            return self.types[row]
        else:
            return '{:.5f}'.format(self.coords[row][col-3])
    def SetValue(self, row, col, value):
        row = self.index[row]
        try:
            if col == 0:
                self.labels[row] = value
            elif col == 1:
                self.select[row] = bool(value)
            elif col == 2:
                self.types[row] = value
            else:
                self.coords[row][col-3] = float(value)
        except:
            pass
        if self.onChange:
            self.onChange()
    # Display column & row labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self,row):
        return str(row)

    # Implement "row movement" by updating the pointer array
    def MoveRow(self,frm,to):
        grid = self.GetView()
        if grid:
            move = self.index[frm]
            del self.index[frm]
            if frm > to:
                self.index.insert(to,move)
            else:
                self.index.insert(to-1,move)
            
            # Notify the grid
            grid.BeginBatch()
            msg = wg.GridTableMessage(
                    self, wg.GRIDTABLE_NOTIFY_ROWS_DELETED, frm, 1
                    )
            grid.ProcessTableMessage(msg)
            msg = wg.GridTableMessage(
                    self, wg.GRIDTABLE_NOTIFY_ROWS_INSERTED, to, 1
                    )
            grid.ProcessTableMessage(msg)
            grid.EndBatch()
        if self.onChange:
            self.onChange()

# def MakeDrawAtom(data,atom):
#     'Convert atom to format needed to draw it'
#     generalData = data['General']
#     deftype = G2obj.validateAtomDrawType(
#         GSASIIpath.GetConfigValue('DrawAtoms_default'),generalData)
#     if generalData['Type'] in ['nuclear','faulted',]:
#         atomInfo = [atom[:2]+atom[3:6]+['1']+[deftype]+
#                     ['']+[[255,255,255]]+atom[9:]+[[],[]]][0]
#     ct,cs = [1,8]         #type & color
#     atNum = generalData['AtomTypes'].index(atom[ct])
#     atomInfo[cs] = list(generalData['Color'][atNum])
#     return atomInfo

class ConstraintDialog(wx.Dialog):
    '''Window to edit Constraint values
    '''
    def __init__(self,parent,title,text,data,separator='*',varname="",varyflag=False):
        wx.Dialog.__init__(self,parent,-1,'Edit '+title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.data = data[:]
        self.newvar = [varname,varyflag]
        panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        topLabl = wx.StaticText(panel,-1,text)
        mainSizer.Add((10,10),1)
        mainSizer.Add(topLabl,0,wx.LEFT,10)
        mainSizer.Add((10,10),1)
        dataGridSizer = wx.FlexGridSizer(cols=3,hgap=2,vgap=2)
        self.OkBtn = wx.Button(panel,wx.ID_OK)
        for id in range(len(self.data)):
            lbl1 = lbl = str(self.data[id][1])
            if lbl[-1] != '=': lbl1 = lbl + ' ' + separator + ' '
            name = wx.StaticText(panel,wx.ID_ANY,lbl1,style=wx.ALIGN_RIGHT)
            scale = G2G.ValidatedTxtCtrl(panel,self.data[id],0,OKcontrol=self.DisableOK)
            dataGridSizer.Add(name,0,wx.LEFT|wx.RIGHT|WACV,5)
            dataGridSizer.Add(scale,0,wx.RIGHT,3)
            if ':' in lbl:
                dataGridSizer.Add(
                    wx.StaticText(panel,-1,G2obj.fmtVarDescr(lbl)),
                    0,wx.RIGHT|WACV,3)
            else:
                dataGridSizer.Add((-1,-1))
        if title == 'New Variable':
            name = wx.StaticText(panel,wx.ID_ANY,"New variable's\nname (optional)",
                style=wx.ALIGN_CENTER)
            scale = G2G.ValidatedTxtCtrl(panel,self.newvar,0,notBlank=False)
            dataGridSizer.Add(name,0,wx.LEFT|wx.RIGHT|WACV,5)
            dataGridSizer.Add(scale,0,wx.RIGHT|WACV,3)
            self.refine = wx.CheckBox(panel,label='Refine?')
            self.refine.SetValue(self.newvar[1]==True)
            self.refine.Bind(wx.EVT_CHECKBOX, self.OnCheckBox)
            dataGridSizer.Add(self.refine,0,wx.RIGHT|WACV,3)
            
        mainSizer.Add(dataGridSizer,0,wx.EXPAND)
        self.OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        self.OkBtn.SetDefault()
        cancelBtn = wx.Button(panel,wx.ID_CANCEL)
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(self.OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)

        mainSizer.Add(btnSizer,0,wx.EXPAND, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()
        self.Fit()
        self.CenterOnParent()
        
    def DisableOK(self,setting):
        for id in range(len(self.data)):  # coefficient cannot be zero
            try:
                if abs(self.data[id][0]) < 1.e-20:
                    setting  = False
                    break
            except:
                pass
        if setting:
            self.OkBtn.Enable()
        else:
            self.OkBtn.Disable()

    def OnCheckBox(self,event):
        self.newvar[1] = self.refine.GetValue()

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)              

    def GetData(self):
        return self.data
        
#####  Constraints ################################################################################           
def CheckConstraints(G2frame,Phases,Histograms,data,newcons=[],reqVaryList=None,seqhst=None,seqmode='use-all'):
    '''Load constraints & check them for errors. 

    N.B. Equivalences based on symmetry (etc.)
    are generated by running :func:`GSASIIstrIO.GetPhaseData`.

    When reqVaryList is included (see WarnConstraintLimit) then 
    parameters with limits are checked against constraints and a 
    warning is shown.
    '''
    G2mv.InitVars()
    #Find all constraints
    constrDict = []
    for key in data:
        if key.startswith('_'): continue
        constrDict += data[key]
    if newcons:
        constrDict = constrDict + newcons
    constrDict, fixedList, ignored = G2mv.ProcessConstraints(constrDict, seqhst=seqhst, seqmode=seqmode)
    parmDict = {}
    # generate symmetry constraints to check for conflicts
    rigidbodyDict = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame, G2frame.root, 'Rigid bodies'))
    rbIds = rigidbodyDict.get('RBIds', {'Vector': [], 'Residue': [],'Spin':[]})
    rbVary, rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict, Print=False)
    parmDict.update(rbDict)
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave) = \
        G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False) # generates atom symmetry constraints
    parmDict.update(phaseDict)
    # get Hist and HAP info
    hapVary, hapDict, controlDict = G2stIO.GetHistogramPhaseData(Phases, Histograms, Print=False, resetRefList=False)
    parmDict.update(hapDict)
    histVary, histDict, controlDict = G2stIO.GetHistogramData(Histograms, Print=False)
    parmDict.update(histDict)
    
    # TODO: twining info needed?
    #TwConstr,TwFixed = G2stIO.makeTwinFrConstr(Phases,Histograms,hapVary)
    #constrDict += TwConstr
    #fixedList += TwFixed
    varyList = rbVary+phaseVary+hapVary+histVary

    msg = G2mv.EvaluateMultipliers(constrDict,parmDict)
    if msg:
        return 'Unable to interpret multiplier(s): '+msg,''
    if reqVaryList:
        varyList = reqVaryList[:]
    errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict) # changes varyList
    
    impossible = []
    if reqVaryList:
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        G2mv.Map2Dict(parmDict,varyList)   # changes parmDict & varyList
        # check for limits on dependent vars
        consVars = [i for i in reqVaryList if i not in varyList]
        impossible = set(
                [str(i) for i in Controls['parmMinDict'] if i in consVars] + 
                [str(i) for i in Controls['parmMaxDict'] if i in consVars])
        if impossible:
            msg = ''
            for i in sorted(impossible):
                if msg: msg += ', '
                msg += i
            msg =  ' &'.join(msg.rsplit(',',1))
            msg = ('Note: limits on variable(s) '+msg+
            ' will be ignored because they are constrained.')
            G2G.G2MessageBox(G2frame,msg,'Limits ignored for constrained vars')
    else:
        G2mv.Map2Dict(parmDict,varyList)   # changes varyList
    return errmsg,warnmsg

def UpdateConstraints(G2frame, data, selectTab=None, Clear=False):
    '''Called when Constraints tree item is selected.
    Displays the constraints in the data window
    '''
        
    def FindEquivVarb(name,nameList):
        'Creates a list of variables appropriate to constrain with name'
        outList = []
        #phlist = []
        items = name.split(':')
        namelist = [items[2],]
        if 'dA' in name:
            namelist = ['dAx','dAy','dAz']
        elif 'AU' in name:
            namelist = ['AUiso','AU11','AU22','AU33','AU12','AU13','AU23']
        elif 'Akappa' in name:
            namelist = ['Akappa%d'%i for i in range(6)]
        elif 'ANe' in name:
            namelist = ['ANe0','ANe1']
        elif 'AD(1' in name:
            orb = name.split(':')[2][-1]
            namelist = ['AD(1,0)'+orb,'AD(1,1)'+orb,'AD(1,-1)'+orb]
        elif 'AD(2' in name:
            orb = name.split(':')[2][-1]
            namelist = ['AD(2,0)'+orb,'AD(2,1)'+orb,'AD(2,-1)'+orb,'AD(2,2)'+orb,'AD(2,-2)'+orb]
        elif 'AD(4' in name:
            orb = name.split(':')[2][-1]
            namelist = ['AD(4,0)'+orb,'AD(4,1)'+orb,'AD(4,-1)'+orb,'AD(4,2)'+orb,'AD(4,-2)'+orb,
                'AD(4,3)'+orb,'AD(4,-3)'+orb,'AD(4,4)'+orb,'AD(4,-4)'+orb]
        elif 'AM' in name:
            namelist = ['AMx','AMy','AMz']
        elif items[-1] in ['A0','A1','A2','A3','A4','A5']:
            namelist = ['A0','A1','A2','A3','A4','A5']
        elif items[-1] in ['D11','D22','D33','D12','D13','D23']:
            namelist = ['D11','D22','D33','D12','D13','D23']
        elif 'Tm' in name:
            namelist = ['Tmin','Tmax']
        elif 'MX' in name or 'MY' in name or 'MZ' in name:
            namelist = ['MXcos','MYcos','MZcos','MXsin','MYsin','MZsin']
        elif 'mV' in name:
            namelist = ['mV0','mV1','mV2']
        elif 'Debye' in name or 'BkPk' in name:   #special cases for Background fxns
            dbname = name.split(';')[0].split(':')[2]
            return [item for item in nameList if dbname in item]
        elif 'RB' in name:
            rbfx = 'RB'+items[2][2]
            if 'T' in name and 'Tr' not in name:
                namelist = [rbfx+'T11',rbfx+'T22',rbfx+'T33',rbfx+'T12',rbfx+'T13',rbfx+'T23']
            if 'L' in name:
                namelist = [rbfx+'L11',rbfx+'L22',rbfx+'L33',rbfx+'L12',rbfx+'L13',rbfx+'L23']
            if 'S' in name:
                namelist = [rbfx+'S12',rbfx+'S13',rbfx+'S21',rbfx+'S23',rbfx+'S31',rbfx+'S32',rbfx+'SAA',rbfx+'SBB']
            if 'U' in name:
                namelist = [rbfx+'U',]
            if 'Tr' in name:
                namelist = [rbfx+'Tr',]

        for item in nameList:
            keys = item.split(':')
            #if keys[0] not in phlist:
            #    phlist.append(keys[0])
            if items[1] == '*' and keys[2] in namelist: # wildcard -- select only sequential options
                keys[1] = '*'
                mitem = ':'.join(keys)
                if mitem == name: continue
                if mitem not in outList: outList.append(mitem)
            elif (keys[2] in namelist or keys[2].split(';')[0] in namelist) and item != name:
                outList.append(item)
        return outList
        
    def SelectVarbs(page,FrstVarb,varList,legend,constType):
        '''Select variables used in constraints after one variable has
        been selected. This routine determines the appropriate variables to be
        used based on the one that has been selected and asks for more to be added.

        It then creates the constraint and adds it to the constraints list.
        
        Called from OnAddEquivalence, OnAddFunction & OnAddConstraint (all but
        OnAddHold)

        :param list page: defines calling page -- type of variables to be used
        :parm GSASIIobj.G2VarObj FrstVarb: reference to first selected variable
        :param list varList: list of other appropriate variables to select from
        :param str legend: header for selection dialog
        :param str constType: type of constraint to be generated
        :returns: a constraint, as defined in
          :ref:`GSASIIobj <Constraint_definitions_table>`
        '''
        choices = [[i]+list(G2obj.VarDescr(i)) for i in varList]
        meaning = G2obj.getDescr(FrstVarb.name)
        if not meaning:
            meaning = "(no definition found!)"
        l = str(FrstVarb).split(':')
        # make lists of phases & histograms to iterate over
        phaselist = [l[0]]
        if l[0]:
            phaselbl = ['phase #'+l[0]]
            if len(Phases) > 1:
                phaselist += ['all'] 
                phaselbl += ['all phases']
        else:
            phaselbl = ['']
        histlist = [l[1]]
        if l[1] == '*':
            pass
        elif l[1]:
            histlbl = ['histogram #'+l[1]]
            if len(Histograms) > 1:
                histlist += ['all']
                histlbl += ['all histograms']
                typ = Histograms[G2obj.LookupHistName(l[1])[0]]['Instrument Parameters'][0]['Type'][1]
                i = 0
                for hist in Histograms:
                    if Histograms[hist]['Instrument Parameters'][0]['Type'][1] == typ: i += 1
                if i > 1:
                    histlist += ['all='+typ]
                    histlbl += ['all '+typ+' histograms']
        else:
            histlbl = ['']
        # make a list of equivalent parameter names
        nameList = [FrstVarb.name]
        for var in varList:
            nam = var.split(":")[2]
            if nam not in nameList: nameList += [nam]
        # add "wild-card" names to the list of variables
        if l[1] == '*':
            pass
        elif page[1] == 'phs':
            if 'RB' in FrstVarb.name:
                pass
            elif FrstVarb.atom is None:
                for nam in nameList:
                    for ph,plbl in zip(phaselist,phaselbl):
                        if plbl: plbl = 'For ' + plbl
                        var = ph+"::"+nam
                        if var == str(FrstVarb) or var in varList: continue
                        varList += [var]
                        choices.append([var,plbl,meaning])
            else:
                for nam in nameList:
                    for ph,plbl in zip(phaselist,phaselbl):
                        if plbl: plbl = ' in ' + plbl
                        for atype in ['']+TypeList:
                            if atype:
                                albl = "For "+atype+" atoms"
                                akey = "all="+atype                        
                            else:
                                albl = "For all atoms"
                                akey = "all"
                            var = ph+"::"+nam+":"+akey
                            if var == str(FrstVarb) or var in varList: continue
                            varList += [var]
                            choices.append([var,albl+plbl,meaning])
        elif page[1] == 'hap':
            if FrstVarb.name == "Scale":
                meaning = "Phase fraction"
            for nam in nameList:
                for ph,plbl in zip(phaselist,phaselbl):
                    if plbl: plbl = 'For ' + plbl
                    for hst,hlbl in zip(histlist,histlbl):
                        if hlbl:
                            if plbl:
                                hlbl = ' in ' + hlbl
                            else:
                                hlbl = 'For ' + hlbl                                
                            var = ph+":"+hst+":"+nam
                            if var == str(FrstVarb) or var in varList: continue
                            varList += [var]
                            choices.append([var,plbl+hlbl,meaning])
        elif page[1] == 'hst':
            if FrstVarb.name == "Scale":
                meaning = "Scale factor"
            for nam in nameList:
                for hst,hlbl in zip(histlist,histlbl):
                    if hlbl:
                        hlbl = 'For ' + hlbl                                
                        var = ":"+hst+":"+nam
                        if var == str(FrstVarb) or var in varList: continue
                        varList += [var]
                        choices.append([var,hlbl,meaning])
        elif page[1] == 'glb' or page[1] == 'sym':
            pass
        else:
            raise Exception('Unknown constraint page '+ page[1])                    
        if len(choices):
            l1 = l2 = 1
            for i1,i2,i3 in choices:
                l1 = max(l1,len(i1))
                l2 = max(l2,len(i2))
            fmt = "{:"+str(l1)+"s} {:"+str(l2)+"s} {:s}"
            atchoices = [fmt.format(*i1) for i1 in choices] # reformat list as str with columns
            dlg = G2G.G2MultiChoiceDialog(
                G2frame,legend,
                'Constrain '+str(FrstVarb)+' with...',atchoices,
                toggle=False,size=(625,400),monoFont=True)
            dlg.CenterOnParent()
            res = dlg.ShowModal()
            Selections = dlg.GetSelections()[:]
            dlg.Destroy()
            if res != wx.ID_OK: return []
            if len(Selections) == 0:
                dlg = wx.MessageDialog(
                    G2frame,
                    'No variables were selected to include with '+str(FrstVarb),
                    'No variables')
                dlg.CenterOnParent()
                dlg.ShowModal()
                dlg.Destroy()
                return []
        else:
            dlg = wx.MessageDialog(
                G2frame,
                'There are no appropriate variables to include with '+str(FrstVarb),
                'No variables')
            dlg.CenterOnParent()
            dlg.ShowModal()
            dlg.Destroy()
            return []
        # now process the variables provided by the user
        varbs = [str(FrstVarb),] # list of selected variables
        for sel in Selections:
            var = varList[sel]
            # phase(s) included
            l = var.split(':')
            if l[0] == "all":
                phlist = [str(Phases[phase]['pId']) for phase in Phases]
            else:
                phlist = [l[0]]
            # histogram(s) included
            if l[1] == "all":
                hstlist = [str(Histograms[hist]['hId']) for hist in Histograms]
            elif '=' in l[1]:
                htyp = l[1].split('=')[1]
                hstlist = [str(Histograms[hist]['hId']) for hist in Histograms if 
                           Histograms[hist]['Instrument Parameters'][0]['Type'][1] == htyp]
            else:
                hstlist = [l[1]]
            if len(l) == 3:
                for ph in phlist:
                    for hst in hstlist:
                        var = ph + ":" + hst + ":" + l[2]
                        if var in varbs: continue
                        varbs.append(var)
            else: # constraints with atoms or rigid bodies
                if len(l) == 5: # rigid body parameter
                    var = ':'.join(l)
                    if var in varbs: continue
                    varbs.append(var)
                elif l[3] == "all":
                    for ph in phlist:
                        key = G2obj.LookupPhaseName(ph)[0]
                        for hst in hstlist: # should be blank
                            for iatm,at in enumerate(Phases[key]['Atoms']):
                                var = ph + ":" + hst + ":" + l[2] + ":" + str(iatm)
                                if var in varbs: continue
                                varbs.append(var)
                elif '=' in l[3]:
                    for ph in phlist:
                        key = G2obj.LookupPhaseName(ph)[0]
                        cx,ct,cs,cia = Phases[key]['General']['AtomPtrs']
                        for hst in hstlist: # should be blank
                            atyp = l[3].split('=')[1]
                            for iatm,at in enumerate(Phases[key]['Atoms']):
                                if at[ct] != atyp: continue
                                var = ph + ":" + hst + ":" + l[2] + ":" + str(iatm)
                                if var in varbs: continue
                                varbs.append(var)
                else:
                    for ph in phlist:
                        key = G2obj.LookupPhaseName(ph)[0]
                        for hst in hstlist: # should be blank
                            var = ph + ":" + hst + ":" + l[2] + ":" + l[3]
                            if var in varbs: continue
                            varbs.append(var)
        if len(varbs) >= 1 or 'constraint' in constType:
            constr = [[1.0,FrstVarb]]
            for item in varbs[1:]:
                constr += [[1.0,G2obj.G2VarObj(item)]]
            if 'equivalence' in constType:
                return [constr+[None,None,'e']]
            elif 'function' in constType:
                return [constr+[None,False,'f']]
            elif 'constraint' in constType:
                return [constr+[1.0,None,'c']]
            else:
                raise Exception('Unknown constraint type: '+str(constType))
        else:
            dlg = wx.MessageDialog(
                G2frame,
                'There are no selected variables to include with '+str(FrstVarb),
                'No variables')
            dlg.CenterOnParent()
            dlg.ShowModal()
            dlg.Destroy()
        return []
                
    def CheckAddedConstraint(newcons):
        '''Check a new constraint that has just been input.
        If there is an error display a message and discard the last entry

        Since the varylist is not available, no warning messages
        should be generated here

        :returns: True if constraint should be added
        '''
        
        errmsg,warnmsg = CheckConstraints(G2frame,Phases,Histograms,data,newcons,seqhst=seqhistnum,seqmode=seqmode)
        if errmsg:
            G2frame.ErrorDialog('Constraint Error',
                'Error with newly added constraint:\n'+errmsg+
                '\nIgnoring newly added constraint',parent=G2frame)
            # reset error status
            errmsg,warnmsg = CheckConstraints(G2frame,Phases,Histograms,data,seqhst=seqhistnum,seqmode=seqmode)
            if errmsg:
                print (errmsg)
                print (G2mv.VarRemapShow([],True))
            return False
        elif warnmsg:
            print ('Warning after constraint addition:\n'+warnmsg)
            txt = 'Warning noted after adding constraint (this may be OK):\n'
            txt += warnmsg + '\n\nKeep this addition?'
            ans = G2G.ShowScrolledInfo(header='Constraint Warning',txt=txt,
                buttonlist=[wx.ID_YES,wx.ID_NO],parent=G2frame,height=250)
            if ans == wx.ID_NO: return False
        return True

    def WarnConstraintLimit():
        '''Check if constraints reference variables with limits.
        Displays a warning message, but does nothing
        '''
        parmDict,reqVaryList = G2frame.MakeLSParmDict()
        try:
            errmsg,warnmsg = CheckConstraints(G2frame,Phases,Histograms,data,[],reqVaryList,seqhst=seqhistnum,seqmode=seqmode)
        except Exception as msg:
            if GSASIIpath.GetConfigValue('debug'): 
                import traceback
                print (traceback.format_exc())
            return 'CheckConstraints error retrieving parameter\nError='+str(msg),''
        return errmsg,warnmsg

    def CheckChangedConstraint():
        '''Check all constraints after an edit has been made.
        If there is an error display a message and reject the change.

        Since the varylist is not available, no warning messages
        should be generated.
        
        :returns: True if the edit should be retained
        '''
        errmsg,warnmsg = CheckConstraints(G2frame,Phases,Histograms,data,seqhst=seqhistnum,seqmode=seqmode)
        if errmsg:
            G2frame.ErrorDialog('Constraint Error',
                'Error after editing constraint:\n'+errmsg+
                '\nDiscarding last constraint edit',parent=G2frame)
            # reset error status
            errmsg,warnmsg = CheckConstraints(G2frame,Phases,Histograms,data,seqhst=seqhistnum,seqmode=seqmode)
            if errmsg:
                print (errmsg)
                print (G2mv.VarRemapShow([],True))
            return False
        elif warnmsg:
            print ('Warning after constraint edit:\n'+warnmsg)
            ans = G2G.ShowScrolledInfo(header='Constraint Warning',
                    txt='Warning noted after last constraint edit:\n'+warnmsg+
                    '\n\nKeep this change?',
                    buttonlist=[wx.ID_YES,wx.ID_NO],parent=G2frame,height=250)
            if ans == wx.ID_NO: return False
        return True
             
    def PageSelection(page):
        'Decode page reference'
        if page[1] == "phs":
            vartype = "phase"
            varList = G2obj.removeNonRefined(phaseList)  # remove any non-refinable prms from list
            constrDictEnt = 'Phase'
        elif page[1] == "hap":
            vartype = "Histogram*Phase"
            varList = G2obj.removeNonRefined(hapList)  # remove any non-refinable prms from list
            constrDictEnt = 'HAP'
        elif page[1] == "hst":
            vartype = "Histogram"
            varList = G2obj.removeNonRefined(histList)  # remove any non-refinable prms from list
            constrDictEnt = 'Hist'
        elif page[1] == "glb":
            vartype = "Global"
            varList = G2obj.removeNonRefined(globalList)   # remove any non-refinable prms from list

            constrDictEnt = 'Global'
        elif page[1] == "sym":
            return None,None,None
        else:
            raise Exception('Should not happen!')
        return vartype,varList,constrDictEnt

    def OnAddHold(event):
        '''Create a new Hold constraint

        Hold constraints allows the user to select one variable (the list of available
        variables depends on which tab is currently active). 
        '''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        varList = G2obj.SortVariables(varList)
        title1 = "Hold "+vartype+" variable"
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
        l2 = l1 = 1
        for i in varList:
            l1 = max(l1,len(i))
            loc,desc = G2obj.VarDescr(i)
            l2 = max(l2,len(loc))
        fmt = "{:"+str(l1)+"s} {:"+str(l2)+"s} {:s}"
        varListlbl = [fmt.format(i,*G2obj.VarDescr(i)) for i in varList]
        #varListlbl = ["("+i+") "+G2obj.fmtVarDescr(i) for i in varList]
        legend = "Select variables to hold (Will not be varied, even if vary flag is set)"
        dlg = G2G.G2MultiChoiceDialog(G2frame,
            legend,title1,varListlbl,toggle=False,size=(625,400),monoFont=True)
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            for sel in dlg.GetSelections():
                Varb = varList[sel]
                VarObj = G2obj.G2VarObj(Varb)
                newcons = [[[0.0,VarObj],None,None,'h']]
                if CheckAddedConstraint(newcons):
                    data[constrDictEnt] += newcons
        dlg.Destroy()
        #wx.CallAfter(OnPageChanged,None)
        wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True)
        
    def OnAddEquivalence(event):
        '''add an Equivalence constraint'''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        title1 = "Create equivalence constraint between "+vartype+" variables"
        title2 = "Select additional "+vartype+" variable(s) to be equivalent with "
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
#        legend = "Select variables to make equivalent (only one of the variables will be varied when all are set to be varied)"
        GetAddVars(page,title1,title2,varList,constrDictEnt,'equivalence')
        
    def OnAddAtomEquiv(event):
        ''' Add equivalences between all parameters on atoms '''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        title1 = "Setup equivalent atom variables"
        title2 = "Select additional atoms(s) to be equivalent with "
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
#        legend = "Select atoms to make equivalent (only one of the atom variables will be varied when all are set to be varied)"
        GetAddAtomVars(page,title1,title2,varList,constrDictEnt,'equivalence')
        
    def OnAddRiding(event):
        ''' Add riding equivalences between all parameters on atoms  - not currently used'''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        title1 = "Setup riding atoms "
        title2 = "Select additional atoms(s) to ride on "
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
#        legend = "Select atoms to ride (only one of the atom variables will be varied when all are set to be varied)"
        GetAddAtomVars(page,title1,title2,varList,constrDictEnt,'riding')
   
    def OnAddFunction(event):
        '''add a Function (new variable) constraint'''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        title1 = "Setup new variable based on "+vartype+" variables"
        title2 = "Include additional "+vartype+" variable(s) to be included with "
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
#        legend = "Select variables to include in a new variable (the new variable will be varied when all included variables are varied)"
        GetAddVars(page,title1,title2,varList,constrDictEnt,'function')
                        
    def OnAddConstraint(event):
        '''add a constraint equation to the constraints list'''
        page = G2frame.Page
        vartype,varList,constrDictEnt = PageSelection(page)
        if vartype is None: return
        title1 = "Creating constraint on "+vartype+" variables"
        title2 = "Select additional "+vartype+" variable(s) to include in constraint with "
        if not varList:
            G2frame.ErrorDialog('No variables','There are no variables of type '+vartype,
                parent=G2frame)
            return
#        legend = "Select variables to include in a constraint equation (the values will be constrainted to equal a specified constant)"
        GetAddVars(page,title1,title2,varList,constrDictEnt,'constraint')

    def GetAddVars(page,title1,title2,varList,constrDictEnt,constType):
        '''Get the variables to be added for OnAddEquivalence, OnAddFunction,
        and OnAddConstraint. Then create and check the constraint.
        '''
        #varListlbl = ["("+i+") "+G2obj.fmtVarDescr(i) for i in varList]
        if constType == 'equivalence':
            omitVars = G2mv.GetDependentVars()
        else:
            omitVars = []
        varList = G2obj.SortVariables([i for i in varList if i not in omitVars])
        l2 = l1 = 1
        for i in varList:
            l1 = max(l1,len(i))
            loc,desc = G2obj.VarDescr(i)
            l2 = max(l2,len(loc))
        fmt = "{:"+str(l1)+"s} {:"+str(l2)+"s} {:s}"
        varListlbl = [fmt.format(i,*G2obj.VarDescr(i)) for i in varList]
        dlg = G2G.G2SingleChoiceDialog(G2frame,'Select 1st variable:',
            title1,varListlbl,monoFont=True,size=(625,400))
        dlg.CenterOnParent()
        if dlg.ShowModal() == wx.ID_OK:
            if constType == 'equivalence':
                omitVars = G2mv.GetDependentVars() + G2mv.GetIndependentVars()
            sel = dlg.GetSelection()
            FrstVarb = varList[sel]
            VarObj = G2obj.G2VarObj(FrstVarb)
            moreVarb = G2obj.SortVariables(FindEquivVarb(FrstVarb,[i for i in varList if i not in omitVars]))
            newcons = SelectVarbs(page,VarObj,moreVarb,title2+FrstVarb,constType)
            if len(newcons) > 0:
                if CheckAddedConstraint(newcons):
                    data[constrDictEnt] += newcons
        dlg.Destroy()
        WarnConstraintLimit()
#        wx.CallAfter(OnPageChanged,None)
        wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True)
                        
    def FindNeighbors(phase,FrstName,AtNames):
        General = phase['General']
        cx,ct,cs,cia = General['AtomPtrs']
        Atoms = phase['Atoms']
        atNames = [atom[ct-1] for atom in Atoms]
        Cell = General['Cell'][1:7]
        Amat,Bmat = G2lat.cell2AB(Cell)
        atTypes = General['AtomTypes']
        Radii = np.array(General['BondRadii'])
        AtInfo = dict(zip(atTypes,Radii)) #or General['BondRadii']
        Orig = atNames.index(FrstName.split()[1])
        OType = Atoms[Orig][ct]
        XYZ = G2mth.getAtomXYZ(Atoms,cx)        
        Neigh = []
        Dx = np.inner(Amat,XYZ-XYZ[Orig]).T
        dist = np.sqrt(np.sum(Dx**2,axis=1))
        sumR = AtInfo[OType]+0.5    #H-atoms only!
        IndB = ma.nonzero(ma.masked_greater(dist-0.85*sumR,0.))
        for j in IndB[0]:
            if j != Orig:
                Neigh.append(AtNames[j])
        return Neigh
        
    def GetAddAtomVars(page,title1,title2,varList,constrDictEnt,constType):
        '''Get the atom variables to be added for OnAddAtomEquiv. Then create and 
        check the constraints. Riding for H atoms only.
        '''
        Atoms = {G2obj.VarDescr(i)[0]:[] for i in varList if 'Atom' in G2obj.VarDescr(i)[0]}
        for item in varList:
            atName = G2obj.VarDescr(item)[0]
            if atName in Atoms:
                Atoms[atName].append(item)
        AtNames = list(Atoms.keys())
        AtNames.sort()
        dlg = G2G.G2SingleChoiceDialog(G2frame,'Select 1st atom:',
            title1,AtNames,monoFont=True,size=(625,400))
        dlg.CenterOnParent()
        FrstAtom = ''
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            FrstAtom = AtNames[sel]
            if 'riding' in constType:
                phaseName = (FrstAtom.split(' in ')[1]).strip()
                phase = Phases[phaseName]
                AtNames = FindNeighbors(phase,FrstAtom,AtNames)
            else:
                AtNames.remove(FrstAtom)
        dlg.Destroy()
        if FrstAtom == '':
            print ('no atom selected')
            return
        dlg = G2G.G2MultiChoiceDialog(
            G2frame,title2+FrstAtom,
            'Constrain '+str(FrstAtom)+' with...',AtNames,
            toggle=False,size=(625,400),monoFont=True)
        if dlg.ShowModal() == wx.ID_OK:
            Selections = dlg.GetSelections()[:]
        else:
            print ('no target atom selected')
            dlg.Destroy()
            return
        dlg.Destroy()
        for name in Atoms[FrstAtom]:
            newcons = []
            constr = []
            if 'riding' in constType:
                if 'AUiso' in name:
                    constr = [[1.0,G2obj.G2VarObj(name)]]
                elif 'AU11' in name:
                    pass
                elif 'AU' not in name:
                    constr = [[1.0,G2obj.G2VarObj(name)]]
            else:
                constr = [[1.0,G2obj.G2VarObj(name)]]
            pref = ':'+name.rsplit(':',1)[0].split(':',1)[1]    #get stuff between phase id & atom id
            for sel in Selections:
                name2 = Atoms[AtNames[sel]][0]
                pid = name2.split(':',1)[0]                     #get phase id for 2nd atom
                id = name2.rsplit(':',1)[-1]                    #get atom no. for 2nd atom
                if 'riding' in constType:
                    pref = pid+pref
                    if 'AUiso' in pref:
                        parts = pref.split('AUiso')
                        constr += [[1.2,G2obj.G2VarObj('%s:%s'%(parts[0]+'AUiso',id))]]
                    elif 'AU' not in pref:
                        constr += [[1.0,G2obj.G2VarObj('%s:%s'%(pref,id))]]
                else:
                    constr += [[1.0,G2obj.G2VarObj('%s:%s'%(pid+pref,id))]]
            if not constr:
                continue
            if 'frac' in pref and 'riding' not in constType:
                newcons = [constr+[1.0,None,'c']]
            else:
                newcons = [constr+[None,None,'e']]
            if len(newcons) > 0:
                if CheckAddedConstraint(newcons):
                    data[constrDictEnt] += newcons
        WarnConstraintLimit()
#        wx.CallAfter(OnPageChanged,None)
        wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True)
                        
    def MakeConstraintsSizer(name,panel):
        '''Creates a sizer displaying all of the constraints entered of
        the specified type.

        :param str name: the type of constraints to be displayed ('HAP',
          'Hist', 'Phase', 'Global', 'Sym-Generated')
        :param wx.Panel panel: parent panel for sizer
        :returns: wx.Sizer created by method
        '''
        if name == 'Sym-Generated':         #show symmetry generated constraints
            Sizer1 =  wx.BoxSizer(wx.VERTICAL)
            if symHolds:
                Sizer1.Add(wx.StaticText(panel,wx.ID_ANY,
                    'Position variables fixed by space group symmetry'))
                Sizer1.Add((-1,5))
                Sizer = wx.FlexGridSizer(0,3,0,0)
                Sizer1.Add(Sizer)
                for var in symHolds:
                    varMean = G2obj.fmtVarDescr(var)
                    helptext = "Prevents variable: "+ var + " ("+ varMean + ")\nfrom being changed"
                    ch = G2G.HelpButton(panel,helptext)
                    Sizer.Add(ch)
                    Sizer.Add(wx.StaticText(panel,label='FIXED'),0,WACV|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT,2)
                    Sizer.Add(wx.StaticText(panel,label=var),0,WACV|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT)
            else:
                Sizer1.Add(wx.StaticText(panel,label='No holds generated'))
            Sizer1.Add((-1,10))
            symGen,SymErr,SymHelp = G2mv.GetSymEquiv(seqmode,seqhistnum)
            symGenD,SymErrD,SymHelpD = G2mv.GetDroppedSym(seqmode,seqhistnum)
            
            if len(symGen) == 0:
                Sizer1.Add(wx.StaticText(panel,label='No equivalences generated'))
                return Sizer1
            Sizer1.Add(wx.StaticText(panel,label='Equivalences generated based on cell/space group input'))
            Sizer1.Add((-1,5))
            Sizer = wx.FlexGridSizer(0,5,0,0)
            Sizer1.Add(Sizer)
            helptext = ''
            for sym,(warnmsg,note),helptext in zip(symGen,SymErr,SymHelp):
                if warnmsg:
                    if helptext: helptext += '\n\n'
                    helptext += warnmsg
                if helptext:
                    ch = G2G.HelpButton(panel,helptext)
                    Sizer.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
                else:
                    Sizer.Add((-1,-1))
                Sizer.Add(wx.StaticText(panel,wx.ID_ANY,'EQUIV'),
                    0,WACV|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT,2)
                Sizer.Add(wx.StaticText(panel,wx.ID_ANY,sym),0,WACV|wx.ALIGN_LEFT|wx.RIGHT|wx.LEFT,2)
                if note:
                    Sizer.Add(wx.StaticText(panel,wx.ID_ANY,note),0,WACV|wx.ALIGN_LEFT|wx.RIGHT|wx.LEFT,2)
                else:
                    Sizer.Add((-1,-1))
                Sizer.Add((-1,-1))
            for sym,(warnmsg,note),helptext in zip(symGenD,SymErrD,SymHelpD):
                if warnmsg:
                    if helptext: helptext += '\n\n'
                    helptext += warnmsg
                if helptext:
                    ch = G2G.HelpButton(panel,helptext)
                    Sizer.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
                else:
                    Sizer.Add((-1,-1))
                Sizer.Add(wx.StaticText(panel,wx.ID_ANY,'EQUIV'),
                    0,WACV|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT,2)
                Sizer.Add(wx.StaticText(panel,wx.ID_ANY,sym),0,WACV|wx.ALIGN_LEFT|wx.RIGHT|wx.LEFT,2)
                if note:
                    Sizer.Add(wx.StaticText(panel,wx.ID_ANY,note),0,WACV|wx.ALIGN_LEFT|wx.RIGHT|wx.LEFT,2)
                else:
                    Sizer.Add((-1,-1))
                Sizer.Add((-1,-1))
            return Sizer1
        constSizer = wx.FlexGridSizer(0,8,0,0)
        maxlen = 50 # characters before wrapping a constraint
        for Id,item in enumerate(data[name]):
            refineflag = False
            helptext = ""
            eqString = ['',]
            problemItem, warnmsg, note = G2mv.getConstrError(item,seqmode,seqhistnum)
            #badVar = False
            #for term in item[:-3]:
            #    if str(term[1]) in G2mv.problemVars:
            #        problemItem = True
            if item[-1] == 'h': # Hold on variable
                constSizer.Add((-1,-1),0)              # blank space for edit button
                typeString = 'FIXED'
                #var = str(item[0][1])
                var,explain,note,warnmsg = item[0][1].fmtVarByMode(seqmode,note,warnmsg)
                #if '?' in var: badVar = True
                varMean = G2obj.fmtVarDescr(var)
                eqString[-1] =  var +'   '
                helptext = "Prevents variable:\n"+ var + " ("+ varMean + ")\nfrom being changed"
            elif item[-1] == 'f' or item[-1] == 'e' or item[-1] == 'c': # not true on original-style (2011?) constraints
                constEdit = wx.Button(panel,label='Edit',style=wx.BU_EXACTFIT)
                constEdit.Bind(wx.EVT_BUTTON,OnConstEdit)
                constSizer.Add(constEdit)            # edit button
                Indx[constEdit.GetId()] = [Id,name]
                if item[-1] == 'f':
                    helptext = "A new variable"
                    if item[-3]:
                        helptext += " named "+str(item[-3])
                    helptext += " is created from a linear combination of the following variables:\n"
                    for term in item[:-3]:
                        m = term[0]
                        if np.isclose(m,0): continue
                        #var = str(term[1])
                        var,explain,note,warnmsg = term[1].fmtVarByMode(seqmode,note,warnmsg)
                        #if '?' in var: badVar = True
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        if eqString[-1] != '':
                            if m >= 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                                m = abs(m)
                        if m == 1:
                            eqString[-1] += '{:} '.format(var)
                        else:
                            eqString[-1] += '{:.3g}*{:} '.format(m,var)
                        varMean = G2obj.fmtVarDescr(var)
                        helptext += '\n  {:.5g} * {:} '.format(m,var) + " ("+ varMean + ")"
                    # Add extra notes about this constraint (such as from ISODISTORT)
                    if '_Explain' in data:
                        hlptxt = None
                        try:
                            hlptxt = data['_Explain'].get(item[-3])
                        except TypeError:
                            # Patch fixed Nov 2021. Older projects have phase RanId in 
                            # item[-3].phase rather than a properly formed G2VarObj
                            hlptxt = data['_Explain'].get(str(item[-3].phase)+item[-3].name)
                        if hlptxt:
                            helptext += '\n\n'+ hlptxt
                    if item[-3]:
                        typeString = str(item[-3]) + ' ='
                        if note: note += ', '
                        note += '(NEW VAR)'
                    else:
                        typeString = 'New Variable = '
                    #print 'refine',item[-2]
                    refineflag = True
                elif item[-1] == 'c':
                    helptext = "The following variables are constrained to equal a constant:"
                    for term in item[:-3]:
                        #var = str(term[1])
                        var,explain,note,warnmsg = term[1].fmtVarByMode(seqmode,note,warnmsg)
                        #if '?' in var: badVar = True
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        m = term[0]
                        if eqString[-1] != '':
                            if term[0] > 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                                m = -term[0]
                        if m == 1:
                            eqString[-1] += '{:} '.format(var)
                        else:
                            eqString[-1] += '{:.3g}*{:} '.format(m,var)
                        varMean = G2obj.fmtVarDescr(var)
                        helptext += '\n  {:.5g} * {:} '.format(m,var) + " ("+ varMean + ")"
                        helptext += explain
                    typeString = 'CONST '
                    eqString[-1] += ' = '+str(item[-3])
                elif item[-1] == 'e' and len(item[:-3]) == 2:
                    if item[0][0] == 0: item[0][0] = 1.0
                    if item[1][0] == 0: item[1][0] = 1.0
                    #var = str(item[0][1])
                    var,explain,note,warnmsg = item[0][1].fmtVarByMode(seqmode,note,warnmsg)
                    #if '?' in var: badVar = True
                    helptext = 'Variable {:} '.format(var) + " ("+ G2obj.fmtVarDescr(var) + ")"
                    helptext += "\n\nis equivalent to "
                    m = item[0][0]/item[1][0]
                    #var1 = str(item[1][1])
                    var1,explain,note,warnmsg = item[1][1].fmtVarByMode(seqmode,note,warnmsg)
                    helptext += '\n  {:.5g} * {:} '.format(m,var1) + " ("+ G2obj.fmtVarDescr(var1) + ")"
                    eqString[-1] += '{:} = {:}'.format(var1,var)
                    if m != 1:
                        eqString[-1] += ' / ' + str(m)
                    typeString = 'EQUIV '
                elif item[-1] == 'e':
                    helptext = "The following variable:"
                    normval = item[0][0]
                    indepterm = item[0][1]
                    for i,term in enumerate(item[:-3]):
                        #var = str(term[1])
                        var,explain,note,warnmsg = term[1].fmtVarByMode(seqmode,note,warnmsg)
                        #if '?' in var: badVar = True
                        if term[0] == 0: term[0] = 1.0
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        varMean = G2obj.fmtVarDescr(var)
                        if i == 0: # move independent variable to end, as requested by Bob
                            helptext += '\n{:} '.format(var) + " ("+ varMean + ")"
                            helptext += "\n\nis equivalent to the following, noting multipliers:"
                            continue
                        elif eqString[-1] != '':
                            eqString[-1] += ' = '
                        #m = normval/term[0]
                        m = term[0]/normval
                        if m == 1:
                            eqString[-1] += '{:}'.format(var)
                        else:
                            eqString[-1] += '{:.3g}*{:} '.format(m,var)
                        helptext += '\n  {:.5g} * {:} '.format(m,var) + " ("+ varMean + ")"
                    eqString[-1] += ' = {:} '.format(indepterm)
                    typeString = 'EQUIV '
                else:
                    print ('Unexpected constraint'+item)
                
            else:
                print ('Removing old-style constraints')
                data[name] = []
                return constSizer
            if warnmsg:
                if helptext: helptext += '\n\nNote warning:\n'
                helptext += warnmsg
            if helptext:
                ch = G2G.HelpButton(panel,helptext)
                constSizer.Add(ch)
            else:
                constSizer.Add((-1,-1))
            constDel = wx.CheckBox(panel,label='sel ')
            constSizer.Add(constDel) # delete selection
            panel.delBtn.checkboxList.append([constDel,Id,name])
            if refineflag:
                refresh = lambda event: wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True)
                ch = G2G.G2CheckBox(panel,'vary ',item,-2,OnChange=refresh)
                constSizer.Add(ch)
            else:
                constSizer.Add((-1,-1))
            if typeString.strip().endswith('='):
                constSizer.Add(wx.StaticText(panel,label=typeString),0,wx.EXPAND,1)
            else:
                constSizer.Add(wx.StaticText(panel,label=typeString),0,wx.EXPAND,1)
            #if badVar: eqString[-1] += ' -- Error: variable removed'
            #if note: eqString[-1] += '  (' + note + ')'
            if len(eqString) > 1:
                Eq = wx.BoxSizer(wx.VERTICAL)
                for s in eqString:
                    line = wx.StaticText(panel,label=s)
                    if problemItem:
                        line.SetForegroundColour(wx.BLACK)
                        line.SetBackgroundColour(wx.YELLOW)
                    Eq.Add(line,0)
                Eq.Add((-1,4))
            else:
                Eq = wx.StaticText(panel,label=eqString[0])
                if problemItem:
                    Eq.SetForegroundColour(wx.BLACK)
                    Eq.SetBackgroundColour(wx.YELLOW)
            constSizer.Add(Eq)
            constSizer.Add((3,3))
            if note:
                Eq = wx.StaticText(panel,label=note,style=WACV)
                if problemItem:
                    Eq.SetForegroundColour(wx.BLACK)
                    Eq.SetBackgroundColour(wx.YELLOW)
            else:
                Eq = (-1,-1)
            constSizer.Add(Eq,1,wx.EXPAND,3)
        if panel.delBtn.checkboxList:
            panel.delBtn.Enable(True)
        else:
            panel.delBtn.Enable(False)
        return constSizer

    def OnConstDel(event):
        'Delete selected constraints'
        sel = G2frame.constr.GetSelection()
        selList = event.GetEventObject().checkboxList
        selList.reverse()
        for obj,Id,name in event.GetEventObject().checkboxList:
            if obj.GetValue(): del data[name][Id]
        wx.CallAfter(UpdateConstraints,G2frame,data,sel,True)

    def OnConstEdit(event):
        '''Called to edit an individual contraint in response to a
        click on its Edit button
        '''
        Obj = event.GetEventObject()
        sel = G2frame.constr.GetSelection()
        Id,name = Indx[Obj.GetId()]
        if data[name][Id][-1] == 'f':
            items = data[name][Id][:-3]
            constType = 'New Variable'
            if data[name][Id][-3]:
                varname = str(data[name][Id][-3])
            else:
                varname = ""
            lbl = 'Enter multiplier for each parameter in the New Var expression'
            dlg = ConstraintDialog(G2frame,constType,lbl,items,
                varname=varname,varyflag=data[name][Id][-2])
        elif data[name][Id][-1] == 'c':
            items = data[name][Id][:-3]+[
                [str(data[name][Id][-3]),'fixed value =']]
            constType = 'Constraint'
            lbl = 'Edit value for each term in constant constraint sum'
            dlg = ConstraintDialog(G2frame,constType,lbl,items)
        elif data[name][Id][-1] == 'e':
            items = data[name][Id][:-3]
            constType = 'Equivalence'
            lbl = 'The following terms are set to be equal:'
            dlg = ConstraintDialog(G2frame,constType,lbl,items,'*')
        else:
            return
        try:
            prev = copy.deepcopy(data[name][Id])
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                for i in range(len(data[name][Id][:-3])):
                    if type(data[name][Id][i]) is tuple: # fix non-mutable construct
                        data[name][Id][i] = list(data[name][Id][i])
                    data[name][Id][i][0] = result[i][0]
                if data[name][Id][-1] == 'c':
                    data[name][Id][-3] = str(result[-1][0])
                elif data[name][Id][-1] == 'f':
                    data[name][Id][-2] = dlg.newvar[1]
                    if dlg.newvar[0]:
                        # process the variable name to put in global form (::var)
                        varname = str(dlg.newvar[0]).strip().replace(' ','_')
                        if varname.startswith('::'):
                            varname = varname[2:]
                        varname = varname.replace(':',';')
                        if varname:
                            data[name][Id][-3] = varname
                        else:
                            data[name][Id][-3] = ''
                if not CheckChangedConstraint():
                    data[name][Id] = prev
            else:
                data[name][Id] = prev
        except:
            import traceback
            print (traceback.format_exc())
        finally:
            dlg.Destroy()
#        wx.CallAfter(OnPageChanged,None)
        G2frame.dataWindow.ClearData() 
        wx.CallAfter(UpdateConstraints,G2frame,data,sel,True)
   
    def UpdateConstraintPanel(panel,typ):
        '''Update the contents of the selected Constraint
        notebook tab. Called in :func:`OnPageChanged`
        '''
        if panel.GetSizer(): panel.GetSizer().Clear(True)
        Siz = wx.BoxSizer(wx.VERTICAL)
        Siz.Add((5,5),0)
        if typ != 'Sym-Generated':
            butSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(panel, wx.ID_ANY, 'Show Errors')
            btn.Bind(wx.EVT_BUTTON,lambda event: G2G.ShowScrolledInfo(panel,errmsg,header='Error info'))
            butSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
            btn.Enable(len(errmsg) > 0)
            btn = wx.Button(panel, wx.ID_ANY, 'Show Warnings')
            butSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
            btn.Bind(wx.EVT_BUTTON,lambda event: G2G.ShowScrolledInfo(panel,warnmsg.replace('&','&&')))
            btn.Enable(len(warnmsg) > 0)
            btn = wx.Button(panel, wx.ID_ANY, 'Show generated constraints')
            butSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
            txt = G2mv.VarRemapShow(linelen=999).replace('&','&&')
            btn.Bind(wx.EVT_BUTTON,lambda event:
                         G2G.ShowScrolledColText(panel,
                        '*** Constraints after processing ***'+txt,
                         header='Generated constraints',col1len=80))
            panel.delBtn = wx.Button(panel, wx.ID_ANY, 'Delete selected')
            butSizer.Add(panel.delBtn,0,wx.ALIGN_CENTER_VERTICAL)
            panel.delBtn.Bind(wx.EVT_BUTTON,OnConstDel)
            panel.delBtn.checkboxList = []
            butSizer.Add((-1,-1),1,wx.EXPAND,1)
            butSizer.Add(G2G.HelpButton(panel,helpIndex='Constraints'))
            Siz.Add(butSizer,0,wx.EXPAND)
            if G2frame.testSeqRefineMode():
                butSizer = wx.BoxSizer(wx.HORIZONTAL)
                butSizer.Add(wx.StaticText(panel,wx.ID_ANY,'  Sequential Ref. Settings.  Wildcard use: '),0,WACV)
                btn = G2G.EnumSelector(panel, data, '_seqmode',
                        ['Set hist # to *', 'Ignore unless hist=*', 'Use as supplied'],
                        ['auto-wildcard',   'wildcards-only',       'use-all'],
                        lambda x: wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True))
                butSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
                butSizer.Add(G2G.HelpButton(panel,helpIndex='Constraints-SeqRef'))
                butSizer.Add(wx.StaticText(panel,wx.ID_ANY,'  Selected histogram: '),0,WACV)
                btn = G2G.EnumSelector(panel, data, '_seqhist',
                        list(seqHistList),list(range(len(seqHistList))),
                        lambda x: wx.CallAfter(UpdateConstraints, G2frame, data, G2frame.constr.GetSelection(), True))
                butSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
                Siz.Add(butSizer,0)
            G2G.HorizontalLine(Siz,panel)
#            Siz.Add((5,5),0)
        Siz.Add(MakeConstraintsSizer(typ,panel),1,wx.EXPAND)
        panel.SetSizer(Siz,True)
        Size = Siz.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        panel.SetSize(Size)
        panel.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
        panel.Show()

    def OnPageChanged(event,selectTab=None):
        '''Called when a tab is pressed or when a "select tab" menu button is
        used (see RaisePage), or to refresh the current tab contents (event=None)
        '''
        if event:       #page change event!
            page = event.GetSelection()
        elif selectTab: #reload previous
            page = selectTab
        else: # called directly, get current page
            try:
                page = G2frame.constr.GetSelection()
            except:
                if GSASIIpath.GetConfigValue('debug'): print('DBG_gpx open error:C++ Run time error - skipped')
                return
        G2frame.constr.ChangeSelection(page)
        text = G2frame.constr.GetPageText(page)
        G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_EQUIVALANCEATOMS,False)
#        G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_ADDRIDING,False)
        if text == 'Histogram/Phase':
            enableEditCons = [False]+4*[True]
            G2frame.Page = [page,'hap']
            UpdateConstraintPanel(HAPConstr,'HAP')
        elif text == 'Histogram':
            enableEditCons = [False]+4*[True]
            G2frame.Page = [page,'hst']
            UpdateConstraintPanel(HistConstr,'Hist')
        elif text == 'Phase':
            enableEditCons = 5*[True]
            G2frame.Page = [page,'phs']
            G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_EQUIVALANCEATOMS,True)
#            G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_ADDRIDING,True)
            if 'DELETED' in str(PhaseConstr):   #seems to be no other way to do this (wx bug)
                if GSASIIpath.GetConfigValue('debug'):
                    print ('DBG_wx error: PhaseConstr not cleanly deleted after Refine')
                return
            UpdateConstraintPanel(PhaseConstr,'Phase')
        elif text == 'Global':
            enableEditCons = [False]+4*[True]
            G2frame.Page = [page,'glb']
            UpdateConstraintPanel(GlobalConstr,'Global')
        else:
            enableEditCons = 5*[False]
            G2frame.Page = [page,'sym']
            UpdateConstraintPanel(SymConstr,'Sym-Generated')
        # remove menu items when not allowed
        for obj,flag in zip(G2frame.dataWindow.ConstraintEdit.GetMenuItems(),enableEditCons): 
            obj.Enable(flag)
        G2frame.dataWindow.SetDataSize()

    def RaisePage(event):
        'Respond to a "select tab" menu button'
        try:
            i = (G2G.wxID_CONSPHASE,
                 G2G.wxID_CONSHAP,
                 G2G.wxID_CONSHIST,
                 G2G.wxID_CONSGLOBAL,
                 G2G.wxID_CONSSYM,
                ).index(event.GetId())
            G2frame.constr.SetSelection(i)
            wx.CallAfter(OnPageChanged,None)
        except ValueError:
            print('Unexpected event in RaisePage')

    def SetStatusLine(text):
        G2frame.GetStatusBar().SetStatusText(text,1)
        
    def OnShowISODISTORT(event):
        ShowIsoDistortCalc(G2frame)

    #### UpdateConstraints execution starts here ##############################
    if Clear:
        G2frame.dataWindow.ClearData() 
    if not data:  # usually created in CheckNotebook
        data.update({'Hist':[],'HAP':[],'Phase':[],'Global':[],
                         '_seqmode':'auto-wildcard', '_seqhist':0})       #empty dict - fill it
    if 'Global' not in data:                                            #patch
        data['Global'] = []
    seqHistList = G2frame.testSeqRefineMode()
    # patch: added ~version 5030 -- new mode for wild-card use in seq refs
    # note that default for older sequential fits is 'wildcards-only' but in new GPX is 'auto-wildcard'
    if seqHistList:
        data['_seqmode'] = data.get('_seqmode','wildcards-only')
    else:
        data['_seqmode'] = data.get('_seqmode','auto-wildcard')
    data['_seqhist'] = data.get('_seqhist',0)
    # end patch
    
    # DEBUG code #=====================================
    #import GSASIIconstrGUI
    #reload(GSASIIconstrGUI)
    #reload(G2obj)
    #reload(G2stIO)
    #import GSASIIstrMain
    #reload(GSASIIstrMain)    
    #reload(G2mv)
    #reload(G2gd)
    #===================================================
    seqmode = 'use-all'
    seqhistnum = None
    if seqHistList:  # Selections used with sequential refinements
        seqmode = data.get('_seqmode','wildcards-only')
        seqhistnum = min(data.get('_seqhist',0),len(seqHistList)-1)
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_SHOWISO,True)
#removed this check as it prevents examination of ISODISTORT constraints without data
    # if not len(Phases) or not len(Histograms):        
    #     dlg = wx.MessageDialog(G2frame,'You need both phases and histograms to see Constraints',
    #         'No phases or histograms')
    #     dlg.CenterOnParent()
    #     dlg.ShowModal()
    #     dlg.Destroy()
    #     return
    # for p in Phases:
    #     if 'ISODISTORT' in Phases[p] and 'G2VarList' in Phases[p]['ISODISTORT']:
    #         G2frame.dataWindow.ConstraintEdit.Enable(G2G.wxID_SHOWISO,True)
    #         break
    ###### patch: convert old-style (str) variables in constraints to G2VarObj objects #####
    for key,value in data.items():
        if key.startswith('_'): continue
        j = 0
        for cons in value:
            #print cons             # DEBUG
            for i in range(len(cons[:-3])):
                if type(cons[i][1]) is str:
                    cons[i][1] = G2obj.G2VarObj(cons[i][1])
                    j += 1
        if j:
            print (str(key) + ': '+str(j)+' variable(s) as strings converted to objects')
    ##### end patch #############################################################################
    rigidbodyDict = G2frame.GPXtree.GetItemPyData(
        G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[],'Spin':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
    badPhaseParms = ['Ax','Ay','Az','Amul','AI/A','Atype','SHorder','AwaveType','FwaveType','PwaveType','MwaveType','Vol','isMag',]
    globalList = list(rbDict.keys())
    globalList.sort()
    try:
        AtomDict = dict([Phases[phase]['pId'],Phases[phase]['Atoms']] for phase in Phases)
    except KeyError:
        G2frame.ErrorDialog('Constraint Error','Constraints cannot be set until a cycle of least squares'+
                            ' has been run.\nWe suggest you refine a scale factor.')
        return

    # create a list of the phase variables
    symHolds = []
    (Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,EFtable,ORBtables,BLtable,MFtable,maxSSwave) = \
        G2stIO.GetPhaseData(Phases,rbIds=rbIds,Print=False,symHold=symHolds)
    phaseList = []
    for item in phaseDict:
        if item.split(':')[2] not in badPhaseParms:
            phaseList.append(item)
    phaseList.sort()
    phaseAtNames = {}
    phaseAtTypes = {}
    TypeList = []
    for item in phaseList:
        Split = item.split(':')
        if Split[2][:2] in ['AU','Af','dA','AM']:
            Id = int(Split[0])
            phaseAtNames[item] = AtomDict[Id][int(Split[3])][0]
            phaseAtTypes[item] = AtomDict[Id][int(Split[3])][1]
            if phaseAtTypes[item] not in TypeList:
                TypeList.append(phaseAtTypes[item])
        else:
            phaseAtNames[item] = ''
            phaseAtTypes[item] = ''
             
    # create a list of the hist*phase variables
    if seqHistList: # for sequential refinement, only process selected histgram in list
        histDict = {seqHistList[seqhistnum]:Histograms[seqHistList[seqhistnum]]}
    else:
        histDict = Histograms
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,histDict,Print=False,resetRefList=False)
    hapList = sorted([i for i in hapDict.keys() if i.split(':')[2] not in ('Type',)])
    # derive list of variables 
    if seqHistList: # convert histogram # to wildcard
        wildList = [] # list of variables with "*" for histogram number
        for i in hapList:
            s = i.split(':')
            if s[1] == "": continue
            s[1] = '*'
            sj = ':'.join(s)
            if sj not in wildList: wildList.append(sj)
        if seqmode == 'use-all':
            hapList += wildList        
        else:
            hapList = wildList        
    histVary,histDict,controlDict = G2stIO.GetHistogramData(histDict,Print=False)
    histList = list(histDict.keys())
    histList.sort()
    if seqHistList: # convert histogram # to wildcard
        wildList = [] # list of variables with "*" for histogram number
        for i in histList:
            s = i.split(':')
            if s[1] == "": continue
            s[1] = '*'
            sj = ':'.join(s)
            if sj not in wildList: wildList.append(sj)
        if seqmode == 'use-all':
            histList += wildList        
        else:
            histList = wildList        

    Indx = {}
    G2frame.Page = [0,'phs']
    
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.ConstraintMenu)
    SetStatusLine('')
    
    G2frame.Bind(wx.EVT_MENU, OnAddConstraint, id=G2G.wxID_CONSTRAINTADD)
    G2frame.Bind(wx.EVT_MENU, OnAddFunction, id=G2G.wxID_FUNCTADD)
    G2frame.Bind(wx.EVT_MENU, OnAddEquivalence, id=G2G.wxID_EQUIVADD)
    G2frame.Bind(wx.EVT_MENU, OnAddHold, id=G2G.wxID_HOLDADD)
    G2frame.Bind(wx.EVT_MENU, OnAddAtomEquiv, id=G2G.wxID_EQUIVALANCEATOMS)
#    G2frame.Bind(wx.EVT_MENU, OnAddRiding, id=G2G.wxID_ADDRIDING)
    G2frame.Bind(wx.EVT_MENU, OnShowISODISTORT, id=G2G.wxID_SHOWISO)
    # tab commands
    for id in (G2G.wxID_CONSPHASE,
               G2G.wxID_CONSHAP,
               G2G.wxID_CONSHIST,
               G2G.wxID_CONSGLOBAL,
               G2G.wxID_CONSSYM,
               ):
        G2frame.Bind(wx.EVT_MENU, RaisePage,id=id)

    #G2frame.constr = G2G.GSNoteBook(parent=G2frame.dataWindow,size=G2frame.dataWindow.GetClientSize())
    G2frame.constr = G2G.GSNoteBook(parent=G2frame.dataWindow)
    G2frame.dataWindow.GetSizer().Add(G2frame.constr,1,wx.ALL|wx.EXPAND)
    # note that order of pages is hard-coded in RaisePage
    PhaseConstr = wx.ScrolledWindow(G2frame.constr)
    G2frame.constr.AddPage(PhaseConstr,'Phase')
    HAPConstr = wx.ScrolledWindow(G2frame.constr)
    G2frame.constr.AddPage(HAPConstr,'Histogram/Phase')
    HistConstr = wx.ScrolledWindow(G2frame.constr)
    G2frame.constr.AddPage(HistConstr,'Histogram')
    GlobalConstr = wx.ScrolledWindow(G2frame.constr)
    G2frame.constr.AddPage(GlobalConstr,'Global')
    SymConstr = wx.ScrolledWindow(G2frame.constr)
    G2frame.constr.AddPage(SymConstr,'Sym-Generated')
    wx.CallAfter(OnPageChanged,None,selectTab)
    G2frame.constr.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    errmsg,warnmsg = WarnConstraintLimit()  # check limits & constraints
    if errmsg:
        G2frame.ErrorDialog('Constraint Error',
            'Error in constraints.\nCheck console output for more information'+
            ' or press "Show Errors" & "Show Warnings" buttons',
            parent=G2frame)
        if seqhistnum is None:
            print ('\nError message(s):\n',errmsg)
        else:
            print ('\nError message(s) for histogram #{}:\n{}'.format(seqhistnum,errmsg))
        if warnmsg: print ('\nAlso note these constraint warning(s):\n'+warnmsg)
    #elif GSASIIpath.GetConfigValue('debug'):
    #    print ('Generated constraints\n',G2mv.VarRemapShow())
        
###### check scale & phase fractions, create constraint if needed #############
def CheckAllScalePhaseFractions(G2frame,refine=True):
    '''Check if scale factor and all phase fractions are refined without a constraint
    for all used histograms, if so, offer the user a chance to create a constraint
    on the sum of phase fractions

    :returns: False if refinement should be continued
    '''
    histograms, phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    cId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints') 
    Constraints = G2frame.GPXtree.GetItemPyData(cId)

    problems = []
    for h,hist in enumerate(histograms):
        if CheckScalePhaseFractions(G2frame,hist,histograms,phases,Constraints):
            problems.append((h,hist))
    if len(problems) == 0: return
    msg = 'You are refining the scale factor and all phase fractions for histogram(s) #'
    for i,(h,hist) in enumerate(problems):
        if i: msg += ', '
        msg += str(h)
    msg += '. This is not recommended as it will produce an unstable refinement. Do you want to create constrain(s) on the sum of phase fractions to address this?'
    if refine:
        msg += '\n\nRefinement may cause a significant shift in scaling, so it may be best to refine only a few other parameters (Press "No" to continue refinement without, "Cancel" to stop.)'
        opts = wx.YES|wx.NO|wx.CANCEL
    else:
        opts = wx.YES|wx.NO
    dlg = wx.MessageDialog(G2frame,msg,'Warning: Constraint Needed',opts)
    ans = dlg.ShowModal()
    dlg.Destroy()
    if ans == wx.ID_YES:
        for h,hist in problems:
            constr = []
            for p in phases:
                if hist not in phases[p]['Histograms']: continue
                if not phases[p]['Histograms'][hist]['Use']: continue
                constr += [[1.0,G2obj.G2VarObj(':'.join(
                    (str(phases[p]['pId']),
                    str(histograms[hist]['hId']),
                    'Scale')
                    ))]]
            Constraints['HAP'].append(constr+[1.0,None,'c'])
        wx.CallAfter(G2frame.GPXtree.SelectItem,cId) # should call SelectDataTreeItem
        UpdateConstraints(G2frame,Constraints,1,True) # repaint with HAP tab
        return False
    elif ans == wx.ID_NO:
        return False
    return True

def CheckScalePhaseFractions(G2frame,hist,histograms,phases,Constraints):
    '''Check if scale factor and all phase fractions are refined without a constraint
    for histogram hist, if so, offer the user a chance to create a constraint
    on the sum of phase fractions
    '''
    if G2frame.testSeqRefineMode():  # more work needed due to seqmode
        return False
#        histStr = '*'
    else: 
        histStr = str(histograms[hist]['hId'])
    # Is this powder? 
    if not hist.startswith('PWDR '): return False
    # do this only if the scale factor is varied
    if not histograms[hist]['Sample Parameters']['Scale'][1]: return False
    # are all phase fractions varied in all used histograms?
    phaseCount = 0
    for p in phases:
        if hist not in phases[p]['Histograms']: continue
        if phases[p]['Histograms'][hist]['Use'] and not phases[p]['Histograms'][hist]['Scale'][1]:
            return False
        else:
            phaseCount += 1
    
    # all phase fractions and scale factor varied, now scan for a constraint
    for c in Constraints.get('HAP',[]):
        if c[-1] != 'c': continue
        if not c[-3]: continue
        if len(c[:-3]) != phaseCount: continue
        # got a constraint equation with right number of terms, is it on phase fractions for
        # the correct histogram?
        if all([(i[1].name == 'Scale' and i[1].varname().split(':')[1] == histStr) for i in c[:-3]]):
            # got a constraint, this is OK
            return False
    return True
        
#### Make nuclear/magnetic phase transition constraints - called by OnTransform in G2phsGUI ##########
def TransConstraints(G2frame,oldPhase,newPhase,Trans,Vec,atCodes):
    '''Add constraints for new magnetic phase created via transformation of old
    nuclear one
    NB: A = [G11,G22,G33,2*G12,2*G13,2*G23]
    '''
    
    def SetUniqAj(pId,iA,SGData):
        SGLaue = SGData['SGLaue']
        SGUniq = SGData['SGUniq']
        if SGLaue in ['m3','m3m']:
            if iA in [0,1,2]:
                parm = '%d::%s'%(pId,'A0')
            else:
                parm = None
        elif SGLaue in ['4/m','4/mmm']:
            if iA in [0,1]:
                parm = '%d::%s'%(pId,'A0')
            elif iA == 2:
                parm = '%d::%s'%(pId,'A2')
            else:
                parm = None
        elif SGLaue in ['6/m','6/mmm','3m1', '31m', '3']:
            if iA in [0,1,3]:
                parm = '%d::%s'%(pId,'A0')
            elif iA == 2:
                parm = '%d::%s'%(pId,'A2')
            else:
                parm = None
        elif SGLaue in ['3R', '3mR']:
            if ia in [0,1,2]:
                parm = '%d::%s'%(pId,'A0')
            else:
                parm = '%d::%s'%(pId,'A3')
        elif SGLaue in ['mmm',]:
            if iA in [0,1,2]:
                parm = '%d::A%s'%(pId,iA)
            else:
                parm = None
        elif SGLaue == '2/m':
            if iA in [0,1,2]:
                parm = '%d::A%s'%(pId,iA)
            elif iA == 3 and SGUniq == 'c':
                parm = '%d::A%s'%(pId,iA)
            elif iA == 4 and SGUniq == 'b':
                parm = '%d::A%s'%(pId,iA)
            elif iA == 5 and SGUniq == 'a':
                parm = '%d::A%s'%(pId,iA)
            else:
                parm = None            
        else:
            parm = '%d::A%s'%(pId,iA)
        return parm
    
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    UseList = newPhase['Histograms']
    detTrans = np.abs(nl.det(Trans))
    opId = oldPhase['pId']
    npId = newPhase['pId']
    cx,ct,cs,cia = newPhase['General']['AtomPtrs']
    nAtoms = newPhase['Atoms']
    nSGData = newPhase['General']['SGData']
    #oAcof = G2lat.cell2A(oldPhase['General']['Cell'][1:7])
    #nAcof = G2lat.cell2A(newPhase['General']['Cell'][1:7])
    item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints')
    if not item:
        print('Error: no constraints in Data Tree')
        return
    constraints = G2frame.GPXtree.GetItemPyData(item)
    xnames = ['dAx','dAy','dAz']
    # constraints on matching atom params between phases
    for ia,code in enumerate(atCodes):
        atom = nAtoms[ia]
        if not ia and atom[cia] == 'A':
            wx.MessageDialog(G2frame,
                'Anisotropic thermal motion constraints are not developed at the present time',
                'Anisotropic thermal constraint?',style=wx.ICON_INFORMATION).ShowModal()
        siteSym = G2spc.SytSym(atom[cx:cx+3],nSGData)[0]
        CSX = G2spc.GetCSxinel(siteSym)
#        CSU = G2spc.GetCSuinel(siteSym)
        item = code.split('+')[0]
        iat,opr = item.split(':')
        Nop = abs(int(opr))%100-1
        if '-' in opr:
            Nop *= -1
        Opr = oldPhase['General']['SGData']['SGOps'][abs(Nop)][0]
        if Nop < 0:         #inversion
            Opr *= -1
        XOpr = np.inner(Opr,Trans)
        invOpr = nl.inv(XOpr)
        for i,ix in enumerate(list(CSX[0])):
            if not ix:
                continue
            name = xnames[i]
            IndpCon = [1.0,G2obj.G2VarObj('%d::%s:%d'%(npId,name,ia))]
            DepCons = []
            for iop,opval in enumerate(invOpr[i]):
                if abs(opval) > 1e-6:
                    DepCons.append([opval,G2obj.G2VarObj('%d::%s:%s'%(opId,xnames[iop],iat))])
            if len(DepCons) == 1:
                constraints['Phase'].append([DepCons[0],IndpCon,None,None,'e'])
            elif len(DepCons) > 1:
                IndpCon[0] = -1.
                constraints['Phase'].append([IndpCon]+DepCons+[0.0,None,'c'])
        for name in ['Afrac','AUiso']:
            IndpCon = [1.0,G2obj.G2VarObj('%d::%s:%d'%(npId,name,ia))]
            DepCons = [1.0,G2obj.G2VarObj('%d::%s:%s'%(opId,name,iat))]
            constraints['Phase'].append([DepCons,IndpCon,None,None,'e'])
            
        # unfinished Anisotropic constraint generation
#        Uids = [[0,0,'AU11'],[1,1,'AU22'],[2,2,'AU33'],[0,1,'AU12'],[0,2,'AU13'],[1,2,'AU23']]
#        DepConsDict = dict(zip(Us,[[],[],[],[],[],[]]))
#        for iu,Uid in enumerate(Uids):
#            UMT = np.zeros((3,3))
#            UMT[Uid[0],Uid[1]] = 1
#            nUMT = G2lat.prodMGMT(UMT,invTrans)
#            nUT = G2lat.UijtoU6(nUMT)
#            for iu,nU in enumerate(nUT):
#                if abs(nU) > 1.e-8:
#                    parm = '%d::%s;%s'%(opId,Us[iu],iat)
#                    DepConsDict[Uid[2]].append([abs(nU%1.),G2obj.G2VarObj(parm)])
#        nUcof = atom[iu:iu+6]
#        conStrings = []
#        for iU,Usi in enumerate(Us):
#            parm = '%d::%s;%d'%(npId,Usi,ia)
#            parmDict[parm] = nUcof[iU]
#            varyList.append(parm)
#            IndpCon = [1.0,G2obj.G2VarObj(parm)]
#            conStr = str([IndpCon,DepConsDict[Usi]])
#            if conStr in conStrings:
#                continue
#            conStrings.append(conStr)
#            if len(DepConsDict[Usi]) == 1:
#                if DepConsDict[Usi][0]:
#                    constraints['Phase'].append([IndpCon,DepConsDict[Usi][0],None,None,'e'])
#            elif len(DepConsDict[Usi]) > 1:        
#                for Dep in DepConsDict[Usi]:
#                    Dep[0] *= -1
#                constraints['Phase'].append([IndpCon]+DepConsDict[Usi]+[0.0,None,'c'])
            
        #how do I do Uij's for most Trans?

    # constraints on lattice parameters between phases
    Aold = G2lat.cell2A(oldPhase['General']['Cell'][1:7])
    if True: # debug
        constraints['Phase'] += G2lat.GenCellConstraints(Trans,opId,npId,Aold,
                                oldPhase['General']['SGData'],nSGData,True)
        print('old A*',G2lat.cell2A(oldPhase['General']['Cell'][1:7]))
        print('new A*',G2lat.cell2A(newPhase['General']['Cell'][1:7]))
        print('old cell',oldPhase['General']['Cell'][1:7])
        print('new cell',newPhase['General']['Cell'][1:7])
    else:
        constraints['Phase'] += G2lat.GenCellConstraints(Trans,opId,npId,Aold,
                                oldPhase['General']['SGData'],nSGData,True)
    # constraints on HAP Scale, etc.
    for hId,hist in enumerate(UseList):    #HAP - seems OK
        ohapkey = '%d:%d:'%(opId,hId)
        nhapkey = '%d:%d:'%(npId,hId)
        IndpCon = [1.0,G2obj.G2VarObj(ohapkey+'Scale')]
        DepCons = [detTrans,G2obj.G2VarObj(nhapkey+'Scale')]
        constraints['HAP'].append([DepCons,IndpCon,None,None,'e'])
        for name in ['Size;i','Mustrain;i']:
            IndpCon = [1.0,G2obj.G2VarObj(ohapkey+name)]
            DepCons = [1.0,G2obj.G2VarObj(nhapkey+name)]
            constraints['HAP'].append([IndpCon,DepCons,None,None,'e'])
        
#### Rigid bodies #############################################################
resRBsel = None
def UpdateRigidBodies(G2frame,data):
    '''Called when Rigid bodies tree item is selected.
    Displays the rigid bodies in the data window
    '''  
    def OnPageChanged(event):
        global resList
        resList = []
        if event:       #page change event!
            page = event.GetSelection()
        else:
            try:
                page = G2frame.rbBook.GetSelection()
            except:
                if GSASIIpath.GetConfigValue('debug'): print('DBG_gpx open error:C++ Run time error - skipped')
                return
        G2frame.rbBook.ChangeSelection(page)
        text = G2frame.rbBook.GetPageText(page)
        if text == 'Vector rigid bodies':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.VectorBodyMenu)
            G2frame.Bind(wx.EVT_MENU, AddVectorRB, id=G2G.wxID_VECTORBODYADD)
            G2frame.Bind(wx.EVT_MENU, ExtractPhaseRB, id=G2G.wxID_VECTORBODYIMP)
            G2frame.Bind(wx.EVT_MENU, AddVectTrans, id=G2G.wxID_VECTORBODYEXTD)
            G2frame.Bind(wx.EVT_MENU, SaveVectorRB, id=G2G.wxID_VECTORBODYSAV)
            G2frame.Bind(wx.EVT_MENU, ReadVectorRB, id=G2G.wxID_VECTORBODYRD)
            G2frame.Page = [page,'vrb']
            UpdateVectorRB()
        elif text == 'Residue rigid bodies':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodyMenu)
            G2frame.Bind(wx.EVT_MENU, AddResidueRB, id=G2G.wxID_RIGIDBODYADD)
            G2frame.Bind(wx.EVT_MENU, ExtractPhaseRB, id=G2G.wxID_RIGIDBODYIMP)
            G2frame.Bind(wx.EVT_MENU, OnImportRigidBody, id=G2G.wxID_RIGIDBODYIMPORT)
            G2frame.Bind(wx.EVT_MENU, OnSaveRigidBody, id=G2G.wxID_RIGIDBODYSAVE)
            G2frame.Bind(wx.EVT_MENU, OnDefineTorsSeq, id=G2G.wxID_RESIDUETORSSEQ) #enable only if residue RBs exist?
            G2frame.Bind(wx.EVT_MENU, DumpVectorRB, id=G2G.wxID_RESBODYSAV)
            G2frame.Bind(wx.EVT_MENU, LoadVectorRB, id=G2G.wxID_RESBODYRD)
            G2frame.Page = [page,'rrb']
            UpdateResidueRB()
        elif text == 'Spinning rigid bodies':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SpinBodyMenu)
            G2frame.Bind(wx.EVT_MENU, AddSpinRB, id=G2G.wxID_SPINBODYADD)
            G2frame.Page = [page,'srb']
            UpdateSpinRB()
        else:
            G2gd.SetDataMenuBar(G2frame)
            #G2frame.Page = [page,'rrb']
            
    def getMacroFile(macName):
        defDir = os.path.join(GSASIIpath.path2GSAS2,'inputs','GSASIImacros')
        if not os.path.exists(defDir):  # patch 3/2024 for svn dir organization
            defDir = os.path.join(GSASIIpath.path2GSAS2,'GSASIImacros')
        if not os.path.exists(defDir):
            print('Warning: GSASIImacros directory not found')
            return []
        dlg = wx.FileDialog(G2frame,message='Choose '+macName+' rigid body macro file',
            defaultDir=defDir,defaultFile="",wildcard="GSAS-II macro file (*.mac)|*.mac",
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                macfile = dlg.GetPath()
                macro = open(macfile,'r')
                head = macro.readline()
                if macName not in head:
                    print (head)
                    print ('**** ERROR - wrong restraint macro file selected, try again ****')
                    macro = []
            else: # cancel was pressed
                macro = []
        finally:
            dlg.Destroy()
        return macro        #advanced past 1st line
        
    def getTextFile():
        dlg = wx.FileDialog(G2frame,'Choose rigid body text file', G2frame.LastGPXdir, '',
            "GSAS-II text file (*.txt)|*.txt|XYZ file (*.xyz)|*.xyz|"
            "Sybyl mol2 file (*.mol2)|*.mol2|PDB file (*.pdb;*.ent)|*.pdb;*.ent",
            wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                txtfile = dlg.GetPath()
                ext = os.path.splitext(txtfile)[1]
                text = open(txtfile,'r')
            else: # cancel was pressed
                ext = ''
                text = []
        finally:
            dlg.Destroy()
        if 'ent' in ext:
            ext = '.pdb'
        return text,ext.lower()
        
    def OnImportRigidBody(event):
        page = G2frame.rbBook.GetSelection()
        if 'Vector' in G2frame.rbBook.GetPageText(page):
            pass
        elif 'Residue' in G2frame.rbBook.GetPageText(page):
            try:
                ImportResidueRB()
            except Exception as msg:
                print('Error reading .xyz file\n  Error msg:',msg)
                if GSASIIpath.GetConfigValue('debug'): 
                    import traceback
                    print (traceback.format_exc())
                    
    def OnSaveRigidBody(event):
        page = G2frame.rbBook.GetSelection()
        if 'Vector' in G2frame.rbBook.GetPageText(page):
            pass
        elif 'Residue' in G2frame.rbBook.GetPageText(page):
            SaveResidueRB()
            
    def DumpVectorRB(event):
        global resRBsel
        if resRBsel not in data['Residue']:
            return
        rbData = data['Residue'][resRBsel]
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose file to save residue rigid body',
            pth, '', 'RRB files (*.resbody)|*.resbody',
            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.resbody'  # set extension
                fp = open(filename,'w')
                fp.write('Name: '+rbData['RBname']+'\n')
                fp.write('atNames: ')
                for i in rbData['atNames']:
                    fp.write(str(i)+" ") 
                fp.write('\n')
                for item in rbData['rbSeq']:
                    fp.write('rbSeq: ') 
                    fp.write('{:d} {:d} {:.1f}: '.format(*item[:3]))
                    for num in item[3]:
                        fp.write('{:d} '.format(num))
                    fp.write('\n')
                for i,sym in enumerate(rbData['rbTypes']):
                    fp.write("{:3s}".format(sym))
                    fp.write('{:8.5f}{:9.5f}{:9.5f}   '
                            .format(*rbData['rbXYZ'][i]))
                    fp.write('\n')
                fp.close()
                print ('Vector rigid body saved to: '+filename)
                
        finally:
            dlg.Destroy()
        
    def LoadVectorRB(event):
        AtInfo = data['Residue']['AtInfo']
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose file to read vector rigid body',
            pth, '', 'RRB files (*.resbody)|*.resbody',
            wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.resbody'  # set extension
                fp = open(filename,'r')
                l = fp.readline().strip()
                if 'Name' not in l:
                    fp.close()
                    G2frame.ErrorDialog('Read Error',
                        'File '+filename+' does not start with Name\nFirst line ='
                        +l+'\ninvalid file',parent=G2frame)
                    return
                name = l.split(':')[1].strip()
                line = fp.readline().strip().split(':')[1].split()
                atNames = [i for i in line]
                types = []
                coords = []
                l = fp.readline().strip()
                rbSeq = []
                while 'rbSeq' in l:
                    tag,vals,lst = l.split(':')
                    seq = []
                    for t,v in zip((int,int,float),vals.split()): 
                        seq.append(t(v))
                    seq.append([])
                    for num in lst.split():
                        seq[-1].append(int(num))
                    rbSeq.append(seq)
                    l = fp.readline().strip()                    
                while l:
                    nums = l.strip().split()
                    types.append(nums.pop(0))
                    t = types[-1]
                    if t not in AtInfo:
                        Info = G2elem.GetAtomInfo(t)
                        AtInfo[t] = [Info['Drad'],Info['Color']]
                    coords.append([float(nums.pop(0)) for j in range(3)])
                    l = fp.readline().strip()
                fp.close()
            else:
                return        
        finally:
            dlg.Destroy()
        coords = np.array(coords)
        rbid = ran.randint(0,sys.maxsize)
        namelist = [data['Residue'][key]['RBname'] for key in data['Residue']
                        if 'RBname' in data['Residue'][key]]
        name = G2obj.MakeUniqueLabel(name,namelist)
        data['Residue'][rbid] = {'RBname':name,
                'rbXYZ': coords,
                'rbRef':[0,1,2,False],
                'rbTypes':types, 'atNames':atNames,
                'useCount':0,
                'rbSeq':rbSeq, 'SelSeq':[0,0],}
        data['RBIds']['Residue'].append(rbid)
        UpdateResidueRB()
    
    def AddVectorRB(event):
        'Create a new vector rigid body'
        AtInfo = data['Vector']['AtInfo']
        dlg = G2G.MultiIntegerDialog(G2frame,'New Rigid Body',['No. atoms','No. translations'],[3,1])
        if dlg.ShowModal() == wx.ID_OK:
            nAtoms,nTrans = dlg.GetValues()
            if nAtoms < 3:
                dlg.Destroy()
                G2G.G2MessageBox(G2frame,'A vector rigid body must have 3 or more atoms')
                return
            rbid = ran.randint(0,sys.maxsize)
            vecMag = [1.0 for i in range(nTrans)]
            vecRef = [False for i in range(nTrans)]
            vecVal = [np.zeros((nAtoms,3)) for j in range(nTrans)]
            rbTypes = ['C' for i in range(nAtoms)]
            Info = G2elem.GetAtomInfo('C')
            AtInfo['C'] = [Info['Drad'],Info['Color']]
            name = 'UNKRB'
            namelist = [data['Vector'][key]['RBname'] for key in data['Vector']
                        if 'RBname' in data['Vector'][key]]
            name = G2obj.MakeUniqueLabel(name,namelist)
            data['Vector'][rbid] = {'RBname':name,'VectMag':vecMag,'rbXYZ':np.zeros((nAtoms,3)),
                'rbRef':[0,1,2,False],'VectRef':vecRef,'rbTypes':rbTypes,'rbVect':vecVal,'useCount':0}
            data['RBIds']['Vector'].append(rbid)
        dlg.Destroy()
        UpdateVectorRB()

    def ExtractPhaseRB(event):
        'Extract a rigid body from a file with a phase'
        def SetupDrawing(atmData):
            '''Add the dicts needed for G2plt.PlotStructure to work to the
            reader .Phase object
            '''
            generalData = atmData['General']
            generalData['BondRadii'] = []

            G2phG.SetDrawingDefaults(atmData['Drawing'])
            atmData['Drawing'].update(
                {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':150.,
                     'viewDir':[0,0,1],'atomPtrs': [2, 1, 6, 17],
                     })
            atmData['Drawing']['showRigidBodies'] = False
            generalData['Map'] = {'MapType':False, 'rho':[]}
            generalData['AtomTypes'] = []
            generalData['BondRadii'] = []
            generalData['AngleRadii'] = []
            generalData['vdWRadii'] = []
            generalData['Color'] = []
            generalData['Isotopes'] = {}
            generalData['Isotope'] = {}
            cx,ct,cs,cia = generalData['AtomPtrs']
            generalData['Mydir'] = G2frame.dirname
            for iat,atom in enumerate(atmData['Atoms']):
                atom[ct] = atom[ct].lower().capitalize()      #force elem symbol to standard form
                if atom[ct] not in generalData['AtomTypes'] and atom[ct] != 'UNK':
                    Info = G2elem.GetAtomInfo(atom[ct])
                    if not Info:
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
                    else:
                        generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                        if 'Nat. Abund.' not in generalData['Isotopes'][atom[ct]]:
                            isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                            generalData['Isotope'][atom[ct]] = isotope
                    generalData['Color'].append(Info['Color'])
                    # if generalData['Type'] == 'magnetic':
                    #     if len(landeg) < len(generalData['AtomTypes']):
                    #         landeg.append(2.0)
            atmData['Drawing']['Atoms'] = []
            for atom in atmData['Atoms']:
                atmData['Drawing']['Atoms'].append(G2mth.MakeDrawAtom(atmData,atom))

        def onCancel(event,page=0):
            'complete or bail out from RB define, cleaning up'
            G2frame.rbBook.DeletePage(G2frame.rbBook.FindPage(pagename))
            G2frame.rbBook.SetSelection(page)

        def Page1():
            '''Show the GUI for first stage of the rigid body with all atoms in 
            phase in crystal coordinates. Select the atoms to go onto the
            next stage
            '''
            def ShowSelection(selections):
                'respond to change in atom selections'
                ct,cs = [1,8]
                generalData = rd.Phase['General']
                for i,atom in enumerate(rd.Phase['Drawing']['Atoms']):
                    if i in selections:
                        factor = 1
                    else:
                        factor = 2.5
                    atNum = generalData['AtomTypes'].index(atom[ct]) 
                    atom[cs] = list(np.array(generalData['Color'][atNum])//factor) 
                draw(*drawArgs)
            def onPage1OK(event):
                '1st section has been completed, move onto next'
                G2frame.G2plotNB.Delete(rd.Phase['General']['Name'])
                GetCoords(atmsel)
                Page2()

            if 'macromolecular' == rd.Phase['General']['Type']:
                # for PDB imports, lets see if a quick reformat of atoms list will work
                rd.Phase['Atoms'] = [a[3:] for a in rd.Phase['Atoms']]
                rd.Phase['General']['AtomPtrs']  = [i-3 for i in rd.Phase['General']['AtomPtrs']]
                rd.Phase['General']['Type'] = 'nuclear'
            SetupDrawing(rd.Phase) # add information to reader object to allow plotting
            atomlist = [atom[0] for atom in rd.Phase['Atoms']]
            atmsel = list(range(len(rd.Phase['Atoms'])))
            # broken -- # why no bonds?
            #for atm in rd.Phase['Drawing']['Atoms']:
            #    atm[6] = 'balls & sticks' 

            draw,drawArgs = G2plt.PlotStructure(G2frame,rd.Phase,True)
            ShowSelection(atmsel)

            if G2frame.rbBook.FindPage(pagename) is not None:
                G2frame.rbBook.DeletePage(G2frame.rbBook.FindPage(pagename))

            RBImp = wx.ScrolledWindow(G2frame.rbBook)
            RBImpPnl = wx.Panel(RBImp)
            G2frame.rbBook.AddPage(RBImp,pagename)
            G2frame.rbBook.SetSelection(G2frame.rbBook.FindPage(pagename))

            HelpInfo = '''
This window shows all the atoms that were read from the 
selected phase file. Select the atoms that will be used in the
rigid body processing (this may include atoms needed to 
define an axis or origin that will not be included in the 
eventual rigid body.) Note that in the plot window, 
unselected atoms appear much darker than selected atoms.
'''
            mainSizer = G2G.G2MultiChoiceWindow(RBImpPnl,
                            'Select atoms to import',
                            atomlist,atmsel,OnChange=ShowSelection,
                            helpText=HelpInfo)

            # OK/Cancel buttons        
            btnsizer = wx.StdDialogButtonSizer()
            OKbtn = wx.Button(RBImpPnl, wx.ID_OK, 'Continue')
            OKbtn.SetDefault()
            btnsizer.AddButton(OKbtn)
            OKbtn.Bind(wx.EVT_BUTTON,onPage1OK)
            btn = wx.Button(RBImpPnl, wx.ID_CANCEL)
            btn.Bind(wx.EVT_BUTTON,onCancel)
            btnsizer.AddButton(btn)
            btnsizer.Realize()
            mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER,50)

            RBImpPnl.SetSizer(mainSizer,True)

            mainSizer.Layout()    
            Size = mainSizer.GetMinSize()
            Size[0] += 40
            Size[1] = max(Size[1],G2frame.GetSize()[1]-200) + 20
            RBImpPnl.SetSize(Size)
            RBImp.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
            RBImp.Scroll(0,0)

        def Page2():
            '''Show the GUI for the second stage, where selected atoms are
            now in Cartesian space, manipulate the axes and export selected 
            atoms to a vector or residue rigid body.
            '''
            def UpdateDraw(event=None):
                'Called when info changes in grid, replots'
                UpdateVectorBody(rbData)
                DrawCallback()
                
            def onSetAll(event):
                'Set all atoms as selected'
                grid.completeEdits()
                for i in range(len(rd.Phase['RBselection'])):
                    rd.Phase['RBselection'][i] = 1 # table needs 0/1 for T/F
                grid.ForceRefresh()
                UpdateDraw()
                
            def onToggle(event):
                'Toggles selection state for all atoms'
                grid.completeEdits()
                for i in range(len(rd.Phase['RBselection'])):
                    rd.Phase['RBselection'][i] = int(not rd.Phase['RBselection'][i])
                grid.ForceRefresh()
                UpdateDraw()
                
            def onSetOrigin(event):
                'Resets origin to midpoint between all selected atoms'
                grid.completeEdits()
                center = np.array([0.,0.,0.])
                count = 0
                for i in range(len(rd.Phase['RBselection'])):
                    if rd.Phase['RBselection'][i]:
                        count += 1
                        center += rd.Phase['RBcoords'][i]
                if count:
                    rd.Phase['RBcoords'] -= center/count
                grid.ForceRefresh()
                UpdateDraw()
                
            def onSetX(event):
                grid.completeEdits()
                center = np.array([0.,0.,0.])
                count = 0
                for i in range(len(rd.Phase['RBselection'])):
                    if rd.Phase['RBselection'][i]:
                        count += 1
                        center += rd.Phase['RBcoords'][i]
                if not count:
                    G2G.G2MessageBox(G2frame,'No atoms selected',
                                    'Selection required')
                    return
                XYZP = center/count
                if np.sqrt(sum(XYZP**2)) < 0.1:
                    G2G.G2MessageBox(G2frame,
                            'The selected atom(s) are too close to the origin',
                            'near origin')
                    return
                if bntOpts['direction'] == 'y':
                    YP = XYZP / np.sqrt(np.sum(XYZP**2))
                    ZP = np.cross((1,0,0),YP)
                    if sum(ZP*ZP) < .1: # pathological condition: Y' along X
                        ZP = np.cross((0,0,1),YP)
                    XP = np.cross(YP,ZP)
                elif bntOpts['direction'] == 'z':
                    ZP = XYZP / np.sqrt(np.sum(XYZP**2))
                    XP = np.cross((0,1,0),ZP)
                    if sum(XP*XP) < .1: # pathological condition: X' along Y
                        XP = np.cross((0,0,1),ZP)
                    YP = np.cross(ZP,XP)
                else:
                    XP = XYZP / np.sqrt(np.sum(XYZP**2))
                    YP = np.cross((0,0,1),XP)
                    if sum(YP*YP) < .1: # pathological condition: X' along Z
                        YP = np.cross((0,1,0),XP)
                    ZP = np.cross(XP,YP)
                trans = np.array((XP,YP,ZP))
                # update atoms in place
                rd.Phase['RBcoords'][:] = np.inner(trans,rd.Phase['RBcoords']).T
                grid.ForceRefresh()
                UpdateDraw()

            def onSetPlane(event): 
                '''Finds eigen vector/matrix for best "ellipsoid" about atoms; 
                rotate atoms so that smallest axis is along choice. 
                '''
                grid.completeEdits()
                selList = [i==1 for i in rd.Phase['RBselection']]
                XYZ = rd.Phase['RBcoords'][selList]
                Natoms = len(XYZ)
                if Natoms < 3: 
                    G2G.G2MessageBox(G2frame,'A plane requires three or more atoms','Need more atoms')
                    return
                Zmat = np.zeros((3,3))
                for xyz in XYZ:
                    Zmat += np.outer(xyz.T,xyz)
                Evec,Emat = nl.eig(Zmat)
                Order = np.argsort(np.nan_to_num(Evec))     #short-long order
                if bntOpts['plane'] == 'xy':        #short along z
                    trans = np.array([Emat[Order[2]],Emat[Order[1]],Emat[Order[0]]])
                elif bntOpts['plane'] == 'yz':      #short along x
                    trans = np.array([Emat[Order[0]],Emat[Order[2]],Emat[Order[1]]])
                elif bntOpts['plane'] == 'xz':      #short along y
                    trans = np.array([Emat[Order[1]],Emat[Order[0]],Emat[Order[2]]])
                else:
                    print('unexpected plane',bntOpts['plane'])
                    return
                # update atoms in place
                rd.Phase['RBcoords'][:] = np.inner(trans,rd.Phase['RBcoords']).T
                grid.ForceRefresh()
                UpdateDraw()

            def onWriteXYZ(event):
                '''Writes selected atoms in a .xyz file for use in Avogadro, etc.
                '''
                grid.completeEdits()
                center = np.array([0.,0.,0.])
                count = 0
                for i in range(len(rd.Phase['RBselection'])):
                    if rd.Phase['RBselection'][i]:
                        count += 1
                        center += rd.Phase['RBcoords'][i]
                if count:
                    center /= count
                else:
                    print('nothing selected')
                    return
                obj = G2IO.ExportBaseclass(G2frame,'XYZ','.xyz')
                #obj.InitExport(None)
                if obj.ExportSelect():    # set export parameters; ask for file name
                    return
                obj.OpenFile()
                obj.Write(str(count))
                obj.Write('')
                for i in range(len(rd.Phase['RBselection'])):
                    if rd.Phase['RBselection'][i]:
                        line = ' ' + rd.Phase['RBtypes'][i]
                        for xyz in rd.Phase['RBcoords'][i]:
                            line += ' ' + str(xyz)
                        obj.Write(line)
                obj.CloseFile()
                #GSASIIpath.IPyBreak()
                
            def onAddVector(event):
                '''Adds selected atoms as a new vector rigid body.
                Closes out the importer tab when done. 
                '''
                grid.completeEdits()
                name = os.path.splitext(os.path.split(filename)[1])[0]
                namelist = [data['Vector'][key]['RBname'] for key in
                        data['Vector'] if 'RBname' in data['Vector'][key]]
                name = G2obj.MakeUniqueLabel(name,namelist)
                rb = MakeVectorBody(name)
                UpdateVectorBody(rb,True)
                if len(rb['rbTypes']) < 3: return # must have at least 3 atoms
                rbid = ran.randint(0,sys.maxsize)
                data['Vector'][rbid] = rb
                data['RBIds']['Vector'].append(rbid)
                for t in rb['rbTypes']:
                    if t in data['Vector']['AtInfo']: continue
                    Info = G2elem.GetAtomInfo(t)
                    data['Vector']['AtInfo'][t] = [Info['Drad'],Info['Color']]
                G2frame.G2plotNB.Delete('Rigid body')
                onCancel(event,0)
                
            def onAddResidue(event):
                '''Adds selected atoms as a new residue rigid body.
                Closes out the importer tab when done. 
                '''
                grid.completeEdits()
                name = os.path.split(filename)[1]
                rbXYZ = []
                rbTypes = []
                atNames = []
                for i in rd.Phase['RBindex']:
                    if rd.Phase['RBselection'][i]:
                        rbXYZ.append(rd.Phase['RBcoords'][i])
                        rbTypes.append(rd.Phase['RBtypes'][i])
                        atNames.append(rd.Phase['RBlbls'][i])
                if len(rbTypes) < 3: return # must have at least 3 atoms
                rbXYZ = np.array(rbXYZ)
                rbid = ran.randint(0,sys.maxsize)
                namelist = [data['Residue'][key]['RBname'] for key in
                        data['Residue'] if 'RBname' in data['Residue'][key]]
                name = G2obj.MakeUniqueLabel(name,namelist)
                data['Residue'][rbid] = {'RBname':name,'rbXYZ':rbXYZ,
                    'rbTypes':rbTypes,'atNames':atNames,'rbRef':[0,1,2,False],
                    'rbSeq':[],'SelSeq':[0,0],'useCount':0}
                data['RBIds']['Residue'].append(rbid)
                for t in rbTypes:
                    if t in data['Residue']['AtInfo']: continue
                    Info = G2elem.GetAtomInfo(t)
                    data['Residue']['AtInfo'][t] = [Info['Drad'],Info['Color']]

                print ('Rigid body added')
                G2frame.G2plotNB.Delete('Rigid body')
                onCancel(event,1)

            if G2frame.rbBook.FindPage(pagename) is not None:
                G2frame.rbBook.DeletePage(G2frame.rbBook.FindPage(pagename))
            RBImp = wx.ScrolledWindow(G2frame.rbBook)
            RBImpPnl = wx.Panel(RBImp)
            G2frame.rbBook.AddPage(RBImp,pagename)
            G2frame.rbBook.SetSelection(G2frame.rbBook.FindPage(pagename))
            AtInfo = {}
            for t in rd.Phase['RBtypes']:
                if t in AtInfo: continue
                Info = G2elem.GetAtomInfo(t)
                AtInfo[t] = [Info['Drad'],Info['Color']]
            plotDefaults = {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':30.,'viewDir':[0,0,1],}

            rd.Phase['RBindex'] = list(range(len(rd.Phase['RBtypes'])))
            rd.Phase['RBselection'] = len(rd.Phase['RBtypes']) * [1]
            name = 'UNKRB'
            namelist = [data['Vector'][key]['RBname'] for key in
                        data['Vector'] if 'RBname' in data['Vector'][key]]
            name = G2obj.MakeUniqueLabel(name,namelist)
            rbData = MakeVectorBody()
            DrawCallback = G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,rbData,plotDefaults)

            mainSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer = wx.BoxSizer(wx.VERTICAL)
            helpText = '''
In this window, if wanted,
one can select one or more atoms and use them
to define an origin, a specified axis or place the selected atoms into 
a selected plane. (Different sets of atoms can be used for each
operation.)
%%Once that is done, atoms can be selected and can be exported in a
"XYZ" file for use in a program such as Avogadro or can be used to 
create a Vector or Residue rigid body. 
'''
            btnSizer.Add(G2G.HelpButton(RBImpPnl,helpText,wrap=400),
                             0,wx.ALIGN_RIGHT)
            btnSizer.Add(wx.StaticText(RBImpPnl,wx.ID_ANY,'Reorder atoms by dragging'),0,wx.ALL)
            btnSizer.Add((-1,15))
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'Set All')
            btn.Bind(wx.EVT_BUTTON,onSetAll)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'Toggle')
            btn.Bind(wx.EVT_BUTTON,onToggle)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)
            btnSizer.Add((-1,15))
            btnSizer.Add(wx.StaticText(RBImpPnl,wx.ID_ANY,'Reorient using selected\natoms...'),0,wx.ALL)
            btnSizer.Add((-1,5))
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'Set origin')
            btn.Bind(wx.EVT_BUTTON,onSetOrigin)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)

            bntOpts = {'plane':'xy','direction':'x'}
            inSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'Place in plane')
            btn.Bind(wx.EVT_BUTTON,onSetPlane)
            inSizer.Add(btn)
            inSizer.Add(G2G.G2ChoiceButton(RBImpPnl,('xy','yz','xz'),None,None,bntOpts,'plane'))
            btnSizer.Add(inSizer,0,wx.ALIGN_CENTER)
            
            inSizer = wx.BoxSizer(wx.HORIZONTAL)
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'Define as')
            btn.Bind(wx.EVT_BUTTON,onSetX)
            inSizer.Add(btn)
            inSizer.Add(G2G.G2ChoiceButton(RBImpPnl,('x','y','z'),None,None,bntOpts,'direction'))
            btnSizer.Add(inSizer,0,wx.ALIGN_CENTER)
            
            btnSizer.Add((-1,15))
            btnSizer.Add(wx.StaticText(RBImpPnl,wx.ID_ANY,'Use selected atoms to\ncreate...'),0,wx.ALL)
            btnSizer.Add((-1,5))
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'export as xyz')
            btn.Bind(wx.EVT_BUTTON,onWriteXYZ)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)
            btnSizer.Add((-1,10))
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'a Vector Body')
            btn.Bind(wx.EVT_BUTTON,onAddVector)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)
            btn = wx.Button(RBImpPnl, wx.ID_ANY, 'a Residue Body')
            btn.Bind(wx.EVT_BUTTON,onAddResidue)
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)
            btn = wx.Button(RBImpPnl, wx.ID_CANCEL)
            btn.Bind(wx.EVT_BUTTON,onCancel)
            btnSizer.Add((-1,10))
            btnSizer.Add(btn,0,wx.ALIGN_CENTER)

            mainSizer.Add(btnSizer)
            mainSizer.Add((5,5))
            grid = DragableRBGrid(RBImpPnl,rd.Phase,UpdateDraw)
            mainSizer.Add(grid)
            RBImpPnl.SetSizer(mainSizer,True)
            mainSizer.Layout()    
            Size = mainSizer.GetMinSize()
            Size[0] += 40
            Size[1] = max(Size[1],G2frame.GetSize()[1]-200) + 20
            RBImpPnl.SetSize(Size)
            RBImp.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
            RBImp.Scroll(0,0)

        def GetCoords(atmsel):
            '''Create orthogonal coordinates for selected atoms.
            Place the origin at the center of the body
            '''
            atms = rd.Phase['Atoms']
            cell = rd.Phase['General']['Cell'][1:7]
            Amat,Bmat = G2lat.cell2AB(cell)
            rd.Phase['RBcoords'] = np.array([np.inner(Amat,atms[i][3:6]) for i in atmsel])
            rd.Phase['RBcoords'] -= rd.Phase['RBcoords'].mean(axis=0)  # origin to middle
            rd.Phase['RBtypes'] = [atms[i][1] for i in atmsel]
            rd.Phase['RBlbls'] = [atms[i][0] for i in atmsel]

        def UpdateVectorBody(rb,useSelection=False):
            '''Put the atoms in order to pass for plotting or for storage as 
            a vector rigid body. 

            :param dict rb: rigid body contents created in :func:`MakeVectorBody`
            :param bool useSelection: True if the rd.Phase['RBselection']
              values will be used to select which atoms are included in the 
              rigid body. If False (default) they are included in rb
              and are used for plotting.          
            '''
            coordlist = []
            typeslist = []
            sellist = []
            for i in rd.Phase['RBindex']:
                use = True
                if useSelection and not rd.Phase['RBselection'][i]: use = False
                if use:
                    coordlist.append(rd.Phase['RBcoords'][i])
                    typeslist.append(rd.Phase['RBtypes'][i])
                    sellist.append(rd.Phase['RBselection'][i])
            coordlist = np.array(coordlist)
            rb['rbXYZ'] = coordlist
            rb['rbVect'] = [coordlist]
            rb['rbTypes'] = typeslist
            if not useSelection:
                rb['Selection'] = sellist
            elif 'Selection' in rb:
                del rb['Selection']

        def MakeVectorBody(name=''):
            '''Make the basic vector rigid body dict (w/o coordinates) used for
            export and for plotting
            '''
            vecMag = [1.0]
            vecRef = [False]
            rb = {'RBname':name,'VectMag':vecMag,
                    'rbRef':[0,1,2,False],'VectRef':vecRef,
                    'useCount':0}
            UpdateVectorBody(rb)
            return rb

        # too lazy to figure out why wx crashes
        if wx.__version__.split('.')[0] != '4':
            wx.MessageBox('Sorry, wxPython 4.x is required to run this command',
                                  caption='Update Python',
                                  style=wx.ICON_EXCLAMATION)
            return
        if platform.python_version()[:1] == '2':
            wx.MessageBox('Sorry, Python >=3.x is required to run this command',
                                  caption='Update Python',
                                  style=wx.ICON_EXCLAMATION)
            return

        # get importer type and a phase file of that type
        G2sc.LoadG2fil()
        choices = [rd.formatName for  rd in G2sc.Readers['Phase']] 
        dlg = G2G.G2SingleChoiceDialog(G2frame,'Select the format of the file',
                                     'select format',choices)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                col = dlg.GetSelection()
            else:
                col = None
                return
        finally:
            dlg.Destroy()
        reader = G2sc.Readers['Phase'][col]

        choices = reader.formatName + " file ("
        w = ""
        for extn in reader.extensionlist:
            if w != "": w += ";"
            w += "*" + extn
        choices += w + ")|" + w
        #choices += "|zip archive (.zip)|*.zip"
        if not reader.strictExtension:
            choices += "|any file (*.*)|*.*"
        typ = '( type '+reader.formatName+')'
        filelist = G2G.GetImportFile(G2frame,
                        message="Choose phase input file"+typ,
                        defaultFile="",wildcard=choices,style=wx.FD_OPEN)
        if len(filelist) != 1: return

        # read in the phase file
        filename = filelist[0]
        rd = reader
        with open(filename, 'r'):
            rd.ReInitialize()
            rd.errors = ""
            if not rd.ContentsValidator(filename):   # Report error
                G2fl.G2Print("Warning: File {} has a validation error".format(filename))
                return
            if len(rd.selections) > 1:
                print("File {} has {} phases. This is unexpected."
                                    .format(filename,len(rd.selections)))
                return

            rd.objname = os.path.basename(filename)
            try:
                rd.Reader(filename)
            except Exception as msg:
                G2fl.G2Print("Warning: read of file {} failed\n{}".format(
                    filename,rd.errors))
                if GSASIIpath.GetConfigValue('debug'):
                    print(msg)
                    import traceback
                    print (traceback.format_exc())
                    GSASIIpath.IPyBreak()
                return

        pagename = 'Rigid body importer'
        Page1()
        return

    def AddVectTrans(event):
        'Add a translation to an existing vector rigid body'
        choices = []
        rbIdlist = []
        for rbid in data['RBIds']['Vector']:
            if rbid != 'AtInfo':
                rbIdlist.append(rbid)
                choices.append(data['Vector'][rbid]['RBname'])
        if len(choices) == 0:
            G2G.G2MessageBox(G2frame,'No Vector Rigid Bodies found',
                                 'No VR Bodies')
            return
        elif len(choices) == 1:
            rbid = rbIdlist[0]
        else:
            dlg = G2G.G2SingleChoiceDialog(G2frame,'Select the rigid body to save',
                                  'select format',choices)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    rbid = rbIdlist[dlg.GetSelection()]
                else:
                    return
            finally:
                dlg.Destroy()
        data['Vector'][rbid]['VectMag'] += [1.0]
        data['Vector'][rbid]['VectRef'] += [False]
        nAtoms = len(data['Vector'][rbid]['rbXYZ'])
        data['Vector'][rbid]['rbVect'] += [np.zeros((nAtoms,3))]
        UpdateVectorRB()
        
    def SaveVectorRB(event):
        choices = []
        rbIdlist = []
        for rbid in data['RBIds']['Vector']:
            if rbid != 'AtInfo':
                rbIdlist.append(rbid)
                choices.append(data['Vector'][rbid]['RBname'])
        if len(choices) == 0:
            G2G.G2MessageBox(G2frame,'No Vector Rigid Bodies found',
                                 'No VR Bodies')
            return
        elif len(choices) == 1:
            rbid = rbIdlist[0]
        else:
            dlg = G2G.G2SingleChoiceDialog(G2frame,'Select the rigid body to save',
                                  'select format',choices)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    rbid = rbIdlist[dlg.GetSelection()]
                else:
                    return
            finally:
                dlg.Destroy()
             
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose file to save vector rigid body',
            pth, '', 'VRB files (*.vecbody)|*.vecbody',
            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.vecbody'  # set extension
                fp = open(filename,'w')
                fp.write('Name: '+data['Vector'][rbid]['RBname']+'\n')
                fp.write('Trans: ')
                for i in data['Vector'][rbid]['VectMag']:
                    fp.write(str(i)+" ") 
                fp.write('\n')
                ntrans = len(data['Vector'][rbid]['VectMag'])
                for i,sym in enumerate(data['Vector'][rbid]['rbTypes']):
                    fp.write("{:3s}".format(sym))
                    for j in range(ntrans):
                        fp.write('{:8.5f}{:9.5f}{:9.5f}   '
                            .format(*data['Vector'][rbid]['rbVect'][j][i]))
                    fp.write('\n')
                fp.close()
                print ('Vector rigid body saved to: '+filename)
        finally:
            dlg.Destroy()
            
    def ReadVectorRB(event):
        AtInfo = data['Vector']['AtInfo']
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose file to read vector rigid body',
            pth, '', 'VRB files (*.vecbody)|*.vecbody',
            wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.vecbody'  # set extension
                fp = open(filename,'r')
                l = fp.readline().strip()
                if 'Name' not in l:
                    fp.close()
                    G2frame.ErrorDialog('Read Error',
                        'File '+filename+' does not start with Name\nFirst line ='
                        +l+'\ninvalid file',parent=G2frame)
                    return
                name = l.split(':')[1].strip()
                trans = fp.readline().strip().split(':')[1].split()
                vecMag = [float(i) for i in trans]
                ntrans = len(trans)
                vecs = [[] for i in range(ntrans)]
                types = []
                l = fp.readline().strip()
                while l:
                    nums = l.strip().split()
                    types.append(nums.pop(0))
                    t = types[-1]
                    if t not in AtInfo:
                        Info = G2elem.GetAtomInfo(t)
                        AtInfo[t] = [Info['Drad'],Info['Color']]
                    for i in range(ntrans):
                        vecs[i].append([float(nums.pop(0)) for j in range(3)])
                    l = fp.readline().strip()
                fp.close()
            else:
                return        
        finally:
            dlg.Destroy()
        natoms = len(types)
        vecs = [np.array(vecs[i]) for i in range(ntrans)]
        rbid = ran.randint(0,sys.maxsize)
        namelist = [data['Vector'][key]['RBname'] for key in data['Vector']
                        if 'RBname' in data['Vector'][key]]
        name = G2obj.MakeUniqueLabel(name,namelist)
        data['Vector'][rbid] = {'RBname':name,'VectMag':vecMag,
                'rbXYZ':np.zeros((natoms,3)),
                'rbRef':[0,1,2,False],'VectRef':ntrans*[False],
                'rbTypes':types,
                'rbVect':vecs,'useCount':0}
        data['RBIds']['Vector'].append(rbid)
        UpdateVectorRB()
        
    def AddResidueRB(event):
        global resRBsel
        AtInfo = data['Residue']['AtInfo']
        macro = getMacroFile('rigid body')
        if not macro:
            return
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'I' == items[0]:
                resRBsel = ran.randint(0,sys.maxsize)
                rbName = items[1]
                rbTypes = []
                rbXYZ = []
                rbSeq = []
                atNames = []
                nAtms,nSeq,nOrig,mRef,nRef = [int(items[i]) for i in [2,3,4,5,6]]
                for iAtm in range(nAtms):
                    macStr = macro.readline().split()
                    atName = macStr[0]
                    atType = macStr[1]
                    atNames.append(atName)
                    rbXYZ.append([float(macStr[i]) for i in [2,3,4]])
                    rbTypes.append(atType)
                    if atType not in AtInfo:
                        Info = G2elem.GetAtomInfo(atType)
                        AtInfo[atType] = [Info['Drad'],Info['Color']]
                rbXYZ = np.array(rbXYZ)-np.array(rbXYZ[nOrig-1])
                for iSeq in range(nSeq):
                    macStr = macro.readline().split()
                    mSeq = int(macStr[0])
                    for jSeq in range(mSeq):
                        macStr = macro.readline().split()
                        iBeg = int(macStr[0])-1
                        iFin = int(macStr[1])-1
                        angle = 0.0
                        nMove = int(macStr[2])
                        iMove = [int(macStr[i])-1 for i in range(3,nMove+3)]
                        rbSeq.append([iBeg,iFin,angle,iMove])
                namelist = [data['Residue'][key]['RBname'] for key in
                           data['Residue'] if 'RBname' in data['Residue'][key]]
                rbName = G2obj.MakeUniqueLabel(rbName,namelist)
                data['Residue'][resRBsel] = {'RBname':rbName,'rbXYZ':rbXYZ,'rbTypes':rbTypes,
                    'atNames':atNames,'rbRef':[nOrig-1,mRef-1,nRef-1,True],'rbSeq':rbSeq,
                    'SelSeq':[0,0],'useCount':0,'molCent':None}
                data['RBIds']['Residue'].append(resRBsel)
                print ('Rigid body '+rbName+' added')
            macStr = macro.readline()
        macro.close()
        UpdateResidueRB()
        
    def ImportResidueRB():
        global resRBsel
        AtInfo = data['Residue']['AtInfo']
        text,ext = getTextFile()
        if not text:
            return
        resRBsel = ran.randint(0,sys.maxsize)
        rbTypes = []
        rbXYZ = []
        atNames = []
        txtStr = text.readline()
        if 'xyz' in ext:
            txtStr = text.readline()
            txtStr = text.readline()
        elif 'mol2' in ext:
            while 'ATOM' not in txtStr:
                txtStr = text.readline()
            txtStr = text.readline()
        elif 'pdb' in ext:
            while 'ATOM' not in txtStr[:6] and 'HETATM' not in txtStr[:6]:
                txtStr = text.readline()
        items = txtStr.split()
        nat = 1
        while len(items):
            if 'txt' in ext:
                atName = items[0]
                atType = items[1]
                rbXYZ.append([float(items[i]) for i in [2,3,4]])
            elif 'xyz' in ext:
                atType = items[0]
                rbXYZ.append([float(items[i]) for i in [1,2,3]])
                atName = '%s%d'%(atType,nat)
            elif 'mol2' in ext:
                atType = items[1]
                atName = items[1]+items[0]
                rbXYZ.append([float(items[i]) for i in [2,3,4]])
            elif 'pdb' in ext:
                atType = items[-1]
                if not items[2][-1].isnumeric():
                    atName = '%s%d'%(items[2],nat)
                else:
                    atName = '5s'%items[2]
                xyz = txtStr[30:55].split()                    
                rbXYZ.append([float(x) for x in xyz])
            atNames.append(atName)
            rbTypes.append(atType)
            if atType not in AtInfo:
                Info = G2elem.GetAtomInfo(atType)
                AtInfo[atType] = [Info['Drad'],Info['Color']]
            txtStr = text.readline()
            if 'mol2' in ext and 'BOND' in txtStr:
                break
            if 'pdb' in ext and ('ATOM' not in txtStr[:6] and 'HETATM' not in txtStr[:6]):
                break
            items = txtStr.split()
            nat += 1
        if len(atNames) < 3:
            G2G.G2MessageBox(G2frame,'Not enough atoms in rigid body; must be 3 or more')
        else:
            rbXYZ = np.array(rbXYZ)-np.array(rbXYZ[0])
            Xxyz = rbXYZ[1]
            X = Xxyz/np.sqrt(np.sum(Xxyz**2))
            Yxyz = rbXYZ[2]
            Y = Yxyz/np.sqrt(np.sum(Yxyz**2))
            Mat = G2mth.getRBTransMat(X,Y)
            rbXYZ = np.inner(Mat,rbXYZ).T
            name = 'UNKRB'
            namelist = [data['Residue'][key]['RBname'] for key in data['Residue']
                        if 'RBname' in data['Residue'][key]]
            name = G2obj.MakeUniqueLabel(name,namelist)
            data['Residue'][resRBsel] = {'RBname':name,'rbXYZ':rbXYZ,'rbTypes':rbTypes,
                'atNames':atNames,'rbRef':[0,1,2,False],'rbSeq':[],'SelSeq':[0,0],'useCount':0,'molCent':False}
            data['RBIds']['Residue'].append(resRBsel)
            print ('Rigid body UNKRB added')
        text.close()
        UpdateResidueRB()
        
    def SaveResidueRB():
        global resRBsel
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(G2frame, 'Choose PDB file for Atom XYZ', pth, '', 
            'PDB files (*.pdb)|*.pdb',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                filename = os.path.splitext(filename)[0]+'.pdb'  # make extension .pdb
                File = open(filename,'w')       
                rbData =  data['Residue'][resRBsel]
                for iat,xyz in enumerate(rbData['rbXYZ']):
                    File.write('ATOM %6d  %-4s%3s     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'%(
                        iat,rbData['atNames'][iat],rbData['RBname'][:3],xyz[0],xyz[1],xyz[2],rbData['rbTypes'][iat]))
                File.close()
                print ('Atom XYZ saved to: '+filename)
        finally:
            dlg.Destroy()
            
    def AddSpinRB(event):
        
        rbid = ran.randint(0,sys.maxsize)
        atType = 'C'
        rbType = 'Q'
        Natoms = 1
        name = 'UNKRB'
        namelist = [data['Spin'][key]['RBname'] for key in data['Spin']]
        name = G2obj.MakeUniqueLabel(name,namelist)
        atColor = G2elem.GetAtomInfo(atType)['Color']
        data['Spin'][rbid] = {'RBname':name,'Natoms':Natoms,'atType':atType,'rbType':rbType,'atColor':atColor,
            'useCount':0,'nSH':0,'SHC':[{},],'Radius':[1.0,False],'Matrix':np.eye(3),'rbPos':np.zeros(3)}
        data['RBIds']['Spin'].append(rbid)
        UpdateSpinRB()
        
    def FindNeighbors(Orig,XYZ,atTypes,atNames,AtInfo):
        Radii = []
        for Atype in atTypes:
            Radii.append(AtInfo[Atype][0])
        Radii = np.array(Radii)
        Neigh = []
        Dx = XYZ-XYZ[Orig]
        dist = np.sqrt(np.sum(Dx**2,axis=1))
        sumR = Radii[Orig]+Radii
        IndB = ma.nonzero(ma.masked_greater(dist-0.85*sumR,0.))
        for j in IndB[0]:
            if j != Orig and atTypes[j] != 'H':
                Neigh.append(atNames[j])
        return Neigh
        
    def FindAllNeighbors(XYZ,atTypes,atNames,AtInfo):
        NeighDict = {}
        for iat,xyz in enumerate(atNames):
            NeighDict[atNames[iat]] = FindNeighbors(iat,XYZ,atTypes,atNames,AtInfo)
        return NeighDict
        
    def FindRiding(Orig,Pivot,NeighDict):
        riding = [Orig,Pivot]
        iAdd = 1
        new = True
        while new:
            newAtms = NeighDict[riding[iAdd]]
            for At in newAtms:
                new = False
                if At not in riding:
                    riding.append(At)
                    new = True
            iAdd += 1
            if iAdd < len(riding):
                new = True
        return riding[2:]
                        
    def OnDefineTorsSeq(event):
        global resRBsel
        rbData = data['Residue'][resRBsel]
        if not len(rbData):
            return
        atNames = rbData['atNames']
        AtInfo = data['Residue']['AtInfo']
        atTypes = rbData['rbTypes']
        XYZ = rbData['rbXYZ']
        neighDict = FindAllNeighbors(XYZ,atTypes,atNames,AtInfo)
        TargList = []            
        dlg = wx.SingleChoiceDialog(G2frame,'Select origin atom for torsion sequence','Origin atom',rbData['atNames'])
        if dlg.ShowModal() == wx.ID_OK:
            Orig = dlg.GetSelection()
            TargList = neighDict[atNames[Orig]]
        dlg.Destroy()
        if not len(TargList):
            return
        dlg = wx.SingleChoiceDialog(G2frame,'Select pivot atom for torsion sequence','Pivot atom',TargList)
        if dlg.ShowModal() == wx.ID_OK:
            Piv = atNames.index(TargList[dlg.GetSelection()])
            riding = FindRiding(atNames[Orig],atNames[Piv],neighDict)
            Riding = []
            for atm in riding:
                Riding.append(atNames.index(atm))
            rbData['rbSeq'].append([Orig,Piv,0.0,Riding])            
        dlg.Destroy()
        UpdateResidueRB()
        

    def UpdateVectorRB(Scroll=0):
        '''Display & edit a selected Vector RB
        '''
        global resRBsel
        def rbNameSizer(rbid,rbData):

            def OnRBName(event):
                event.Skip()
                Obj = event.GetEventObject()
                name = Obj.GetValue()
                if name == rbData['RBname']: return # no change
                namelist = [data['Vector'][key]['RBname'] for key in
                        data['Vector'] if 'RBname' in data['Vector'][key]]
                name = G2obj.MakeUniqueLabel(name,namelist)
                rbData['RBname'] = name
                wx.CallAfter(UpdateVectorRB)
                
            def OnDelRB(event):
                Obj = event.GetEventObject()
                rbid = Indx[Obj.GetId()]
                if rbid in data['Vector']:
                    del data['Vector'][rbid]
                    data['RBIds']['Vector'].remove(rbid)
                    rbData['useCount'] -= 1
                wx.CallAfter(UpdateVectorRB)
                
            def OnPlotRB(event):
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,rbData,plotDefaults)
            
            # start of rbNameSizer
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(VectorRBDisplay,-1,'Rigid body name: '),0,WACV)
            RBname = wx.TextCtrl(VectorRBDisplay,-1,rbData['RBname'],
                                     style=wx.TE_PROCESS_ENTER)
            RBname.Bind(wx.EVT_LEAVE_WINDOW, OnRBName)
            RBname.Bind(wx.EVT_TEXT_ENTER,OnRBName)
            RBname.Bind(wx.EVT_KILL_FOCUS,OnRBName)
            nameSizer.Add(RBname,0,WACV)
            nameSizer.Add((5,0),)
            plotRB =  wx.Button(VectorRBDisplay,wx.ID_ANY,'Plot',style=wx.BU_EXACTFIT)
            plotRB.Bind(wx.EVT_BUTTON, OnPlotRB)
            Indx[plotRB.GetId()] = rbid
            nameSizer.Add(plotRB,0,WACV)
            nameSizer.Add((5,0),)
            if not rbData['useCount']:
                delRB = wx.Button(VectorRBDisplay,wx.ID_ANY,"Delete",style=wx.BU_EXACTFIT)
                delRB.Bind(wx.EVT_BUTTON, OnDelRB)
                Indx[delRB.GetId()] = rbid
                nameSizer.Add(delRB,0,WACV)
            nameSizer.Add((-1,-1),1,wx.EXPAND,1)
            nameSizer.Add(G2G.HelpButton(VectorRBDisplay,helpIndex=G2frame.dataWindow.helpKey))
            return nameSizer
            
        def rbRefAtmSizer(rbid,rbData):
            
            def OnRefSel(event):
                Obj = event.GetEventObject()
                iref = Indx[Obj.GetId()]
                sel = Obj.GetValue()
                rbData['rbRef'][iref] = atNames.index(sel)
                FillRefChoice(rbid,rbData)
            
            refAtmSizer = wx.BoxSizer(wx.HORIZONTAL)
            atNames = [name+str(i) for i,name in enumerate(rbData['rbTypes'])]
            rbRef = rbData.get('rbRef',[0,1,2,False])
            rbData['rbRef'] = rbRef
            if rbData['useCount']:
                refAtmSizer.Add(wx.StaticText(VectorRBDisplay,-1,
                    'Orientation reference atoms A-B-C: %s, %s, %s'%(atNames[rbRef[0]], \
                     atNames[rbRef[1]],atNames[rbRef[2]])),0)
            else:
                refAtmSizer.Add(wx.StaticText(VectorRBDisplay,-1,
                    'Orientation reference atoms A-B-C: '),0,WACV)
                for i in range(3):
                    choices = [atNames[j] for j in refChoice[rbid][i]]
                    refSel = wx.ComboBox(VectorRBDisplay,-1,value='',
                        choices=choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    refSel.SetValue(atNames[rbRef[i]])
                    refSel.Bind(wx.EVT_COMBOBOX, OnRefSel)
                    Indx[refSel.GetId()] = i
                    refAtmSizer.Add(refSel,0,WACV)
                refHelpInfo = '''
* The "Orientation Reference" control defines the Cartesian
axes for rigid bodies with the three atoms, A, B and C. 
The vector from B to A defines the x-axis and the y axis is placed 
in the plane defined by B to A and C to A. A,B,C must not be collinear.
'''
                hlp = G2G.HelpButton(VectorRBDisplay,refHelpInfo,wrap=400)
                refAtmSizer.Add(hlp,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,2)
            return refAtmSizer
                        
        def rbVectMag(rbid,imag,rbData):
            
            def OnRBVectorMag(event):
                event.Skip()
                Obj = event.GetEventObject()
                rbid,imag = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    if val <= 0.:
                        raise ValueError
                    rbData['VectMag'][imag] = val
                except ValueError:
                    pass
                Obj.SetValue('%8.4f'%(val))
                wx.CallAfter(UpdateVectorRB,VectorRB.GetScrollPos(wx.VERTICAL))
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,data['Vector'][rbid],plotDefaults)
                
            def OnRBVectorRef(event):
                Obj = event.GetEventObject()
                rbid,imag = Indx[Obj.GetId()]
                rbData['VectRef'][imag] = Obj.GetValue()
                        
            magSizer = wx.BoxSizer(wx.HORIZONTAL)
            magSizer.Add(wx.StaticText(VectorRBDisplay,-1,'Translation magnitude: '),0,WACV)
            magValue = wx.TextCtrl(VectorRBDisplay,-1,'%8.4f'%(rbData['VectMag'][imag]),
                                       style=wx.TE_PROCESS_ENTER)
            Indx[magValue.GetId()] = [rbid,imag]
            magValue.Bind(wx.EVT_TEXT_ENTER,OnRBVectorMag)
            magValue.Bind(wx.EVT_KILL_FOCUS,OnRBVectorMag)
            magSizer.Add(magValue,0,WACV)
            magSizer.Add((5,0),)
            magref = wx.CheckBox(VectorRBDisplay,label=' Refine?') 
            magref.SetValue(rbData['VectRef'][imag])
            magref.Bind(wx.EVT_CHECKBOX,OnRBVectorRef)
            Indx[magref.GetId()] = [rbid,imag]
            magSizer.Add(magref,0,WACV)
            return magSizer
            
        def rbVectors(rbid,imag,mag,XYZ,rbData):

            def TypeSelect(event):
                AtInfo = data['Vector']['AtInfo']
                r,c = event.GetRow(),event.GetCol()
                if vecGrid.GetColLabelValue(c) == 'Type':
                    PE = G2elemGUI.PickElement(G2frame,oneOnly=True)
                    if PE.ShowModal() == wx.ID_OK:
                        if PE.Elem != 'None':
                            El = PE.Elem.strip().lower().capitalize()
                            if El not in AtInfo:
                                Info = G2elem.GetAtomInfo(El)
                                AtInfo[El] = [Info['Drad'],Info['Color']]
                            rbData['rbTypes'][r] = El
                            vecGrid.SetCellValue(r,c,El)
                    PE.Destroy()
                wx.CallAfter(UpdateVectorRB,VectorRB.GetScrollPos(wx.VERTICAL))

            def ChangeCell(event):
                r,c =  event.GetRow(),event.GetCol()
                if r >= 0 and (0 <= c < 3):
                    try:
                        val = float(vecGrid.GetCellValue(r,c))
                        rbData['rbVect'][imag][r][c] = val
                    except ValueError:
                        pass
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,data['Vector'][rbid],plotDefaults)
                wx.CallAfter(UpdateVectorRB,VectorRB.GetScrollPos(wx.VERTICAL))

            vecSizer = wx.BoxSizer()
            Types = 3*[wg.GRID_VALUE_FLOAT+':10,5',]+[wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,5',]
            colLabels = ['Vector x','Vector y','Vector z','Type','Cart x','Cart y','Cart z']
            table = []
            rowLabels = []
            atNames = []
            for ivec,xyz in enumerate(rbData['rbVect'][imag]):
                table.append(list(xyz)+[rbData['rbTypes'][ivec],]+list(XYZ[ivec]))
                rowLabels.append(str(ivec))
                atNames.append(rbData['rbTypes'][ivec]+str(ivec))
            rbData['atNames'] = atNames
            vecTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            vecGrid = G2G.GSGrid(VectorRBDisplay)
            vecGrid.SetTable(vecTable, True)
            if 'phoenix' in wx.version():
                vecGrid.Bind(wg.EVT_GRID_CELL_CHANGED, ChangeCell)
            else:
                vecGrid.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeCell)
            if not imag:
                vecGrid.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, TypeSelect)
            attr = wx.grid.GridCellAttr()
            attr.IncRef()
            attr.SetEditor(G2G.GridFractionEditor(vecGrid))
            for c in range(3):
                attr.IncRef()
                vecGrid.SetColAttr(c, attr)
            for row in range(vecTable.GetNumberRows()):
                if imag:
                    vecGrid.SetCellStyle(row,3,VERY_LIGHT_GREY,True)                    
                for col in [4,5,6]:
                    vecGrid.SetCellStyle(row,col,VERY_LIGHT_GREY,True)
#            vecGrid.SetScrollRate(0,0)
            vecGrid.AutoSizeColumns(False)
            vecSizer.Add(vecGrid)
            return vecSizer
        
        def FillRefChoice(rbid,rbData):
            choiceIds = [i for i in range(len(rbData['rbTypes']))]
            
            rbRef = rbData.get('rbRef',[-1,-1,-1,False])
            for i in range(3):
                if rbRef[i] in choiceIds: choiceIds.remove(rbRef[i])
            refChoice[rbid] = [choiceIds[:],choiceIds[:],choiceIds[:]]
            for i in range(3):
                refChoice[rbid][i].append(rbRef[i])
                refChoice[rbid][i].sort()
                
        def OnRBSelect(event):
            global resRBsel
            sel = rbSelect.GetSelection()
            if sel == 0: return # 1st entry is blank
            rbname = rbchoice[sel-1]
            resRBsel = RBnames[rbname]
            wx.CallLater(100,UpdateVectorRB)
            
        # beginning of UpdateVectorRB
        AtInfo = data['Vector']['AtInfo']
        refChoice = {}
        #RefObjs = []

        GS = VectorRBDisplay.GetSizer()
        if GS: 
            try:        #get around a c++ error in wx 4.0; doing is again seems to be OK
                GS.Clear(True)
            except:
                GS.Clear(True)
        
        RBnames = {}
        for rbid in data['RBIds']['Vector']:
            RBnames.update({data['Vector'][rbid]['RBname']:rbid,})
        if not RBnames:
            return
        rbchoice = list(RBnames.keys())
        rbchoice.sort()
        if GS: 
            VectorRBSizer = GS
        else:
            VectorRBSizer = wx.BoxSizer(wx.VERTICAL)
        if resRBsel not in data['RBIds']['Vector']:
            resRBsel = RBnames[rbchoice[0]]
        if len(RBnames) > 1:
            selSizer = wx.BoxSizer(wx.HORIZONTAL)
            selSizer.Add(wx.StaticText(VectorRBDisplay,label=' Select rigid body to view:'),0)
            rbSelect = wx.ComboBox(VectorRBDisplay,choices=['']+rbchoice)
            name = data['Vector'][resRBsel]['RBname']
            rbSelect.SetSelection(1+rbchoice.index(name))
            rbSelect.Bind(wx.EVT_COMBOBOX,OnRBSelect)
            selSizer.Add(rbSelect,0)
            VectorRBSizer.Add(selSizer,0)
        rbData = data['Vector'][resRBsel]
        if 'DELETED' in str(G2frame.GetStatusBar()):   #seems to be no other way to do this (wx bug)
            if GSASIIpath.GetConfigValue('debug'):
                print ('DBG_wx error: Rigid Body/Status not cleanly deleted after Refine')
            return
        SetStatusLine(' You may use e.g. "c60" or "s60" for a vector entry')
        FillRefChoice(resRBsel,rbData)
        VectorRBSizer.Add(rbNameSizer(resRBsel,rbData),0,wx.EXPAND)
        VectorRBSizer.Add(rbRefAtmSizer(resRBsel,rbData),0)
        XYZ = np.array([[0.,0.,0.] for Ty in rbData['rbTypes']])
        for imag,mag in enumerate(rbData['VectMag']):
            XYZ += mag*rbData['rbVect'][imag]
            VectorRBSizer.Add(rbVectMag(rbid,imag,rbData),0)
            VectorRBSizer.Add(rbVectors(rbid,imag,mag,XYZ,rbData),0)
        VectorRBSizer.Add((5,5),0)
        data['Vector'][rbid]['rbXYZ'] = XYZ        
        
        VectorRBSizer.Add((5,25),)
        VectorRBSizer.Layout()    
        VectorRBDisplay.SetSizer(VectorRBSizer,True)
        VectorRBDisplay.SetAutoLayout(True)
        Size = VectorRBSizer.GetMinSize()
        VectorRBDisplay.SetSize(Size)
       
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        VectorRB.SetSize(Size)
        VectorRB.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
        G2frame.dataWindow.SendSizeEvent()
        
        VectorRBDisplay.Show()
        
    def UpdateSpinRB():
        
        def OnTypeSel(event):
            Obj = event.GetEventObject()
            ObjId = event.GetId()
            data['Spin'][Indx[ObjId]]['rbType'] = Obj.GetValue()
            
        def OnAtSel(event):
            Obj = event.GetEventObject()
            ObjId = event.GetId()
            PE = G2elemGUI.PickElement(G2frame,oneOnly=False)
            if PE.ShowModal() == wx.ID_OK:
                if PE.Elem != 'None':
                    El = PE.Elem.strip().lower().capitalize()
                    data['Spin'][Indx[ObjId]]['atType'] = El
                    data['Spin'][Indx[ObjId]]['Color'] = G2elem.GetAtomInfo(El)['Color']
                    Obj.ChangeValue(El)
                    if 'Q' in El:
                        wx.CallAfter(UpdateSpinRB)
                    
        def OnElSel(event):
            Obj = event.GetEventObject()
            ObjId = event.GetId()
            PE = G2elemGUI.PickElement(G2frame,oneOnly=False,ifOrbs=True)
            if PE.ShowModal() == wx.ID_OK:
                if PE.Elem != 'None':
                    El = PE.Elem.strip().lower().capitalize()
                    data['Spin'][Indx[ObjId]]['elType'] = El
                    Obj.ChangeValue(El)
                   
        def OnSymSel(event):
            ObjId = event.GetId()
            data['Spin'][Indx[ObjId]]['RBsym'] = simsel.GetValue()
            
        #patch
        if 'Spin' not in data:
            data['Spin'] = {}
            data['RBIds'].update({'Spin':[]})
        # end patch
        if not len(data['Spin']):
            return
        GS = SpinRBDisplay.GetSizer()
        if GS: 
            try:        #get around a c++ error in wx 4.0; doing is again seems to be OK
                GS.Clear(True)
            except:
                GS.Clear(True)
        
        SetStatusLine(' ')
        
        if GS: 
            SpinRBSizer = GS
        else:
            SpinRBSizer = wx.BoxSizer(wx.VERTICAL)
        Indx = {}
        SpinRBSizer.Add(wx.StaticText(SpinRBDisplay,label=' Spinning rigid body shells/nonspherical atoms (Atom=Q & select Orbitals):'))
        nQ = 0
        for spinID in data['Spin']:
            if 'Q' in data['Spin'][spinID]['atType']:
                nQ += 1
        if nQ:
            bodSizer = wx.FlexGridSizer(0,6,5,5)
        else:
            bodSizer = wx.FlexGridSizer(0,5,5,5)
        for item in ['Name','Type','RB sym','Atom','Number']:
            bodSizer.Add(wx.StaticText(SpinRBDisplay,label=item))
        for ibod,spinID in enumerate(data['Spin']):
            if nQ:
                bodSizer.Add(wx.StaticText(SpinRBDisplay,label='Orbitals from'))
            bodSizer.Add(G2G.ValidatedTxtCtrl(SpinRBDisplay,data['Spin'][spinID],'RBname'))
            bodSizer.Add(wx.StaticText(SpinRBDisplay,label='Q'),0)
            data['Spin'][spinID]['rbType'] = 'Q'    #patch
            symchoice = ['53m','m3m','-43m','6/mmm','-6m2','-3m','3m','32','4/mmm','-42m','mmm','2/m','-1','1']
            data['Spin'][spinID]['RBsym'] = data['Spin'][spinID].get('RBsym','53m')
            simsel = wx.ComboBox(SpinRBDisplay,choices=symchoice,value=data['Spin'][spinID]['RBsym'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[simsel.GetId()] = spinID
            simsel.Bind(wx.EVT_COMBOBOX,OnSymSel)
            bodSizer.Add(simsel)
            atSel = wx.TextCtrl(SpinRBDisplay,value=data['Spin'][spinID]['atType'],style=wx.TE_PROCESS_ENTER)
            atSel.Bind(wx.EVT_TEXT_ENTER,OnAtSel)
            Indx[atSel.GetId()] = spinID
            bodSizer.Add(atSel,0)
            bodSizer.Add(G2G.ValidatedTxtCtrl(SpinRBDisplay,data['Spin'][spinID],'Natoms'))
            if 'Q' in data['Spin'][spinID]['atType']:
                data['Spin'][spinID]['elType'] = data['Spin'][spinID].get('elType','C')
                elSel = wx.TextCtrl(SpinRBDisplay,value=data['Spin'][spinID]['elType'],style=wx.TE_PROCESS_ENTER)
                elSel.Bind(wx.EVT_TEXT_ENTER,OnElSel)
                Indx[elSel.GetId()] = spinID
                bodSizer.Add(elSel,0)
            elif nQ:
                bodSizer.Add((5,5))
        
        SpinRBSizer.Add(bodSizer)
        SpinRBSizer.Add((5,25),)
        SpinRBSizer.Layout()    
        SpinRBDisplay.SetSizer(SpinRBSizer,True)
        SpinRBDisplay.SetAutoLayout(True)
        Size = SpinRBSizer.GetMinSize()
        SpinRBDisplay.SetSize(Size)
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        SpinRB.SetSize(Size)
        SpinRB.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
        G2frame.dataWindow.SendSizeEvent()
        
        SpinRBDisplay.Show()
        
    def UpdateResidueRB():
        '''Draw the contents of the Residue Rigid Body tab for Rigid Bodies tree entry
        '''
        global resRBsel
        def rbNameSizer(rbid,rbData):
            
            def OnRBName(event):
                event.Skip()
                Obj = event.GetEventObject()
                name = Obj.GetValue()
                if name == rbData['RBname']: return # no change
                namelist = [data['Residue'][key]['RBname'] for key in
                        data['Residue'] if 'RBname' in data['Residue'][key]]
                name = G2obj.MakeUniqueLabel(name,namelist)
                rbData['RBname'] = name
                wx.CallAfter(UpdateResidueRB)
                
            def OnDelRB(event):
                Obj = event.GetEventObject()
                rbid = Indx[Obj.GetId()]
                if rbid in data['Residue']: 
                    del data['Residue'][rbid]
                    data['RBIds']['Residue'].remove(rbid)
                wx.CallAfter(UpdateResidueRB)
                
            def OnStripH(event):
                Obj = event.GetEventObject()
                rbid = Indx[Obj.GetId()]
                if rbid in data['Residue']:
                    newNames = []
                    newTypes = []
                    newXYZ = []
                    for i,atype in enumerate(rbData['rbTypes']):
                        if atype != 'H':
                            newNames.append(rbData['atNames'][i])
                            newTypes.append(rbData['rbTypes'][i])
                            newXYZ.append(rbData['rbXYZ'][i])
                    rbData['atNames'] = newNames
                    rbData['rbTypes'] = newTypes
                    rbData['rbXYZ'] = newXYZ
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                wx.CallAfter(UpdateResidueRB)

            def OnPlotRB(event):
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                
            # start of rbNameSizer
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(ResidueRBDisplay,-1,'Residue name: '),0,WACV)
            RBname = wx.TextCtrl(ResidueRBDisplay,-1,rbData['RBname'],style=wx.TE_PROCESS_ENTER)
            RBname.Bind(wx.EVT_LEAVE_WINDOW, OnRBName)
            RBname.Bind(wx.EVT_TEXT_ENTER,OnRBName)
            RBname.Bind(wx.EVT_KILL_FOCUS,OnRBName)
            nameSizer.Add(RBname,0,WACV)
            nameSizer.Add((5,0),)
            plotRB =  wx.Button(ResidueRBDisplay,wx.ID_ANY,'Plot',style=wx.BU_EXACTFIT)
            plotRB.Bind(wx.EVT_BUTTON, OnPlotRB)
            Indx[plotRB.GetId()] = rbid
            nameSizer.Add(plotRB,0,WACV)
            nameSizer.Add((5,0),)
            if not rbData['useCount']:
                delRB = wx.Button(ResidueRBDisplay,wx.ID_ANY,"Delete",style=wx.BU_EXACTFIT)
                delRB.Bind(wx.EVT_BUTTON, OnDelRB)
                Indx[delRB.GetId()] = rbid
                nameSizer.Add(delRB,0,WACV)
                if 'H'  in rbData['rbTypes']:
                    stripH = wx.Button(ResidueRBDisplay,wx.ID_ANY,"Strip H-atoms",style=wx.BU_EXACTFIT)
                    stripH.Bind(wx.EVT_BUTTON, OnStripH)
                    Indx[stripH.GetId()] = rbid
                    nameSizer.Add(stripH,0,WACV)
            nameSizer.Add(wx.StaticText(ResidueRBDisplay,-1,'  body type #'+
                                        str(data['RBIds']['Residue'].index(rbid))),0,WACV)
            nameSizer.Add((-1,-1),1,wx.EXPAND,1)
            nameSizer.Add(G2G.HelpButton(ResidueRBDisplay,helpIndex=G2frame.dataWindow.helpKey))
            return nameSizer
            
        def rbResidues(rbid,rbData):
            
            def TypeSelect(event):
                AtInfo = data['Residue']['AtInfo']
                r,c = event.GetRow(),event.GetCol()
                if resGrid.GetColLabelValue(c) == 'Type':
                    PE = G2elemGUI.PickElement(G2frame,oneOnly=True)
                    if PE.ShowModal() == wx.ID_OK:
                        if PE.Elem != 'None':
                            El = PE.Elem.strip().lower().capitalize()
                            if El not in AtInfo:
                                Info = G2elem.GetAtomInfo(El)
                                AtInfo[El] = [Info['Drad'],Info['Color']]
                            rbData['rbTypes'][r] = El
                            resGrid.SetCellValue(r,c,El)
                    PE.Destroy()

            def ChangeCell(event):
                r,c =  event.GetRow(),event.GetCol()
                if c == 0:
                    rbData['atNames'][r] = resGrid.GetCellValue(r,c)
                if r >= 0 and (2 <= c <= 4):
                    try:
                        val = float(resGrid.GetCellValue(r,c))
                        rbData['rbXYZ'][r][c-2] = val
                    except ValueError:
                        pass
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                        
            def OnRefSel(event):
                Obj = event.GetEventObject()
                iref,res,jref = Indx[Obj.GetId()]
                sel = Obj.GetValue()
                ind = atNames.index(sel)
                if rbData['rbTypes'][ind] == 'H':
                    G2G.G2MessageBox(G2frame,'You should not select an H-atom for rigid body orientation')
                rbData['rbRef'][iref] = ind
                FillRefChoice(rbid,rbData)
                for i,ref in enumerate(RefObjs[jref]):
                    ref.SetItems([atNames[j] for j in refChoice[rbid][i]])
                    ref.SetValue(atNames[rbData['rbRef'][i]])                    
                rbXYZ = rbData['rbXYZ']
                if not iref:     #origin change
                    rbXYZ -= rbXYZ[ind]
                Xxyz = rbXYZ[rbData['rbRef'][1]]
                X = Xxyz/np.sqrt(np.sum(Xxyz**2))
                Yxyz = rbXYZ[rbData['rbRef'][2]]
                Y = Yxyz/np.sqrt(np.sum(Yxyz**2))
                Mat = G2mth.getRBTransMat(X,Y)
                rbXYZ = np.inner(Mat,rbXYZ).T
                rbData['rbXYZ'] = rbXYZ
                rbData['molCent'] = False
                res.ClearSelection()
                resTable = res.GetTable()
                for r in range(res.GetNumberRows()):
                    row = resTable.GetRowValues(r)
                    row[2:4] = rbXYZ[r]
                    resTable.SetRowValues(r,row)
                res.ForceRefresh()
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                
            def OnMolCent(event):
                rbData['molCent'] = not rbData['molCent']
                if rbData['molCent']:
                    Obj = event.GetEventObject()
                    res = Indx[Obj.GetId()]
                    rbXYZ = rbData['rbXYZ']
                    rbCent = np.array([np.sum(rbXYZ[:,0]),np.sum(rbXYZ[:,1]),np.sum(rbXYZ[:,2])])/rbXYZ.shape[0]
                    rbXYZ -= rbCent
                    rbData['rbXYZ'] = rbXYZ
                    res.ClearSelection()
                    resTable = res.GetTable()
                    for r in range(res.GetNumberRows()):
                        row = resTable.GetRowValues(r)
                        row[2:4] = rbXYZ[r]
                        resTable.SetRowValues(r,row)
                    res.ForceRefresh()
                    G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                    
            def OnHDswitch(event):
                Obj = event.GetEventObject()
                res = Indx[Obj.GetId()]
                resTable = res.GetTable()
                AtInfo = data['Residue']['AtInfo']
                for r in range(res.GetNumberRows()):
                    row = resTable.GetRowValues(r)
                    if row[1] == 'H':
                        row[0] = 'D'+row[0][1:]
                        row[1] = 'D'
                    elif row[1] == 'D':
                        row[0] = 'H'+row[0][1:]
                        row[1] = 'H'                        
                    rbData['atNames'][r] = row[0]
                    rbData['rbTypes'][r] = row[1]
                    Info = G2elem.GetAtomInfo(row[1])
                    AtInfo[row[1]] = [Info['Drad'],Info['Color']]
                    resTable.SetRowValues(r,row)
                res.ForceRefresh()
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                    
            def OnCycleXYZ(event):
                Obj = event.GetEventObject()
                res = Indx[Obj.GetId()]
                rbXYZ = rbData['rbXYZ']
                resTable = res.GetTable()
                for r in range(res.GetNumberRows()):
                    rbXYZ[r] = np.roll(rbXYZ[r],1)
                    row = resTable.GetRowValues(r)
                    row[2:4] = rbXYZ[r]
                    resTable.SetRowValues(r,row)
                res.ForceRefresh()
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
                
            Types = 2*[wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,5',]
            colLabels = ['Name','Type','Cart x','Cart y','Cart z']
            table = []
            rowLabels = []
            for ivec,xyz in enumerate(rbData['rbXYZ']):
                table.append([rbData['atNames'][ivec],]+[rbData['rbTypes'][ivec],]+list(xyz))
                rowLabels.append(str(ivec))
            vecTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            resGrid = G2G.GSGrid(ResidueRBDisplay)
            Indx[resGrid.GetId()] = rbid
            resList.append(resGrid)
            resGrid.SetTable(vecTable, True)
            if 'phoenix' in wx.version():
                resGrid.Bind(wg.EVT_GRID_CELL_CHANGED, ChangeCell)
            else:
                resGrid.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeCell)
            resGrid.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, TypeSelect)
            for c in range(2,5):
                attr = wx.grid.GridCellAttr()
                attr.IncRef()
                attr.SetEditor(G2G.GridFractionEditor(resGrid))
                resGrid.SetColAttr(c, attr)
            for row in range(vecTable.GetNumberRows()):
                resGrid.SetReadOnly(row,1,True)
            resGrid.AutoSizeColumns(False)
            vecSizer = wx.BoxSizer()
            vecSizer.Add(resGrid)
            
            refAtmSizer = wx.BoxSizer(wx.HORIZONTAL)
            atNames = rbData['atNames']
            rbRef = rbData['rbRef']
            refHelpInfo = '''
* The "Orientation Reference" control defines the Cartesian
axes for rigid bodies with the three atoms, A, B and C. 
The vector from B to A defines the x-axis and the y axis is placed 
in the plane defined by B to A and C to A. A,B,C must not be collinear.
 
%%* The origin is at A unless the "Center RB?" button is pressed.

%%* The 'Cycle XYZ' button will permute the rigid body XYZ coordinates so
XYZ --> ZXY. Repeat if needed.

%%* The "Center RB?" button will shift the origin of the 
rigid body to be the midpoint of all atoms in the body (not mass weighted).
'''
            hlp = G2G.HelpButton(ResidueRBDisplay,refHelpInfo,wrap=400)
            refAtmSizer.Add(hlp,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,2)
            refAtmSizer2 = None
            if rbData['rbRef'][3] or rbData['useCount']:
                refAtmSizer.Add(wx.StaticText(ResidueRBDisplay,-1,
                    'Orientation reference non-H atoms A-B-C: %s, %s, %s'%(atNames[rbRef[0]], \
                     atNames[rbRef[1]],atNames[rbRef[2]])),0)
            else:
                refAtmSizer.Add(wx.StaticText(ResidueRBDisplay,-1,
                    'Orientation reference non-H atoms A-B-C: '),0,WACV)
                refObj = [0,0,0]
                for i in range(3):
                    choices = [atNames[j] for j in refChoice[rbid][i]]
                    refSel = wx.ComboBox(ResidueRBDisplay,-1,value='',
                        choices=choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    refSel.SetValue(atNames[rbRef[i]])
                    refSel.Bind(wx.EVT_COMBOBOX, OnRefSel)
                    Indx[refSel.GetId()] = [i,resGrid,len(RefObjs)]
                    refObj[i] = refSel
                    refAtmSizer.Add(refSel,0,WACV)
                refAtmSizer.Add((50,-1))   # moves the help button out a bit
                RefObjs.append(refObj)
                refAtmSizer2 = wx.BoxSizer(wx.HORIZONTAL)
                cycleXYZ = wx.Button(ResidueRBDisplay,label=' Cycle XYZ?')
                cycleXYZ.Bind(wx.EVT_BUTTON,OnCycleXYZ)
                Indx[cycleXYZ.GetId()] = resGrid
                refAtmSizer2.Add(cycleXYZ,0,WACV)
                if 'molCent' not in rbData: rbData['molCent'] = False           #patch
                molcent = wx.Button(ResidueRBDisplay,label=' Center RB?')
                molcent.Bind(wx.EVT_BUTTON,OnMolCent)
                Indx[molcent.GetId()] = resGrid
                refAtmSizer2.Add(molcent,0,WACV)
                HDbut = wx.Button(ResidueRBDisplay,label='H <-> D?')
                HDbut.Bind(wx.EVT_BUTTON,OnHDswitch)
                Indx[HDbut.GetId()] = resGrid
                refAtmSizer2.Add(HDbut,0,WACV)

            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add(refAtmSizer,0,wx.EXPAND)
            if refAtmSizer2: mainSizer.Add(refAtmSizer2)
            mainSizer.Add(vecSizer)
            return mainSizer
            
        def Add2SeqSizer(seqSizer,angSlide,rbid,iSeq,Seq,atNames):
            
            def ChangeAngle(event):
                event.Skip()
                Obj = event.GetEventObject()
                rbid,Seq = Indx[Obj.GetId()][:2]
                val = Seq[2]
                try:
                    val = float(Obj.GetValue())
                    Seq[2] = val
                except ValueError:
                    pass
                Obj.SetValue('%8.2f'%(val))
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,data['Residue'][rbid],plotDefaults)
                
            def OnRadBtn(event):
                Obj = event.GetEventObject()
                Seq,iSeq,angId = Indx[Obj.GetId()]
                data['Residue'][rbid]['SelSeq'] = [iSeq,angId]
                angSlide.SetValue(int(100*Seq[2]))
                
            def OnDelBtn(event):
                Obj = event.GetEventObject()
                rbid,Seq = Indx[Obj.GetId()]
                data['Residue'][rbid]['rbSeq'].remove(Seq)        
                wx.CallAfter(UpdateResidueRB)
            
            iBeg,iFin,angle,iMove = Seq
            ang = wx.TextCtrl(ResidueRBDisplay,wx.ID_ANY,
                    '%8.2f'%(angle),size=(70,-1),style=wx.TE_PROCESS_ENTER)
            if not iSeq:
                radBt = wx.RadioButton(ResidueRBDisplay,wx.ID_ANY,
                                           '',style=wx.RB_GROUP)
                data['Residue'][rbid]['SelSeq'] = [iSeq,ang.GetId()]
                radBt.SetValue(True)
            else:
                radBt = wx.RadioButton(ResidueRBDisplay,wx.ID_ANY,'')
            radBt.Bind(wx.EVT_RADIOBUTTON,OnRadBtn)                   
            seqSizer.Add(radBt)
            delBt =  wx.Button(ResidueRBDisplay,wx.ID_ANY,'Del',
                                style=wx.BU_EXACTFIT)
            delBt.Bind(wx.EVT_BUTTON,OnDelBtn)
            seqSizer.Add(delBt)
            bond = wx.StaticText(ResidueRBDisplay,wx.ID_ANY,
                        '%s %s'%(atNames[iBeg],atNames[iFin]),size=(50,20))
            seqSizer.Add(bond,0,WACV)
            Indx[radBt.GetId()] = [Seq,iSeq,ang.GetId()]
            Indx[delBt.GetId()] = [rbid,Seq]
            Indx[ang.GetId()] = [rbid,Seq,ang]
            ang.Bind(wx.EVT_TEXT_ENTER,ChangeAngle)
            ang.Bind(wx.EVT_KILL_FOCUS,ChangeAngle)
            seqSizer.Add(ang,0,WACV)
            atms = ''
            for i in iMove:    
                atms += ' %s,'%(atNames[i])
            moves = wx.StaticText(ResidueRBDisplay,wx.ID_ANY,
                            atms[:-1],size=(200,20))
            seqSizer.Add(moves,1,wx.EXPAND|wx.RIGHT)
            return seqSizer
            
        def SlideSizer():
            
            def OnSlider(event):
                Obj = event.GetEventObject()
                rbData = Indx[Obj.GetId()]
                iSeq,angId = rbData['SelSeq']
                val = float(Obj.GetValue())/100.
                rbData['rbSeq'][iSeq][2] = val
                Indx[angId][2].SetValue('%8.2f'%(val))
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
            
            slideSizer = wx.BoxSizer(wx.HORIZONTAL)
            slideSizer.Add(wx.StaticText(ResidueRBDisplay,-1,'Selected torsion angle:'),0)
            iSeq,angId = rbData['SelSeq']
            angSlide = G2G.G2Slider(ResidueRBDisplay,-1,
                int(100*rbData['rbSeq'][iSeq][2]),0,36000,size=(200,20),
                style=wx.SL_HORIZONTAL)
            angSlide.Bind(wx.EVT_SLIDER, OnSlider)
            Indx[angSlide.GetId()] = rbData
            slideSizer.Add(angSlide,0)            
            return slideSizer,angSlide
            
        def FillRefChoice(rbid,rbData):
            choiceIds = [i for i in range(len(rbData['atNames']))]
            for seq in rbData['rbSeq']:
                for i in seq[3]:
                    try:
                        choiceIds.remove(i)
                    except ValueError:
                        pass
            rbRef = rbData['rbRef']
            for i in range(3):
                try:
                    choiceIds.remove(rbRef[i])
                except ValueError:
                    pass
            refChoice[rbid] = [choiceIds[:],choiceIds[:],choiceIds[:]]
            for i in range(3):
                refChoice[rbid][i].append(rbRef[i])
                refChoice[rbid][i].sort()
                
        def OnRBSelect(event):
            global resRBsel
            sel = rbSelect.GetSelection()
            if sel == 0: return # 1st entry is blank
            rbname = rbchoice[sel-1]
            resRBsel = RBnames[rbname]
            rbData = data['Residue'][resRBsel]
            G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
            wx.CallLater(100,UpdateResidueRB)
            
        #----- beginning of UpdateResidueRB -----------------------------------------------
        AtInfo = data['Residue']['AtInfo']
        refChoice = {}
        RefObjs = []

        GS = ResidueRBDisplay.GetSizer()
        if GS: 
            try:        #get around a c++ error in wx 4.0; doing is again seems to be OK
                GS.Clear(True)
            except:
                GS.Clear(True)
        
        RBnames = {}
        for rbid in data['RBIds']['Residue']:
            RBnames.update({data['Residue'][rbid]['RBname']:rbid,})
        if not RBnames:
            return
        rbchoice = list(RBnames.keys())
        rbchoice.sort()
        if GS: 
            ResidueRBSizer = GS
        else:
            ResidueRBSizer = wx.BoxSizer(wx.VERTICAL)
        if resRBsel not in data['RBIds']['Residue']:
            resRBsel = RBnames[rbchoice[0]]
        if len(RBnames) > 1:
            selSizer = wx.BoxSizer(wx.HORIZONTAL)
            selSizer.Add(wx.StaticText(ResidueRBDisplay,
                                label=' Select residue to view:'),0)
            rbSelect = wx.ComboBox(ResidueRBDisplay,choices=['']+rbchoice)
            name = data['Residue'][resRBsel]['RBname']
            rbSelect.SetSelection(1+rbchoice.index(name))
            rbSelect.Bind(wx.EVT_COMBOBOX,OnRBSelect)
            selSizer.Add(rbSelect,0)
            ResidueRBSizer.Add(selSizer,0)
        rbData = data['Residue'][resRBsel]
        FillRefChoice(resRBsel,rbData)
        ResidueRBSizer.Add(rbNameSizer(resRBsel,rbData),0,wx.EXPAND)
        ResidueRBSizer.Add(rbResidues(resRBsel,rbData),0)
        if len(rbData['rbSeq']):
            ResidueRBSizer.Add((-1,15),0)
            slideSizer,angSlide = SlideSizer()
            seqSizer = wx.FlexGridSizer(0,5,4,8)
            for lbl in 'Sel','','Bond','Angle','Riding atoms':
                seqSizer.Add(wx.StaticText(ResidueRBDisplay,wx.ID_ANY,lbl))
            ResidueRBSizer.Add(seqSizer)
#            for iSeq,Seq in enumerate(rbData['rbSeq']):
#                ResidueRBSizer.Add(SeqSizer(angSlide,resRBsel,iSeq,Seq,rbData['atNames']))
            for iSeq,Seq in enumerate(rbData['rbSeq']):
                Add2SeqSizer(seqSizer,angSlide,resRBsel,iSeq,Seq,rbData['atNames'])
            ResidueRBSizer.Add(slideSizer,)

        ResidueRBSizer.Add((5,25),)
        ResidueRBSizer.Layout()    
        ResidueRBDisplay.SetSizer(ResidueRBSizer,True)
        ResidueRBDisplay.SetAutoLayout(True)
        Size = ResidueRBSizer.GetMinSize()
        ResidueRBDisplay.SetSize(Size)
       
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        ResidueRB.SetSize(Size)
        ResidueRB.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1)) # dataframe already scrolls
        G2frame.dataWindow.SendSizeEvent()
        
        ResidueRBDisplay.Show()
        
    def SetStatusLine(text):
        G2frame.GetStatusBar().SetStatusText(text,1)                                      

    #### UpdateRigidBodies starts here =========
    global resList,resRBsel            
    if not data.get('RBIds') or not data:
        data.update({'Vector':{'AtInfo':{}},'Residue':{'AtInfo':{}},
            'RBIds':{'Vector':[],'Residue':[]}})       #empty/bad dict - fill it
    Indx = {}
    resList = []
    plotDefaults = {'oldxy':[0.,0.],'Quaternion':[0.,0.,0.,1.],'cameraPos':30.,'viewDir':[0,0,1],}
    G2frame.rbBook = G2G.GSNoteBook(parent=G2frame.dataWindow)
    G2frame.dataWindow.GetSizer().Add(G2frame.rbBook,1,wx.ALL|wx.EXPAND)
    VectorRB = wx.ScrolledWindow(G2frame.rbBook)
    VectorRBDisplay = wx.Panel(VectorRB)
    G2frame.rbBook.AddPage(VectorRB,'Vector rigid bodies')
    ResidueRB = wx.ScrolledWindow(G2frame.rbBook)
    ResidueRBDisplay = wx.Panel(ResidueRB)
    G2frame.rbBook.AddPage(ResidueRB,'Residue rigid bodies')
    # vector RBs are not too common, so select Residue as the default when one is present
    SpinRB = wx.ScrolledWindow(G2frame.rbBook)
    SpinRBDisplay = wx.Panel(SpinRB)
    G2frame.rbBook.AddPage(SpinRB,'Spinning rigid bodies')
    if len(data['RBIds']['Residue']) > 0 and len(data['RBIds']['Vector']) == 0:
        G2frame.rbBook.ChangeSelection(1)
        OnPageChanged(None)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodyMenu)
    SetStatusLine('')
    UpdateVectorRB()
    G2frame.rbBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    wx.CallAfter(OnPageChanged,None)
    
def ShowIsoDistortCalc(G2frame,phase=None):
    '''Compute the ISODISTORT mode values from the current coordinates.
    Called in response to the (Phase/Atoms tab) AtomCompute or 
    Constraints/Edit Constr. "Show ISODISTORT modes" menu item, which 
    should be enabled only when Phase['ISODISTORT'] is defined. 
    '''
    def _onClose(event):
        dlg.EndModal(wx.ID_CANCEL)

    Phases = G2frame.GetPhaseData()
    isophases = [p for p in Phases if 'G2VarList' in Phases[p]['ISODISTORT']]
    if not isophases:
        G2G.G2MessageBox(G2frame,'no ISODISTORT mode data for any phase')
        return
    if phase and phase not in isophases:
        G2G.G2MessageBox(G2frame,'no ISODISTORT mode data for this phase')
        return
    elif not phase and len(isophases) == 1:
        phase = isophases[0]
    elif not phase:
        dlg = wx.SingleChoiceDialog(G2frame,'Select phase from ISODISTORT phases',
            'Select Phase',isophases)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            phase = isophases[sel]
        else:
            return

    covdata = G2frame.GPXtree.GetItemPyData(
        G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Covariance'))
    # make a lookup table for named NewVar Phase constraints
    sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints') 
    Constraints = G2frame.GPXtree.GetItemPyData(sub)
    constrDict = {}
    for c in Constraints['Phase']:
        if c[-1] != 'f' or not c[-3]: continue
        constrDict[str(c[-3])] = c

    dlg = wx.Dialog(G2frame,wx.ID_ANY,'ISODISTORT mode values',#size=(630,400),
        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    data = Phases[phase]
    ISO = data['ISODISTORT']
    mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
        'ISODISTORT mode computation for coordinates in phase '+str(data['General'].get('Name'))))
    aSizer = wx.BoxSizer(wx.HORIZONTAL)
    panel1 = wxscroll.ScrolledPanel(
        dlg, wx.ID_ANY,#size=(100,200),
        style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
    subSizer1 = wx.FlexGridSizer(cols=3,hgap=5,vgap=2)
    panel2 = wxscroll.ScrolledPanel(
        dlg, wx.ID_ANY,#size=(100,200),
        style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
    subSizer2 = wx.FlexGridSizer(cols=4,hgap=5,vgap=2)
    subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,'ISODISTORT\nname'))
    subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,'GSAS-II\nname'))
    subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,' value'),0,wx.ALIGN_RIGHT)
    for i in range(3): subSizer1.Add((-1,5)) # spacer
    subSizer2.Add((-1,-1))
    subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'ISODISTORT\nMode name'))
    subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'GSAS-II\nname'))
    subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,' value'),0,wx.ALIGN_RIGHT)
    for i in range(4): subSizer2.Add((-1,5))
    # ISODISTORT displacive modes
    if 'G2VarList' in ISO:
        dispVals,dispSUs,modeVals,modeSUs = G2mth.CalcIsoDisp(Phases[phase],covdata=covdata)
        for (lbl,xyz,xyzsig,G2var,
             var,mval,msig,G2mode) in zip(
                ISO['IsoVarList'],dispVals,dispSUs,ISO['G2VarList'],
                ISO['IsoModeList'],modeVals,modeSUs,ISO['G2ModeList'] ):
            if str(G2mode) in constrDict:
                ch = G2G.HelpButton(panel2,fmtHelp(constrDict[str(G2mode)],var))
                subSizer2.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
            else:
                subSizer2.Add((-1,-1))
            subSizer1.Add(wx.StaticText(panel1,label=str(lbl)))
            subSizer1.Add(wx.StaticText(panel1,label=str(G2var)))
            try:
                value = G2mth.ValEsd(xyz,xyzsig)
            except TypeError:
                value = str(xyz)            
            subSizer1.Add(wx.StaticText(panel1,label=value),0,wx.ALIGN_RIGHT)

            subSizer2.Add(wx.StaticText(panel2,label=str(var)))
            subSizer2.Add(wx.StaticText(panel2,label=str(G2mode)))
            try:
                # value = G2mth.ValEsd(mval,msig)
                value = '%.5f'%(mval/2.)      #why /2.0
            except TypeError:
                value = str(mval)
            subSizer2.Add(wx.StaticText(panel2,label=value),0,wx.ALIGN_RIGHT)
    # ISODISTORT occupancy modes
    if 'G2OccVarList' in ISO:
        deltaList = []
        parmDict,varyList = G2frame.MakeLSParmDict()
        for gv,Ilbl in zip(ISO['G2OccVarList'],ISO['OccVarList']):
            var = gv.varname()
            albl = Ilbl[:Ilbl.rfind('_')]
            pval = ISO['BaseOcc'][albl]
            if var in parmDict:
                cval = parmDict[var][0]
            else:
                dlg.EndModal(wx.ID_CANCEL)
                G2frame.ErrorDialog('Atom not found',"No value found for parameter "+str(var))
                return
            deltaList.append(cval-pval)
        modeVals = np.inner(ISO['Var2OccMatrix'],deltaList)
        for lbl,delocc,var,val,norm,G2mode in zip(
                ISO['OccVarList'],deltaList,
                ISO['OccModeList'],modeVals,ISO['OccNormList'],ISO['G2OccModeList']):
            if str(G2mode) in constrDict:
                ch = G2G.HelpButton(panel2,fmtHelp(constrDict[str(G2mode)],var))
                subSizer2.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
            else:
                subSizer2.Add((-1,-1))
            subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,str(lbl)))
            try:
                value = G2fl.FormatSigFigs(delocc)
            except TypeError:
                value = str(delocc)
            subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,value),0,wx.ALIGN_RIGHT)
            #subSizer.Add((10,-1))
            subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,str(var)))
            try:
                value = G2fl.FormatSigFigs(val/norm)
                if 'varyList' in covdata:
                    if str(G2mode) in covdata['varyList']:
                        sig = covdata['sig'][covdata['varyList'].index(str(G2mode))]
                        value = G2mth.ValEsd(val/norm,sig/norm)
            except TypeError:
                value = '?'
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
    dlg.CenterOnParent()
    dlg.ShowModal()
    dlg.Destroy()

def ShowIsoModes(G2frame,phase):
    '''Show details about the ISODISTORT mode and the displacements they 
    translate to.
    '''
    def _onClose(event):
        dlg.EndModal(wx.ID_CANCEL)

    # make a lookup table for named NewVar Phase constraints
    sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints') 
    Constraints = G2frame.GPXtree.GetItemPyData(sub)
    constrDict = {}
    for c in Constraints['Phase']:
        if c[-1] != 'f' or not c[-3]: continue
        constrDict[str(c[-3])] = c
        
    Phases = G2frame.GetPhaseData()
    data = Phases[phase]
    ISO = data['ISODISTORT']
    if 'IsoVarList' in ISO:
        modeExp = {}
        for i,row in enumerate(ISO['Var2ModeMatrix']):
            line = '('
            l = 1
            for j,k in enumerate(row):
                var = ISO['IsoVarList'][j]
                if np.isclose(k,0): continue
                if l > 30:
                    line += '\n'
                    l = 0
                l += 3
                if k < 0 and j > 0:
                    line += ' - '
                    k = -k
                elif j > 0: 
                    line += ' + '
                if np.isclose(k,1):
                    l1 = '%s ' % str(var)
                else:
                    l1 = '%.3f * %s' % (k,str(var))
                line += l1
                l += len(l1)
            line += ') / {:.3g}'.format(ISO['NormList'][i])
            modeExp[str(ISO['G2ModeList'][i])] = line

        crdExp = {}
        for i,(lbl,row) in enumerate(zip(ISO['IsoVarList'],ISO['Mode2VarMatrix'])):
            l = ''
            for j,(k,n) in enumerate(zip(row,ISO['NormList'])):
                if np.isclose(k,0): continue
                l1 = ''
                if j > 0 and k < 0:
                    k = -k
                    l1 = '-'
                elif j > 0:
                    l1 += ' +'
                if np.isclose(k,1):
                    l += '{:} {:4g} * {:}'.format(l1, n, ISO['G2ModeList'][j])
                else:
                    l += '{:} {:3g} * {:4g} * {:}'.format(l1, k, n, ISO['G2ModeList'][j])               
            crdExp[lbl] = l

    dlg = wx.Dialog(G2frame,wx.ID_ANY,'ISODISTORT modes and displacements',#size=(630,400),
        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
        'ISODISTORT modes and displacements in phase '+str(data['General'].get('Name','?'))))
    # ISODISTORT displacive modes
    if 'G2VarList' in ISO:
        panel1 = wxscroll.ScrolledPanel(dlg, style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        subSizer1 = wx.FlexGridSizer(cols=5,hgap=5,vgap=2)
        panel2 = wxscroll.ScrolledPanel(dlg, style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        subSizer2 = wx.FlexGridSizer(cols=6,hgap=5,vgap=2)
        subSizer1.Add(wx.StaticText(panel1,label='ISODISTORT\nvariable'),0,wx.ALIGN_CENTER)
        subSizer1.Add((8,-1)) # spacer
        subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,'GSAS-II\nequiv.'),0,wx.ALIGN_CENTER)
        subSizer1.Add((8,-1)) # spacer
        subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,'expression'),0,wx.ALIGN_CENTER)
        subSizer2.Add((-1,-1))
        subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'ISODISTORT\nMode name'),0,wx.ALIGN_CENTER)
        subSizer2.Add((8,-1)) # spacer
        subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'GSAS-II\nequiv.'),0,wx.ALIGN_CENTER)
        subSizer2.Add((8,-1)) # spacer
        subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,'expression'),0,wx.ALIGN_CENTER)
        for i in range(6): subSizer2.Add((-1,5))
        for i,lbl in enumerate(ISO['IsoVarList']):
            var = str(ISO['G2VarList'][i]).replace('::dA','::A')
            if np.isclose(ISO['G2coordOffset'][i],0):
                G2var = '{:}'.format(var)
            elif ISO['G2coordOffset'][i] < 0:
                G2var = '({:} + {:.3g})'.format(var,-ISO['G2coordOffset'][i])
            else:
                G2var = '({:} - {:.3g})'.format(var,ISO['G2coordOffset'][i])
            subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,str(lbl)),0,WACV)
            subSizer1.Add((-1,-1)) # spacer
            subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,str(G2var)),0,WACV)
            subSizer1.Add((-1,-1)) # spacer
            subSizer1.Add(wx.StaticText(panel1,wx.ID_ANY,crdExp[lbl]),0,WACV)
            
        for (isomode,G2mode) in zip(ISO['IsoModeList'],ISO['G2ModeList']):
            if str(G2mode) in constrDict:
                ch = G2G.HelpButton(panel2,fmtHelp(constrDict[str(G2mode)],isomode))
                subSizer2.Add(ch,0,wx.LEFT|wx.RIGHT|WACV|wx.ALIGN_CENTER,1)
            else:
                subSizer2.Add((-1,-1))
            subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,str(isomode)))
            subSizer2.Add((-1,-1)) # spacer
            subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,str(G2mode)))
            subSizer2.Add((-1,-1)) # spacer
            subSizer2.Add(wx.StaticText(panel2,wx.ID_ANY,modeExp[str(G2mode)]))
            
        # finish up ScrolledPanel
        panel1.SetSizer(subSizer1)
        panel2.SetSizer(subSizer2)
        panel1.SetAutoLayout(1)
        panel1.SetupScrolling()
        panel2.SetAutoLayout(1)
        panel2.SetupScrolling()
        # Allow window to be enlarged but not made smaller
        w1,l1 = subSizer1.GetSize()
        w2,l2 = subSizer2.GetSize()
        panel1.SetMinSize((min(700,w1+20),max(50,l1)))
        panel2.SetMinSize((min(700,w2+20),max(50,l2)))
        mainSizer.Add(panel1,1, wx.ALL|wx.EXPAND,1)
        mainSizer.Add(panel2,1, wx.ALL|wx.EXPAND,1)

    # make OK button 
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(dlg, wx.ID_CLOSE) 
    btn.Bind(wx.EVT_BUTTON,_onClose)
    btnsizer.Add(btn)
    mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)

    dlg.SetSizer(mainSizer)
    mainSizer.Fit(dlg)
    dlg.SetMinSize(dlg.GetSize())
    dlg.CenterOnParent()
    dlg.ShowModal()
    dlg.Destroy()

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
