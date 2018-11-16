# -*- coding: utf-8 -*-
#GSASIIconstrGUI - constraint GUI routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIconstrGUI: Constraint GUI routines*
------------------------------------------

Used to define constraints and rigid bodies.

'''
from __future__ import division, print_function
import sys
import wx
import wx.grid as wg
import random as ran
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI
import GSASIIstrIO as G2stIO
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIlattice as G2lat
import GSASIIdataGUI as G2gd
import GSASIIctrlGUI as G2G
import GSASIIplot as G2plt
import GSASIIobj as G2obj
import GSASIIspc as G2spc
VERY_LIGHT_GREY = wx.Colour(235,235,235)

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
        mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
        mainSizer.Add((10,10),1)
        dataGridSizer = wx.FlexGridSizer(cols=3,hgap=2,vgap=2)
        self.OkBtn = wx.Button(panel,wx.ID_OK)
        for id in range(len(self.data)):
            lbl1 = lbl = str(self.data[id][1])
            if lbl[-1] != '=': lbl1 = lbl + ' ' + separator + ' '
            name = wx.StaticText(panel,wx.ID_ANY,lbl1,
                                 style=wx.ALIGN_RIGHT)
            scale = G2G.ValidatedTxtCtrl(panel,self.data[id],0,
                                          typeHint=float,
                                          OKcontrol=self.DisableOK)
            dataGridSizer.Add(name,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,5)
            dataGridSizer.Add(scale,0,wx.RIGHT,3)
            if ':' in lbl:
                dataGridSizer.Add(
                    wx.StaticText(panel,-1,G2obj.fmtVarDescr(lbl)),
                    0,wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,3)
            else:
                dataGridSizer.Add((-1,-1))
        if title == 'New Variable':
            name = wx.StaticText(panel,wx.ID_ANY,"New variable's\nname (optional)",
                                 style=wx.ALIGN_CENTER)
            scale = G2G.ValidatedTxtCtrl(panel,self.newvar,0,
                                          typeHint=str,notBlank=False)
            dataGridSizer.Add(name,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,5)
            dataGridSizer.Add(scale,0,wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,3)
            self.refine = wx.CheckBox(panel,label='Refine?')
            self.refine.SetValue(self.newvar[1]==True)
            self.refine.Bind(wx.EVT_CHECKBOX, self.OnCheckBox)
            dataGridSizer.Add(self.refine,0,wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,3)
            
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

        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()
        self.Fit()
        self.CenterOnParent()
        
    def DisableOK(self,setting):
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
        
################################################################################
#####  Constraints
################################################################################           
       
def UpdateConstraints(G2frame,data):
    '''Called when Constraints tree item is selected.
    Displays the constraints in the data window
    '''
    if not data:
        data.update({'Hist':[],'HAP':[],'Phase':[],'Global':[]})       #empty dict - fill it
    if 'Global' not in data:                                            #patch
        data['Global'] = []
    # DEBUG code ########################################
    #import GSASIIconstrGUI
    #reload(GSASIIconstrGUI)
    #reload(G2obj)
    #reload(G2stIO)
    #import GSASIIstrMain
    #reload(GSASIIstrMain)    
    #reload(G2mv)
    #reload(G2gd)
    ###################################################
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    if not len(Phases) or not len(Histograms):
        dlg = wx.MessageDialog(G2frame,'You need both phases and histograms to see Constraints',
            'No phases or histograms')
        dlg.CenterOnParent()
        dlg.ShowModal()
        dlg.Destroy()
        return
    G2obj.IndexAllIds(Histograms,Phases)
    ##################################################################################
    # patch: convert old-style (str) variables in constraints to G2VarObj objects
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
    ##################################################################################
    rigidbodyDict = G2frame.GPXtree.GetItemPyData(
        G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
    badPhaseParms = ['Ax','Ay','Az','Amul','AI/A','Atype','SHorder','mV0','mV1','mV2','waveType','Vol','isMag',]
    globalList = list(rbDict.keys())
    globalList.sort()
    try:
        AtomDict = dict([Phases[phase]['pId'],Phases[phase]['Atoms']] for phase in Phases)
    except KeyError:
        G2frame.ErrorDialog('Constraint Error','Constraints cannot be set until a cycle of least squares'+
                            ' has been run.\nWe suggest you refine a scale factor.')
        return

    # create a list of the phase variables
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable,MFtable,maxSSwave = G2stIO.GetPhaseData(Phases,rbIds=rbIds,Print=False)
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
    seqList = G2frame.testSeqRefineMode()
    if seqList: # for sequential refinement, only process 1st histgram in list
        histDict = {seqList[0]:Histograms[seqList[0]]}
    else:
        histDict = Histograms
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,histDict,Print=False,resetRefList=False)
    hapList = sorted([i for i in hapDict.keys() if i.split(':')[2] not in ('Type',)])
    if seqList: # convert histogram # to wildcard
        wildList = [] # list of variables with "*" for histogram number
        for i in hapList:
            s = i.split(':')
            if s[1] == "": continue
            s[1] = '*'
            sj = ':'.join(s)
            if sj not in wildList: wildList.append(sj)
        hapList = wildList
    histVary,histDict,controlDict = G2stIO.GetHistogramData(histDict,Print=False)
    histList = list(histDict.keys())
    histList.sort()
    if seqList: # convert histogram # to wildcard
        wildList = [] # list of variables with "*" for histogram number
        for i in histList:
            s = i.split(':')
            if s[1] == "": continue
            s[1] = '*'
            sj = ':'.join(s)
            if sj not in wildList: wildList.append(sj)
        histList = wildList        
    Indx = {}
    G2frame.Page = [0,'phs']
        
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
        elif 'AM' in name:
            namelist = ['AMx','AMy','AMz']
        elif items[-1] in ['A0','A1','A2','A3','A4','A5']:
            namelist = ['A0','A1','A2','A3','A4','A5']
        elif items[-1] in ['D11','D22','D33','D12','D13','D23']:
            namelist = ['D11','D22','D33','D12','D13','D23']
        elif 'Tm' in name:
            namelist = ['Tmin','Tmax']
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

        for item in nameList:
            keys = item.split(':')
            #if keys[0] not in phlist:
            #    phlist.append(keys[0])
            if items[1] == '*' and keys[2] in namelist: # wildcard -- select only sequential options
                keys[1] = '*'
                mitem = ':'.join(keys)
                if mitem == name: continue
                if mitem not in outList: outList.append(mitem)
            elif keys[2] in namelist and item != name:
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
    
    def FindAllCons(data):
        ''' Find all constraints
        '''
        allcons = []
        for key in data:
            if key.startswith('_'): continue
            allcons += data[key]
        return allcons
        
    def CheckConstraints(constraintSet):
        '''Check for errors in a set of constraints. Constraints based on symmetry (etc.)
        are generated by running :func:`GSASIIstrIO.GetPhaseData`.
        '''
        G2mv.InitVars()    
        constDictList,fixedList,ignored = G2stIO.ProcessConstraints(constraintSet)
        # generate symmetry constraints to check for conflicts
        rigidbodyDict = G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,maxSSwave = G2stIO.GetPhaseData(
            Phases,RestraintDict=None,rbIds=rbIds,Print=False) # generates atom symmetry constraints
        msg = G2mv.EvaluateMultipliers(constDictList,phaseDict)
        if msg: return 'Unable to interpret multiplier(s): '+msg
        return G2mv.CheckConstraints('',constDictList,fixedList)

    def CheckAddedConstraint(newcons):
        '''Check a new constraint that has just been input.
        If there is an error display a message and discard the last entry

        Since the varylist is not available, no warning messages
        should be generated here

        :returns: True if constraint should be added
        '''
        
        allcons1 = FindAllCons(data)
        allcons = allcons1[:]
        allcons += newcons
        if not len(allcons): return True
        errmsg,warnmsg = CheckConstraints(allcons)
        if errmsg:
            G2frame.ErrorDialog('Constraint Error',
                'Error with newly added constraint:\n'+errmsg+
                '\nIgnoring newly added constraint',parent=G2frame)
            # Note no multiplier formulas in GUI, skipping EvaluateMultipliers  
            # reset error status
            errmsg,warnmsg = CheckConstraints(allcons1)
            if errmsg:
                print (errmsg)
                print (G2mv.VarRemapShow([],True))
            return False
        elif warnmsg:
            print ('Unexpected contraint warning:\n'+warnmsg)
        return True

    def CheckChangedConstraint():
        '''Check all constraints after an edit has been made.
        If there is an error display a message and reject the change.

        Since the varylist is not available, no warning messages
        should be generated.
        
        :returns: True if the edit should be retained
        '''
        allcons = FindAllCons(data)
        if not len(allcons): return True
        # Note no multiplier formulas in GUI, skipping EvaluateMultipliers  
        errmsg,warnmsg = CheckConstraints(allcons)
        if errmsg:
            G2frame.ErrorDialog('Constraint Error',
                'Error after editing constraint:\n'+errmsg+
                '\nDiscarding last constraint edit',parent=G2frame)
            # reset error status
            errmsg,warnmsg = CheckConstraints(allcons)
            if errmsg:
                print (errmsg)
                print (G2mv.VarRemapShow([],True))
            return False
        elif warnmsg:
            print ('Unexpected contraint warning:\n'+warnmsg)
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
        wx.CallAfter(OnPageChanged,None)
        
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
        ''' Add riding equivalences between all parameters on atoms '''
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
        wx.CallAfter(OnPageChanged,None)
                        
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
        wx.CallAfter(OnPageChanged,None)
                        
    def MakeConstraintsSizer(name,pageDisplay):
        '''Creates a sizer displaying all of the constraints entered of
        the specified type.

        :param str name: the type of constraints to be displayed ('HAP',
          'Hist', 'Phase', 'Global', 'Sym-Generated')
        :param wx.Panel pageDisplay: parent panel for sizer
        :returns: wx.Sizer created by method
        '''
        if name == 'Sym-Generated':         #show symmetry generated constraints
            Sizer1 =  wx.BoxSizer(wx.VERTICAL)
            Sizer1.Add(wx.StaticText(pageDisplay,wx.ID_ANY,
                                    'Equivalences generated based on cell/space group input'))
            Sizer1.Add((-1,5))
            Sizer = wx.FlexGridSizer(0,2,0,0)
            Sizer1.Add(Sizer)
            for sym in G2mv.GetSymEquiv():
                Sizer.Add(wx.StaticText(pageDisplay,wx.ID_ANY,'EQUIV'),
                           0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT,3)
                Sizer.Add(wx.StaticText(pageDisplay,wx.ID_ANY,sym))
                Sizer.Add((-1,-1))
                Sizer.Add((-1,2))
            return Sizer1
        constSizer = wx.FlexGridSizer(0,6,0,0)
        maxlen = 50 # characters before wrapping a constraint
        for Id,item in enumerate(data[name]):
            refineflag = False
            helptext = ""
            eqString = ['',]
            problemItem = False
            for term in item[:-3]:
                if str(term[1]) in G2mv.problemVars:
                    problemItem = True
            if item[-1] == 'h': # Hold on variable
                constSizer.Add((-1,-1),0)              # blank space for edit button
                typeString = 'FIXED'
                var = str(item[0][1])
                varMean = G2obj.fmtVarDescr(var)
                eqString[-1] =  var +'   '
                helptext = "Prevents variable:\n"+ var + " ("+ varMean + ")\nfrom being changed"
            elif isinstance(item[-1],str): # not true on original-style (2011?) constraints
                constEdit = wx.Button(pageDisplay,wx.ID_ANY,'Edit',style=wx.BU_EXACTFIT)
                constEdit.Bind(wx.EVT_BUTTON,OnConstEdit)
                constSizer.Add(constEdit,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER,1)            # edit button
                Indx[constEdit.GetId()] = [Id,name]
                if item[-1] == 'f':
                    helptext = "A new variable"
                    if item[-3]:
                        helptext += " named "+str(item[-3])
                    helptext += " is created from a linear combination of the following variables:\n"
                    for term in item[:-3]:
                        var = str(term[1])
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        m = term[0]
                        if eqString[-1] != '':
                            if m >= 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                                m = abs(m)
                        eqString[-1] += '%.3f*%s '%(m,var)
                        varMean = G2obj.fmtVarDescr(var)
                        helptext += "\n" + var + " ("+ varMean + ")"
                    if '_Explain' in data:
                        if data['_Explain'].get(item[-3]):
                            helptext += '\n\n'
                            helptext += data['_Explain'][item[-3]]
                    # typeString = 'NEWVAR'
                    # if item[-3]:
                    #     eqString[-1] += ' = '+item[-3]
                    # else:
                    #     eqString[-1] += ' = New Variable'
                    if item[-3]:
                        typeString = item[-3] + ' = '
                    else:
                        typeString = 'New Variable = '
                    #print 'refine',item[-2]
                    refineflag = True
                elif item[-1] == 'c':
                    helptext = "The following variables constrained to equal a constant:"
                    for term in item[:-3]:
                        var = str(term[1])
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        if eqString[-1] != '':
                            if term[0] > 0:
                                eqString[-1] += ' + '
                            else:
                                eqString[-1] += ' - '
                        eqString[-1] += '%.3f*%s '%(abs(term[0]),var)
                        varMean = G2obj.fmtVarDescr(var)
                        helptext += "\n" + var + " ("+ varMean + ")"
                    typeString = 'CONST'
                    eqString[-1] += ' = '+str(item[-3])
                elif item[-1] == 'e':
                    helptext = "The following variables are set to be equivalent, noting multipliers:"
                    normval = item[:-3][1][0]
                    for i,term in enumerate(item[:-3]):
                        var = str(term[1])
                        if term[0] == 0: term[0] = 1.0
                        if len(eqString[-1]) > maxlen:
                            eqString.append(' ')
                        if i == 0: # move independent variable to end
                            indepterm = term
                            continue
                        elif eqString[-1] != '':
                            eqString[-1] += ' = '
                        if normval/term[0] == 1:
                            eqString[-1] += '%s'% var
                        else:
                            eqString[-1] += '%.3f*%s'%(normval/term[0],var)
                        varMean = G2obj.fmtVarDescr(var)
                        helptext += "\n" + var + " ("+ varMean + ")"
                    if normval/indepterm[0] == 1:
                        eqString[-1] += ' = %s'% str(indepterm[1])
                    else:
                        eqString[-1] += ' = %.3f*%s'%(normval/indepterm[0],str(indepterm[1]))
                    typeString = 'EQUIV'
                else:
                    print ('Unexpected constraint'+item)
                
            else:
                print ('Removing old-style constraints')
                data[name] = []
                return constSizer
            constDel = wx.Button(pageDisplay,wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
            constDel.Bind(wx.EVT_BUTTON,OnConstDel)
            Indx[constDel.GetId()] = [Id,name]
            constSizer.Add(constDel,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER,1)             # delete button
            if helptext:
                ch = G2G.HelpButton(pageDisplay,helptext)
                constSizer.Add(ch,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER,1)
            else:
                constSizer.Add((-1,-1))
            if refineflag:
                ch = G2G.G2CheckBox(pageDisplay,'',item,-2)
                constSizer.Add(ch,0,wx.LEFT|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER,1)
            else:
                constSizer.Add((-1,-1))                
            constSizer.Add(wx.StaticText(pageDisplay,wx.ID_ANY,typeString),
                           0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER|wx.RIGHT|wx.LEFT,3)
            if problemItem: eqString[-1] += ' -- Conflict: see console'
            if len(eqString) > 1:
                Eq = wx.BoxSizer(wx.VERTICAL)
                for s in eqString:
                    line = wx.StaticText(pageDisplay,wx.ID_ANY,s)
                    if problemItem: line.SetBackgroundColour(wx.YELLOW)
                    Eq.Add(line,0,wx.ALIGN_CENTER_VERTICAL)
                Eq.Add((-1,4))
            else:
                Eq = wx.StaticText(pageDisplay,wx.ID_ANY,eqString[0])
                if problemItem: Eq.SetBackgroundColour(wx.YELLOW)
            constSizer.Add(Eq,1,wx.ALIGN_CENTER_VERTICAL)
        return constSizer
                
    def OnConstDel(event):
        'Delete a constraint'
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        del data[name][Id]
        allcons = FindAllCons(data)     #should I call CheckChangedConstraint() instead?
        if not len(allcons): return
        CheckConstraints(allcons)
        wx.CallAfter(OnPageChanged,None)
        
    def OnConstEdit(event):
        '''Called to edit an individual contraint in response to a
        click on its Edit button
        '''
        Obj = event.GetEventObject()
        Id,name = Indx[Obj.GetId()]
        if data[name][Id][-1] == 'f':
            items = data[name][Id][:-3]
            constType = 'New Variable'
            if data[name][Id][-3]:
                varname = data[name][Id][-3]
            else:
                varname = ""
            lbl = 'Enter value for each term in constraint; sum = new variable'
            dlg = ConstraintDialog(G2frame,constType,lbl,items,
                                   varname=varname,varyflag=data[name][Id][-2])
        elif data[name][Id][-1] == 'c':
            items = data[name][Id][:-3]+[
                [data[name][Id][-3],'fixed value =']]
            constType = 'Constraint'
            lbl = 'Edit value for each term in constant constraint sum'
            dlg = ConstraintDialog(G2frame,constType,lbl,items)
        elif data[name][Id][-1] == 'e':
            items = data[name][Id][:-3]
            constType = 'Equivalence'
            lbl = 'The following terms are set to be equal:'
            dlg = ConstraintDialog(G2frame,constType,lbl,items,'/')
        else:
            return
        try:
            prev = data[name][Id][:]
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetData()
                for i in range(len(data[name][Id][:-3])):
                    if type(data[name][Id][i]) is tuple: # fix non-mutable construct
                        data[name][Id][i] = list(data[name][Id][i])
                    data[name][Id][i][0] = result[i][0]
                if data[name][Id][-1] == 'c':
                    data[name][Id][-3] = str(result[-1][0])
                elif data[name][Id][-1] == 'f':
                    # process the variable name to put in global form (::var)
                    varname = str(dlg.newvar[0]).strip().replace(' ','_')
                    if varname.startswith('::'):
                        varname = varname[2:]
                    varname = varname.replace(':',';')
                    if varname:
                        data[name][Id][-3] = varname
                    else:
                        data[name][Id][-3] = ''
                    data[name][Id][-2] = dlg.newvar[1]
                if not CheckChangedConstraint():
                    data[name][Id] = prev
        except:
            import traceback
            print (traceback.format_exc())
        finally:
            dlg.Destroy()
        wx.CallAfter(OnPageChanged,None)
   
    def UpdateConstraintPanel(panel,typ):
        '''Update the contents of the selected Constraint
        notebook tab. Called in :func:`OnPageChanged`
        '''
        if panel.GetSizer(): panel.GetSizer().Clear(True)
        Siz = wx.BoxSizer(wx.VERTICAL)
        Siz.Add((5,5),0)
        Siz.Add(MakeConstraintsSizer(typ,panel),1,wx.EXPAND)
        panel.SetSizer(Siz,True)
        Size = Siz.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        panel.SetSize(Size)
        panel.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        panel.Show()

    def OnPageChanged(event):
        '''Called when a tab is pressed or when a "select tab" menu button is
        used (see RaisePage), or to refresh the current tab contents (event=None)
        '''
        if event:       #page change event!
            page = event.GetSelection()
        else: # called directly, get current page
            page = G2frame.constr.GetSelection()
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
        
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.ConstraintMenu)
    SetStatusLine('')
    
    G2frame.Bind(wx.EVT_MENU, OnAddConstraint, id=G2G.wxID_CONSTRAINTADD)
    G2frame.Bind(wx.EVT_MENU, OnAddFunction, id=G2G.wxID_FUNCTADD)
    G2frame.Bind(wx.EVT_MENU, OnAddEquivalence, id=G2G.wxID_EQUIVADD)
    G2frame.Bind(wx.EVT_MENU, OnAddHold, id=G2G.wxID_HOLDADD)
    G2frame.Bind(wx.EVT_MENU, OnAddAtomEquiv, id=G2G.wxID_EQUIVALANCEATOMS)
#    G2frame.Bind(wx.EVT_MENU, OnAddRiding, id=G2G.wxID_ADDRIDING)
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
    wx.CallAfter(OnPageChanged,None)
    G2frame.constr.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    # validate all the constrants -- should not see any errors here normally
    allcons = FindAllCons(data)
    if not len(allcons): return
    errmsg,warnmsg = CheckConstraints(allcons)
    if errmsg:
        G2frame.ErrorDialog('Constraint Error',
                            'Error in constraints:\n'+errmsg+'\nCheck console output for more information',
                            parent=G2frame)
        print (errmsg)
        print (G2mv.VarRemapShow([],True))
    elif warnmsg:
        print ('Unexpected contraint warning:\n'+warnmsg)

################################################################################
# check scale & phase fractions, create constraint if needed
################################################################################
def CheckAllScalePhaseFractions(G2frame):
    '''Check if scale factor and all phase fractions are refined without a constraint
    for all used histograms, if so, offer the user a chance to create a constraint
    on the sum of phase fractions
    '''
    histograms, phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    for i,hist in enumerate(histograms):
        CheckScalePhaseFractions(G2frame,hist,histograms,phases)
        
def CheckScalePhaseFractions(G2frame,hist,histograms,phases):
    '''Check if scale factor and all phase fractions are refined without a constraint
    for histogram hist, if so, offer the user a chance to create a constraint
    on the sum of phase fractions
    '''
    if G2frame.testSeqRefineMode():
        histStr = '*'
    else:
        histStr = str(histograms[hist]['hId'])
    # Is this powder? 
    if not hist.startswith('PWDR '): return
    # do this only if the scale factor is varied
    if not histograms[hist]['Sample Parameters']['Scale'][1]: return
    # are all phase fractions varied in all used histograms?
    phaseCount = 0
    for p in phases:
        if hist not in phases[p]['Histograms']: continue
        if phases[p]['Histograms'][hist]['Use'] and not phases[p]['Histograms'][hist]['Scale'][1]:
            return
        else:
            phaseCount += 1
    
    # all phase fractions and scale factor varied, now scan through constraints
    sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints') 
    Constraints = G2frame.GPXtree.GetItemPyData(sub)
    for c in Constraints.get('HAP',[]):
        if c[-1] != 'c': continue
        if not c[-3]: continue
        if len(c[:-3]) != phaseCount: continue
        # got a constraint equation with right number of terms, is it on phase fractions for
        # the correct histogram?
        if all([(i[1].name == 'Scale' and i[1].varname().split(':')[1] == histStr) for i in c[:-3]]):
            # got a constraint, this is OK
            return
    dlg = wx.MessageDialog(G2frame,'You are refining the scale factor and all phase fractions for histogram #'+
        histStr+'. This will produce an unstable refinement. '+
        'Do you want to constrain the sum of phase fractions?','Create constraint?',wx.OK|wx.CANCEL)
    if dlg.ShowModal() != wx.ID_OK:
        dlg.Destroy()
        return
    dlg.Destroy()

    constr = []
    for p in phases:
        if hist not in phases[p]['Histograms']: continue
        if not phases[p]['Histograms'][hist]['Use']: continue
        constr += [[1.0,G2obj.G2VarObj(':'.join((str(phases[p]['pId']),histStr,'Scale')))]]
    constr += [1.0,None,'c']
    Constraints['HAP'] += [constr]
        
################################################################################
#### Make nuclear/magnetic phase transition constraints - called by OnTransform in G2phsGUI
################################################################################        
        
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
            elif iA == 3 and SGUniq == 'a':
                parm = '%d::A%s'%(pId,iA)
            elif iA == 4 and SGUniq == 'b':
                parm = '%d::A%s'%(pId,iA)
            elif iA == 5 and SGUniq == 'c':
                parm = '%d::A%s'%(pId,iA)
            else:
                parm = None            
        else:
            parm = '%d::A%s'%(pId,iA)
        return parm
    
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    UseList = newPhase['Histograms']
    detTrans = np.abs(nl.det(Trans))
    nAcof = G2lat.cell2A(newPhase['General']['Cell'][1:7])
    
    opId = oldPhase['pId']
    npId = newPhase['pId']
    cx,ct,cs,cia = newPhase['General']['AtomPtrs']
    nAtoms = newPhase['Atoms']
    oSGData = oldPhase['General']['SGData']
    nSGData = newPhase['General']['SGData']
    oAcof = G2lat.cell2A(oldPhase['General']['Cell'][1:7])
    nAcof = G2lat.cell2A(newPhase['General']['Cell'][1:7])
    item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints')
    if not item:
        return
    constraints = G2frame.GPXtree.GetItemPyData(item)
    parmDict = {}
    varyList = []
    xnames = ['dAx','dAy','dAz']
#    invTrans = nl.inv(Trans)
#    Us = ['AU11','AU22','AU33','AU12','AU13','AU23']
#    Uids = [[0,0,'AU11'],[1,1,'AU22'],[2,2,'AU33'],[0,1,'AU12'],[0,2,'AU13'],[1,2,'AU23']]
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
        for ix in list(set(CSX[0])):
            if not ix:
                continue
            name = xnames[ix-1]
            IndpCon = [1.0,G2obj.G2VarObj('%d::%s:%d'%(npId,name,ia))]
            DepCons = []
            for iop,opval in enumerate(XOpr[ix-1]):
                if opval:
                    DepCons.append([opval,G2obj.G2VarObj('%d::%s:%s'%(opId,xnames[iop],iat))])
            if len(DepCons) == 1:
                constraints['Phase'].append([DepCons[0],IndpCon,None,None,'e'])
            elif len(DepCons) > 1:
                for Dep in DepCons:
                    Dep[0] *= -1
                constraints['Phase'].append([IndpCon]+DepCons+[0.0,None,'c'])
        for name in ['Afrac','AUiso']:
            IndpCon = [1.0,G2obj.G2VarObj('%d::%s:%d'%(npId,name,ia))]
            DepCons = [1.0,G2obj.G2VarObj('%d::%s:%s'%(opId,name,iat))]
            constraints['Phase'].append([DepCons,IndpCon,None,None,'e'])
        
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
        
#unfortunately, this doesn't always work!
#    As = ['A0','A1','A2','A3','A4','A5']
#    Aids = [[0,0,'A0'],[1,1,'A1'],[2,2,'A2'],[0,1,'A3'],[0,2,'A4'],[1,2,'A5']]
#    DepConsDict = dict(zip(As,[[],[],[],[],[],[]]))
#    T = Trans
##Symbolic code:
#    '''
#         T00**2*a0  T01**2*a1 T02**2*a2 T00*T01*a3    T00*T02*a4    T01*T02*a5 
#         T10**2*a0  T11**2*a1 T12**2*a2 T10*T11*a3    T10*T12*a4    T11*T12*a5 
#         T20**2*a0  T21**2*a1 T22**2*a2 T20*T21*a3    T20*T22*a4    T21*T22*a5 
#         2*T00*T10*a0      2*T01*T11*a1     2*T02*T12*a2     (T00*T11 + T01*T10)*a3      (T00*T12 + T02*T10)*a4      (T01*T12 + T02*T11)*a5 
#         2*T00*T20*a0      2*T01*T21*a1     2*T02*T22*a2     (T00*T21 + T01*T20)*a3      (T00*T22 + T02*T20)*a4      (T01*T22 + T02*T21)*a5 
#         2*T10*T20*a0      2*T11*T21*a1     2*T12*T22*a2     (T10*T21 + T11*T20)*a3      (T10*T22 + T12*T20)*a4      (T11*T22 + T12*T21)*a5 
#    '''
#    conMat = [
#        [T[0,0]**2,T[0,1]**2,T[0,2]**2,T[0,0]*T[0,1],T[0,0]*T[0,2],T[0,1]*T[0,2]],
#        [T[1,0]**2,T[1,1]**2,T[1,2]**2,T[1,0]*T[1,1],T[1,0]*T[1,2],T[1,1]*T[1,2]],
#        [T[2,0]**2,T[2,1]**2,T[2,2]**2,T[2,0]*T[2,1],T[2,0]*T[2,2],T[2,1]*T[2,2]],
#        [2.*T[0,0]*T[1,0],2.*T[0,1]*T[1,1],2.*T[0,2]*T[1,2],T[0,0]*T[1,1]+T[0,1]*T[1,0],T[0,0]*T[1,2]+T[0,2]*T[1,0],T[0,1]*T[1,2]+T[0,2]*T[1,1]],
#        [2.*T[0,0]*T[2,0],2.*T[0,1]*T[2,1],2.*T[0,2]*T[2,2],T[0,0]*T[2,1]+T[0,1]*T[2,0],T[0,0]*T[2,2]+T[0,2]*T[2,0],T[0,1]*T[2,2]+T[0,2]*T[2,1]],
#        [2.*T[1,0]*T[2,0],2.*T[1,1]*T[2,1],2.*T[1,2]*T[2,2],T[1,0]*T[2,1]+T[1,1]*T[2,0],T[1,0]*T[2,2]+T[1,2]*T[2,0],T[1,1]*T[2,2]+T[1,2]*T[2,1]]]
#    
#    for iA,Aid in enumerate(Aids):
#        if abs(nAcof[iA]) > 1.e-8:
#            for ia in [0,1,2,3,4,5]:
#                cA = conMat[ia][iA]
#                if abs(cA) > 1.e-8:
#                    parm = SetUniqAj(npId,ia,nSGData)
#                    DepConsDict[Aid[2]].append([cA,G2obj.G2VarObj(parm)])
#    conStrings = []
#    for iA,Asi in enumerate(As):
#        parm = SetUniqAj(opId,iA,oSGData)
#        parmDict[parm] = oAcof[iA]
#        varyList.append(parm)
#        IndpCon = [1.0,G2obj.G2VarObj(parm)]
#        conStr = str([IndpCon,DepConsDict[Asi]])
#        if conStr in conStrings:
#            continue
#        conStrings.append(conStr)
#        if len(DepConsDict[Asi]) == 1:
#            if DepConsDict[Asi][0]:
#                constraints['Phase'].append([IndpCon,DepConsDict[Asi][0],None,None,'e'])
#        elif len(DepConsDict[Asi]) > 1:        
#            for Dep in DepConsDict[Asi]:
#                Dep[0] *= -1
#            constraints['Phase'].append([IndpCon]+DepConsDict[Asi]+[0.0,None,'c'])
            
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
        
################################################################################
#### Rigid bodies
################################################################################

def UpdateRigidBodies(G2frame,data):
    '''Called when Rigid bodies tree item is selected.
    Displays the rigid bodies in the data window
    '''
    if not data.get('RBIds') or not data:
        data.update({'Vector':{'AtInfo':{}},'Residue':{'AtInfo':{}},
            'RBIds':{'Vector':[],'Residue':[]}})       #empty/bad dict - fill it
            
    global resList,rbId
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
    
    def OnPageChanged(event):
        global resList
        resList = []
        if event:       #page change event!
            page = event.GetSelection()
        else:
            page = G2frame.rbBook.GetSelection()
        G2frame.rbBook.ChangeSelection(page)
        text = G2frame.rbBook.GetPageText(page)
        if text == 'Vector rigid bodies':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.VectorBodyMenu)
            G2frame.Bind(wx.EVT_MENU, AddVectorRB, id=G2G.wxID_VECTORBODYADD)
            G2frame.Page = [page,'vrb']
            UpdateVectorRB()
        elif text == 'Residue rigid bodies':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodyMenu)
            G2frame.Bind(wx.EVT_MENU, AddResidueRB, id=G2G.wxID_RIGIDBODYADD)
            G2frame.Bind(wx.EVT_MENU, OnImportRigidBody, id=G2G.wxID_RIGIDBODYIMPORT)
            G2frame.Bind(wx.EVT_MENU, OnDefineTorsSeq, id=G2G.wxID_RESIDUETORSSEQ) #enable only if residue RBs exist?
            G2frame.Page = [page,'rrb']
            UpdateResidueRB()
            
    def getMacroFile(macName):
        defDir = os.path.join(os.path.split(__file__)[0],'GSASIImacros')
        dlg = wx.FileDialog(G2frame,message='Choose '+macName+' rigid body macro file',
            defaultDir=defDir,defaultFile="",wildcard="GSAS-II macro file (*.mac)|*.mac",
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                macfile = dlg.GetPath()
                macro = open(macfile,'Ur')
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
        dlg = wx.FileDialog(G2frame,'Choose rigid body text file', '.', '',
            "GSAS-II text file (*.txt)|*.txt|XYZ file (*.xyz)|*.xyz|"
            "Sybyl mol2 file (*.mol2)|*.mol2|PDB file (*.pdb;*.ent)|*.pdb;*.ent",
            wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                txtfile = dlg.GetPath()
                ext = os.path.splitext(txtfile)[1]
                text = open(txtfile,'Ur')
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
            ImportResidueRB()
            
    def AddVectorRB(event):
        AtInfo = data['Vector']['AtInfo']
        dlg = G2G.MultiIntegerDialog(G2frame,'New Rigid Body',['No. atoms','No. translations'],[1,1])
        if dlg.ShowModal() == wx.ID_OK:
            nAtoms,nTrans = dlg.GetValues()
            rbId = ran.randint(0,sys.maxsize)
            vecMag = [1.0 for i in range(nTrans)]
            vecRef = [False for i in range(nTrans)]
            vecVal = [np.zeros((nAtoms,3)) for j in range(nTrans)]
            rbTypes = ['C' for i in range(nAtoms)]
            Info = G2elem.GetAtomInfo('C')
            AtInfo['C'] = [Info['Drad'],Info['Color']]
            data['Vector'][rbId] = {'RBname':'UNKRB','VectMag':vecMag,'rbXYZ':np.zeros((nAtoms,3)),
                'rbRef':[0,1,2,False],'VectRef':vecRef,'rbTypes':rbTypes,'rbVect':vecVal,'useCount':0}
            data['RBIds']['Vector'].append(rbId)
        dlg.Destroy()
        UpdateVectorRB()
        
    def AddResidueRB(event):
        AtInfo = data['Residue']['AtInfo']
        macro = getMacroFile('rigid body')
        if not macro:
            return
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'I' == items[0]:
                rbId = ran.randint(0,sys.maxsize)
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
                data['Residue'][rbId] = {'RBname':rbName,'rbXYZ':rbXYZ,'rbTypes':rbTypes,
                    'atNames':atNames,'rbRef':[nOrig-1,mRef-1,nRef-1,True],'rbSeq':rbSeq,
                    'SelSeq':[0,0],'useCount':0}
                data['RBIds']['Residue'].append(rbId)
                print ('Rigid body '+rbName+' added')
            macStr = macro.readline()
        macro.close()
        UpdateResidueRB()
        
    def ImportResidueRB():
        AtInfo = data['Residue']['AtInfo']
        text,ext = getTextFile()
        if not text:
            return
        rbId = ran.randint(0,sys.maxsize)
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
                #print txtStr
        items = txtStr.split()
        while len(items):
            if 'txt' in ext:
                atName = items[0]
                atType = items[1]
                rbXYZ.append([float(items[i]) for i in [2,3,4]])
            elif 'xyz' in ext:
                atType = items[0]
                rbXYZ.append([float(items[i]) for i in [1,2,3]])
                atName = atType+str(len(rbXYZ))
            elif 'mol2' in ext:
                atType = items[1]
                atName = items[1]+items[0]
                rbXYZ.append([float(items[i]) for i in [2,3,4]])
            elif 'pdb' in ext:
                atType = items[-1]
                atName = items[2]
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
            data['Residue'][rbId] = {'RBname':'UNKRB','rbXYZ':rbXYZ,'rbTypes':rbTypes,
                'atNames':atNames,'rbRef':[0,1,2,False],'rbSeq':[],'SelSeq':[0,0],'useCount':0}
            data['RBIds']['Residue'].append(rbId)
            print ('Rigid body UNKRB added')
        text.close()
        UpdateResidueRB(rbId)
        
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
        global rbId
        rbData = data['Residue'][rbId]
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
        UpdateResidueRB(rbId)

    def UpdateVectorRB(Scroll=0):
        AtInfo = data['Vector']['AtInfo']
        refChoice = {}
        if 'DELETED' in str(G2frame.GetStatusBar()):   #seems to be no other way to do this (wx bug)
            if GSASIIpath.GetConfigValue('debug'):
                print ('DBG_wx error: Rigid Body/Status not cleanly deleted after Refine')
            return
        SetStatusLine(' You may use e.g. "c60" or "s60" for a vector entry')
        def rbNameSizer(rbId,rbData):

            def OnRBName(event):
                event.Skip()
                Obj = event.GetEventObject()
                rbData['RBname'] = Obj.GetValue()
                
            def OnDelRB(event):
                Obj = event.GetEventObject()
                rbId = Indx[Obj.GetId()]
                if rbId in data['Vector']:
                    del data['Vector'][rbId]
                    data['RBIds']['Vector'].remove(rbId)
                    rbData['useCount'] -= 1
                wx.CallAfter(UpdateVectorRB)
                
            def OnPlotRB(event):
                Obj = event.GetEventObject()
                Obj.SetValue(False)
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,rbData,plotDefaults)
            
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(VectorRBDisplay,-1,'Rigid body name: '),
                0,wx.ALIGN_CENTER_VERTICAL)
            RBname = wx.TextCtrl(VectorRBDisplay,-1,rbData['RBname'])
            Indx[RBname.GetId()] = rbId
            RBname.Bind(wx.EVT_TEXT_ENTER,OnRBName)
            RBname.Bind(wx.EVT_KILL_FOCUS,OnRBName)
            nameSizer.Add(RBname,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add((5,0),)
            plotRB = wx.CheckBox(VectorRBDisplay,-1,'Plot?')
            Indx[plotRB.GetId()] = rbId
            plotRB.Bind(wx.EVT_CHECKBOX,OnPlotRB)
            nameSizer.Add(plotRB,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add((5,0),)
            if not rbData['useCount']:
                delRB = wx.CheckBox(VectorRBDisplay,-1,'Delete?')
                Indx[delRB.GetId()] = rbId
                delRB.Bind(wx.EVT_CHECKBOX,OnDelRB)
                nameSizer.Add(delRB,0,wx.ALIGN_CENTER_VERTICAL)
            return nameSizer
            
        def rbRefAtmSizer(rbId,rbData):
            
            def OnRefSel(event):
                Obj = event.GetEventObject()
                iref = Indx[Obj.GetId()]
                sel = Obj.GetValue()
                rbData['rbRef'][iref] = atNames.index(sel)
                FillRefChoice(rbId,rbData)
            
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
                    'Orientation reference atoms A-B-C: '),0,wx.ALIGN_CENTER_VERTICAL)
                for i in range(3):
                    choices = [atNames[j] for j in refChoice[rbId][i]]
                    refSel = wx.ComboBox(VectorRBDisplay,-1,value='',
                        choices=choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    refSel.SetValue(atNames[rbRef[i]])
                    refSel.Bind(wx.EVT_COMBOBOX, OnRefSel)
                    Indx[refSel.GetId()] = i
                    refAtmSizer.Add(refSel,0,wx.ALIGN_CENTER_VERTICAL)
            return refAtmSizer
                        
        def rbVectMag(rbId,imag,rbData):
            
            def OnRBVectorMag(event):
                event.Skip()
                Obj = event.GetEventObject()
                rbId,imag = Indx[Obj.GetId()]
                try:
                    val = float(Obj.GetValue())
                    if val <= 0.:
                        raise ValueError
                    rbData['VectMag'][imag] = val
                except ValueError:
                    pass
                Obj.SetValue('%8.4f'%(val))
                wx.CallAfter(UpdateVectorRB,VectorRB.GetScrollPos(wx.VERTICAL))
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,data['Vector'][rbId],plotDefaults)
                
            def OnRBVectorRef(event):
                Obj = event.GetEventObject()
                rbId,imag = Indx[Obj.GetId()]
                rbData['VectRef'][imag] = Obj.GetValue()
                        
            magSizer = wx.BoxSizer(wx.HORIZONTAL)
            magSizer.Add(wx.StaticText(VectorRBDisplay,-1,'Translation magnitude: '),
                0,wx.ALIGN_CENTER_VERTICAL)
            magValue = wx.TextCtrl(VectorRBDisplay,-1,'%8.4f'%(rbData['VectMag'][imag]))
            Indx[magValue.GetId()] = [rbId,imag]
            magValue.Bind(wx.EVT_TEXT_ENTER,OnRBVectorMag)
            magValue.Bind(wx.EVT_KILL_FOCUS,OnRBVectorMag)
            magSizer.Add(magValue,0,wx.ALIGN_CENTER_VERTICAL)
            magSizer.Add((5,0),)
            magref = wx.CheckBox(VectorRBDisplay,-1,label=' Refine?') 
            magref.SetValue(rbData['VectRef'][imag])
            magref.Bind(wx.EVT_CHECKBOX,OnRBVectorRef)
            Indx[magref.GetId()] = [rbId,imag]
            magSizer.Add(magref,0,wx.ALIGN_CENTER_VERTICAL)
            return magSizer
            
        def rbVectors(rbId,imag,mag,XYZ,rbData):

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
                G2plt.PlotRigidBody(G2frame,'Vector',AtInfo,data['Vector'][rbId],plotDefaults)
                wx.CallAfter(UpdateVectorRB,VectorRB.GetScrollPos(wx.VERTICAL))

            vecSizer = wx.BoxSizer()
            Types = 3*[wg.GRID_VALUE_FLOAT+':10,5',]+[wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,5',]
            colLabels = ['Vector x','Vector y','Vector z','Type','Cart x','Cart y','Cart z']
            table = []
            rowLabels = []
            for ivec,xyz in enumerate(rbData['rbVect'][imag]):
                table.append(list(xyz)+[rbData['rbTypes'][ivec],]+list(XYZ[ivec]))
                rowLabels.append(str(ivec))
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
            attr.SetEditor(G2G.GridFractionEditor(vecGrid))
            for c in range(3):
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
        
        def FillRefChoice(rbId,rbData):
            choiceIds = [i for i in range(len(rbData['rbTypes']))]
            
            rbRef = rbData.get('rbRef',[-1,-1,-1,False])
            for i in range(3):
                choiceIds.remove(rbRef[i])
            refChoice[rbId] = [choiceIds[:],choiceIds[:],choiceIds[:]]
            for i in range(3):
                refChoice[rbId][i].append(rbRef[i])
                refChoice[rbId][i].sort()     
            
        if VectorRB.GetSizer(): VectorRB.GetSizer().Clear(True)
        VectorRBSizer = wx.BoxSizer(wx.VERTICAL)
        for rbId in data['RBIds']['Vector']:
            if rbId != 'AtInfo':
                rbData = data['Vector'][rbId]
                FillRefChoice(rbId,rbData)
                VectorRBSizer.Add(rbNameSizer(rbId,rbData),0)
                VectorRBSizer.Add(rbRefAtmSizer(rbId,rbData),0)
                XYZ = np.array([[0.,0.,0.] for Ty in rbData['rbTypes']])
                for imag,mag in enumerate(rbData['VectMag']):
                    XYZ += mag*rbData['rbVect'][imag]
                    VectorRBSizer.Add(rbVectMag(rbId,imag,rbData),0)
                    VectorRBSizer.Add(rbVectors(rbId,imag,mag,XYZ,rbData),0)
                VectorRBSizer.Add((5,5),0)
                data['Vector'][rbId]['rbXYZ'] = XYZ        
        VectorRBSizer.Layout()    
        VectorRBDisplay.SetSizer(VectorRBSizer,True)
        Size = VectorRBSizer.GetMinSize()
        Size[0] += 40
        Size[1] = max(Size[1],450) + 20
        VectorRBDisplay.SetSize(Size)
        VectorRB.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        VectorRB.Scroll(0,Scroll)
        
    def UpdateResidueRB(rbId=0):
        AtInfo = data['Residue']['AtInfo']
        refChoice = {}
        RefObjs = []

        def rbNameSizer(rbId,rbData):

            def OnRBName(event):
                Obj = event.GetEventObject()
                rbData['RBname'] = Obj.GetValue()
                wx.CallAfter(UpdateResidueRB,rbId)
                
            def OnDelRB(event):
                Obj = event.GetEventObject()
                rbId = Indx[Obj.GetId()]
                if rbId in data['Residue']: 
                    del data['Residue'][rbId]
                    data['RBIds']['Residue'].remove(rbId)
                wx.CallAfter(UpdateResidueRB)
                
            def OnStripH(event):
                Obj = event.GetEventObject()
                rbId = Indx[Obj.GetId()]
                if rbId in data['Residue']:
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
                wx.CallAfter(UpdateResidueRB,rbId)
                    
            def OnPlotRB(event):
                Obj = event.GetEventObject()
                Obj.SetValue(False)
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,rbData,plotDefaults)
            
            nameSizer = wx.BoxSizer(wx.HORIZONTAL)
            nameSizer.Add(wx.StaticText(ResidueRBDisplay,-1,'Residue name: '),
                0,wx.ALIGN_CENTER_VERTICAL)
            RBname = wx.TextCtrl(ResidueRBDisplay,-1,rbData['RBname'])
            Indx[RBname.GetId()] = rbId
            RBname.Bind(wx.EVT_TEXT_ENTER,OnRBName)
            RBname.Bind(wx.EVT_KILL_FOCUS,OnRBName)
            nameSizer.Add(RBname,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add((5,0),)
            plotRB = wx.CheckBox(ResidueRBDisplay,-1,'Plot?')
            Indx[plotRB.GetId()] = rbId
            plotRB.Bind(wx.EVT_CHECKBOX,OnPlotRB)
            nameSizer.Add(plotRB,0,wx.ALIGN_CENTER_VERTICAL)
            nameSizer.Add((5,0),)
            if not rbData['useCount']:
                delRB = wx.CheckBox(ResidueRBDisplay,-1,'Delete?')
                Indx[delRB.GetId()] = rbId
                delRB.Bind(wx.EVT_CHECKBOX,OnDelRB)
                nameSizer.Add(delRB,0,wx.ALIGN_CENTER_VERTICAL)
                if 'H'  in rbData['rbTypes']:
                    stripH = wx.CheckBox(ResidueRBDisplay,-1,'Strip H-atoms?')
                    Indx[stripH.GetId()] = rbId
                    stripH.Bind(wx.EVT_CHECKBOX,OnStripH)
                    nameSizer.Add(stripH,0,wx.ALIGN_CENTER_VERTICAL)
            return nameSizer
            
        def rbResidues(rbId,rbData):
            
            def TypeSelect(event):
                AtInfo = data['Residue']['AtInfo']
                r,c = event.GetRow(),event.GetCol()
                if vecGrid.GetColLabelValue(c) == 'Type':
                    PE = G2elemGUI.PickElement(G2frame,oneOnly=True)
                    if PE.ShowModal() == wx.ID_OK:
                        if PE.Elem != 'None':
                            El = PE.Elem.strip().lower().capitalize()
                            if El not in AtInfo:
                                Info = G2elem.GetAtomInfo(El)
                                AtInfo[El] = [Info['Drad']['Color']]
                            rbData['rbTypes'][r] = El
                            vecGrid.SetCellValue(r,c,El)
                    PE.Destroy()

            def ChangeCell(event):
                r,c =  event.GetRow(),event.GetCol()
                if r >= 0 and (0 <= c < 3):
                    try:
                        val = float(vecGrid.GetCellValue(r,c))
                        rbData['rbXYZ'][r][c] = val
                    except ValueError:
                        pass
                        
            def RowSelect(event):
                r,c =  event.GetRow(),event.GetCol()
                if c < 0:                   #only row clicks
                    for vecgrid in resList:
                        vecgrid.ClearSelection()
                    vecGrid.SelectRow(r,True)

            def OnRefSel(event):
                
                Obj = event.GetEventObject()
                iref,res,jref = Indx[Obj.GetId()]
                sel = Obj.GetValue()
                ind = atNames.index(sel)
                if rbData['rbTypes'][ind] == 'H':
                    G2G.G2MessageBox(G2frame,'You should not select an H-atom for rigid body orientation')
                rbData['rbRef'][iref] = ind
                FillRefChoice(rbId,rbData)
                for i,ref in enumerate(RefObjs[jref]):
                    ref.SetItems([atNames[j] for j in refChoice[rbId][i]])
                    ref.SetValue(atNames[rbData['rbRef'][i]])                    
                rbXYZ = rbData['rbXYZ']
                if not iref:     #origin change
                    rbXYZ -= rbXYZ[ind]
                #TODO - transform all atom XYZ by axis choices
                Xxyz = rbXYZ[rbData['rbRef'][1]]
                X = Xxyz/np.sqrt(np.sum(Xxyz**2))
                Yxyz = rbXYZ[rbData['rbRef'][2]]
                Y = Yxyz/np.sqrt(np.sum(Yxyz**2))
                Mat = G2mth.getRBTransMat(X,Y)
                rbXYZ = np.inner(Mat,rbXYZ).T
                rbData['rbXYZ'] = rbXYZ
                res.ClearSelection()
                resTable = res.GetTable()
                for r in range(res.GetNumberRows()):
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
            vecGrid = G2G.GSGrid(ResidueRBDisplay)
            Indx[vecGrid.GetId()] = rbId
            resList.append(vecGrid)
            vecGrid.SetTable(vecTable, True)
            if 'phoenix' in wx.version():
                vecGrid.Bind(wg.EVT_GRID_CELL_CHANGED, ChangeCell)
            else:
                vecGrid.Bind(wg.EVT_GRID_CELL_CHANGE, ChangeCell)
            vecGrid.Bind(wg.EVT_GRID_CELL_LEFT_DCLICK, TypeSelect)
            vecGrid.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, RowSelect)
            attr = wx.grid.GridCellAttr()
            attr.SetEditor(G2G.GridFractionEditor(vecGrid))
            for c in range(3):
                vecGrid.SetColAttr(c, attr)
            for row in range(vecTable.GetNumberRows()):
                for col in range(5):
                    vecGrid.SetCellStyle(row,col,VERY_LIGHT_GREY,True)
            vecGrid.AutoSizeColumns(False)
            vecSizer = wx.BoxSizer()
            vecSizer.Add(vecGrid)
            
            refAtmSizer = wx.BoxSizer(wx.HORIZONTAL)
            atNames = rbData['atNames']
            rbRef = rbData['rbRef']
            if rbData['rbRef'][3] or rbData['useCount']:
                refAtmSizer.Add(wx.StaticText(ResidueRBDisplay,-1,
                    'Orientation reference non-H atoms A-B-C: %s, %s, %s'%(atNames[rbRef[0]], \
                     atNames[rbRef[1]],atNames[rbRef[2]])),0)
            else:
                refAtmSizer.Add(wx.StaticText(ResidueRBDisplay,-1,
                    'Orientation reference non-H atoms A-B-C: '),0,wx.ALIGN_CENTER_VERTICAL)
                refObj = [0,0,0]
                for i in range(3):
                    choices = [atNames[j] for j in refChoice[rbId][i]]
                    refSel = wx.ComboBox(ResidueRBDisplay,-1,value='',
                        choices=choices,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                    refSel.SetValue(atNames[rbRef[i]])
                    refSel.Bind(wx.EVT_COMBOBOX, OnRefSel)
                    Indx[refSel.GetId()] = [i,vecGrid,len(RefObjs)]
                    refObj[i] = refSel
                    refAtmSizer.Add(refSel,0,wx.ALIGN_CENTER_VERTICAL)
                RefObjs.append(refObj)
            
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add(refAtmSizer)
            mainSizer.Add(vecSizer)
            return mainSizer
            
        def SeqSizer(angSlide,rbId,iSeq,Seq,atNames):
            
            def ChangeAngle(event):
                event.Skip()
                Obj = event.GetEventObject()
                rbId,Seq = Indx[Obj.GetId()][:2]
                val = Seq[2]
                try:
                    val = float(Obj.GetValue())
                    Seq[2] = val
                except ValueError:
                    pass
                Obj.SetValue('%8.2f'%(val))
                G2plt.PlotRigidBody(G2frame,'Residue',AtInfo,data['Residue'][rbId],plotDefaults)
                
            def OnRadBtn(event):
                Obj = event.GetEventObject()
                Seq,iSeq,angId = Indx[Obj.GetId()]
                data['Residue'][rbId]['SelSeq'] = [iSeq,angId]
                angSlide.SetValue(int(100*Seq[2]))
                
            def OnDelBtn(event):
                Obj = event.GetEventObject()
                rbId,Seq = Indx[Obj.GetId()]
                data['Residue'][rbId]['rbSeq'].remove(Seq)        
                wx.CallAfter(UpdateResidueRB,rbId)
            
            seqSizer = wx.FlexGridSizer(0,5,2,2)
            seqSizer.AddGrowableCol(3,0)
            iBeg,iFin,angle,iMove = Seq
            ang = wx.TextCtrl(ResidueRBDisplay,-1,'%8.2f'%(angle),size=(50,20))
            if not iSeq:
                radBt = wx.RadioButton(ResidueRBDisplay,-1,'',style=wx.RB_GROUP)
                data['Residue'][rbId]['SelSeq'] = [iSeq,ang.GetId()]
            else:
                radBt = wx.RadioButton(ResidueRBDisplay,-1,'')
            radBt.Bind(wx.EVT_RADIOBUTTON,OnRadBtn)                   
            seqSizer.Add(radBt)
            delBt = wx.RadioButton(ResidueRBDisplay,-1,'')
            delBt.Bind(wx.EVT_RADIOBUTTON,OnDelBtn)
            seqSizer.Add(delBt)
            bond = wx.TextCtrl(ResidueRBDisplay,-1,'%s %s'%(atNames[iBeg],atNames[iFin]),size=(50,20))
            seqSizer.Add(bond,0,wx.ALIGN_CENTER_VERTICAL)
            Indx[radBt.GetId()] = [Seq,iSeq,ang.GetId()]
            Indx[delBt.GetId()] = [rbId,Seq]
            Indx[ang.GetId()] = [rbId,Seq,ang]
            ang.Bind(wx.EVT_TEXT_ENTER,ChangeAngle)
            ang.Bind(wx.EVT_KILL_FOCUS,ChangeAngle)
            seqSizer.Add(ang,0,wx.ALIGN_CENTER_VERTICAL)
            atms = ''
            for i in iMove:    
                atms += ' %s,'%(atNames[i])
            moves = wx.TextCtrl(ResidueRBDisplay,-1,atms[:-1],size=(200,20))
            seqSizer.Add(moves,1,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.RIGHT)
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
            angSlide = wx.Slider(ResidueRBDisplay,-1,
                int(100*rbData['rbSeq'][iSeq][2]),0,36000,size=(200,20),
                style=wx.SL_HORIZONTAL)
            angSlide.Bind(wx.EVT_SLIDER, OnSlider)
            Indx[angSlide.GetId()] = rbData
            slideSizer.Add(angSlide,0)            
            return slideSizer,angSlide
            
        def FillRefChoice(rbId,rbData):
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
            refChoice[rbId] = [choiceIds[:],choiceIds[:],choiceIds[:]]
            for i in range(3):
                refChoice[rbId][i].append(rbRef[i])
                refChoice[rbId][i].sort()
                
        def OnSelect(event):
            rbname = rbchoice[select.GetSelection()]
            rbId = RBnames[rbname]
            wx.CallLater(100,UpdateResidueRB,rbId)
            
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
        ResidueRBSizer = wx.BoxSizer(wx.VERTICAL)
        if len(RBnames) > 1:
            selSizer = wx.BoxSizer(wx.HORIZONTAL)
            selSizer.Add(wx.StaticText(ResidueRBDisplay,label=' Select residue to view:'),0)
            rbchoice.sort()
            select = wx.ComboBox(ResidueRBDisplay,choices=rbchoice)
            select.Bind(wx.EVT_COMBOBOX,OnSelect)
            selSizer.Add(select,0)
            ResidueRBSizer.Add(selSizer,0)
        if not rbId:
            rbId = RBnames[rbchoice[0]]
        rbData = data['Residue'][rbId]
        FillRefChoice(rbId,rbData)
        ResidueRBSizer.Add(rbNameSizer(rbId,rbData),0)
        ResidueRBSizer.Add(rbResidues(rbId,rbData),0)
        ResidueRBSizer.Add((5,5),0)
        if rbData['rbSeq']:
            slideSizer,angSlide = SlideSizer()
        if len(rbData['rbSeq']):
            ResidueRBSizer.Add(wx.StaticText(ResidueRBDisplay,-1,
                'Sel  Del  Bond             Angle      Riding atoms'),
                0,wx.ALIGN_CENTER_VERTICAL)                       
        for iSeq,Seq in enumerate(rbData['rbSeq']):
            ResidueRBSizer.Add(SeqSizer(angSlide,rbId,iSeq,Seq,rbData['atNames']))
        if rbData['rbSeq']:
            ResidueRBSizer.Add(slideSizer,)

        ResidueRBSizer.Add((5,25),)
        ResidueRBSizer.Layout()    
        ResidueRBDisplay.SetSizer(ResidueRBSizer,True)
        ResidueRBDisplay.SetAutoLayout(True)
        Size = ResidueRBSizer.GetMinSize()
        ResidueRBDisplay.SetSize(Size)
        ResidueRBDisplay.Show()
        
    def SetStatusLine(text):
        G2frame.GetStatusBar().SetStatusText(text,1)                                      

    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RigidBodyMenu)
    SetStatusLine('')
    UpdateVectorRB()
    G2frame.rbBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    wx.CallAfter(OnPageChanged,None)
