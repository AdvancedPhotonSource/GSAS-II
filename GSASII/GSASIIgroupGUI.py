# -*- coding: utf-8 -*-
'''
Routines for working with groups of histograms. 

Groups are defined in Controls entry ['Groups'] which contains three entries:

* Controls['Groups']['groupDict'] 
   a dict where each key is the name of the group and the value is a list of 
   histograms in the group
* Controls['Groups']['notGrouped']
   a count of the number of histograms that are not in any group
* Controls['Groups']['template'] 
   the string used to set the grouping

See SearchGroups in :func:`GSASIIdataGUI.UpdateControls`.


** Parameter Data Table **
Prior to display in the GUI, parameters are organized in a dict where each 
dict entry has contents of form:

 *  'label' : `innerdict`

where `innerdict` can contain the following elements:

 *  'val' : (array, key)
 *  'range' : (float,float)
 *  'ref' : (array, key)
 *  'str' : string
 *  'init' : float
 *  'rowlbl' : (array, key)

One of 'val', 'ref' or 'str' elements will be present. 

 *  The 'str' value is something that cannot be edited; If 'str' is
    present, it will be the only entry in `innerdict`. It is used for
    a parameter value that is typically computed or must be edited in
    the histogram section.

 *  The 'val' tuple provides a reference to the float value for the 
    defined quantity, where array[key] provides r/w access to the
    parameter.

 *  The 'range' list/tuple provides min and max float value for the 
    defined quantity to be defined. Use None for any value that
    should not be enforced. The 'range' values will  

 *  The 'ref' tuple provides a reference to the bool value, where 
    array[key] provides r/w access to the refine flag for the 
    labeled quantity

    Both 'ref' and 'val' are usually defined together, but either may 
    occur alone. These exceptions will be for parameters where a single
    refine flag is used for a group of parameters.

 *  The 'init' value is also something that cannot be edited. 
    These 'init' values are used for Instrument Parameters 
    where there is both a cuurent value for the parameter as 
    well as an initial value (usually read from the instrument 
    parameters file whne the histogram is read. If 'init' is 
    present in `innerdict`, there will also be a 'val' entry 
    in `innerdict` and likely a 'ref' entry as well.

 *  The 'rowlbl' value provides a reference to a str value that 
    will be an editable row label (FreePrmX sample parametric 
    values). 

'''

# import math
# import os
# import re
# import copy
# import platform
# import pickle
# import sys
# import random as ran

#import numpy as np
# import numpy.ma as ma
import wx

# from . import GSASIIpath
from . import GSASIIdataGUI as G2gd
# from . import GSASIIobj as G2obj
from . import GSASIIpwdGUI as G2pdG
# from . import GSASIIimgGUI as G2imG
# from . import GSASIIElem as G2el
# from . import GSASIIfiles as G2fil
# from . import GSASIIctrlGUI as G2G
# from . import GSASIImath as G2mth
# from . import GSASIIElem as G2elem
from . import GSASIIspc as G2spc
from . import GSASIIlattice as G2lat
# from . import GSASIIpwd as G2pwd
from . import GSASIIctrlGUI as G2G
from . import GSASIIpwdplot as G2pwpl
WACV = wx.ALIGN_CENTER_VERTICAL

def UpdateGroup(G2frame,item,plot=True):
    def onDisplaySel(event):
        G2frame.GroupInfo['displayMode'] = dsplType.GetValue()
        wx.CallAfter(UpdateGroup,G2frame,item,False)
        #UpdateGroup(G2frame,item)
    def OnCopyAll(event):
        G2G.G2MessageBox(G2frame,
                    f'Sorry, not fully implemented yet',
                                     'In progress')
        return
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        groupDict = Controls.get('Groups',{}).get('groupDict',{})
        groupName = G2frame.GroupInfo['groupName']
        # make a list of groups of the same length as the current
        curLen = len(groupDict[groupName])
        matchGrps = []
        for g in groupDict:
            if g == groupName: continue
            if curLen != len(groupDict[g]): continue
            matchGrps.append(g)
        if len(matchGrps) == 0:
            G2G.G2MessageBox(G2frame,
                    f'No groups found with {curLen} histograms',
                                     'No matching groups')
            return
        elif len(matchGrps) > 1:
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Copy to which groups?', 'Copy to?', matchGrps)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    selList = [matchGrps[i] for i in dlg.GetSelections()]
            finally:
                dlg.Destroy()
        else:
            selList = matchGrps
        if len(selList) == 0: return
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        if G2frame.GroupInfo['displayMode'].startswith('Sample'):
            prmTable = getSampleVals(G2frame,Histograms)
        elif G2frame.GroupInfo['displayMode'].startswith('Instrument'):
            prmTable = getInstVals(G2frame,Histograms)
        elif G2frame.GroupInfo['displayMode'].startswith('Limits'):
            CopyCtrl = False
            prmTable = getLimitVals(G2frame,Histograms)
        elif G2frame.GroupInfo['displayMode'].startswith('Background'):
            prmTable = getBkgVals(G2frame,Histograms)
            CopyCtrl = False
        else:
            print('Unexpected', G2frame.GroupInfo['displayMode'])
            return
        for h in selList: # group
            for src,dst in zip(groupDict[groupName],groupDict[h]):
                print('copy',src,'to',dst)
                for i in prmTable[src]:
                    #if i not in prmTable[dst]:
                    #    print
                    #    continue
                    if 'val' in  prmTable[src][i]:
                        breakpoint()
        # so what do we copy?
        #breakpoint()
        print('OnCopyAll')
    def OnCopySel(event):
        print('OnCopySel')
        G2G.G2MessageBox(G2frame,
                    f'Sorry, not fully implemented yet',
                                     'In progress')
        return

    if not hasattr(G2frame,'GroupInfo'):
        G2frame.GroupInfo = {}
    G2frame.GroupInfo['displayMode'] = G2frame.GroupInfo.get('displayMode','Sample')
    G2frame.GroupInfo['groupName'] = G2frame.GPXtree.GetItemText(item)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.GroupMenu)
    G2frame.Bind(wx.EVT_MENU, OnCopyAll, id=G2G.wxID_GRPALL)
    G2frame.Bind(wx.EVT_MENU, OnCopySel, id=G2G.wxID_GRPSEL)

    G2frame.dataWindow.ClearData()
    G2frame.dataWindow.helpKey = "Groups/Powder"
    topSizer = G2frame.dataWindow.topBox
    topParent = G2frame.dataWindow.topPanel
    dsplType = wx.ComboBox(topParent,wx.ID_ANY,
                               value=G2frame.GroupInfo['displayMode'],
                               choices=['Hist/Phase','Sample','Instrument',
                                            'Instrument-\u0394',
                                            'Limits','Background'],
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
    dsplType.Bind(wx.EVT_COMBOBOX, onDisplaySel)
    topSizer.Add(dsplType,0,WACV)
    topSizer.Add(wx.StaticText(topParent,
            label=f' parameters for group "{histLabels(G2frame)[0]}"'),
                     0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(topParent,helpIndex=G2frame.dataWindow.helpKey))

    if G2frame.GroupInfo['displayMode'].startswith('Hist'):
        HAPframe(G2frame)
    else:
        HistFrame(G2frame)
    if plot: G2pwpl.PlotPatterns(G2frame,plotType='GROUP')
    G2frame.dataWindow.SetDataSize()
    #wx.CallLater(100,G2frame.SendSizeEvent)
    wx.CallAfter(G2frame.SendSizeEvent)

def histLabels(G2frame):
    '''Find portion of the set of hist names that are the same for all
    histograms in the current group (determined by ``G2frame.GroupInfo['groupName']``)
    and then for each histogram, the characters that are different.

    :Returns: commonltrs, histlbls where 

      * commonltrs is a str containing the letters shared by all 
        histograms in the group and where differing letters are 
        replaced by a square box.
      * histlbls is a list with an str for each histogram listing
        the characters that differ in each histogram.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(
        G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.GroupInfo['groupName']
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    l = max([len(i) for i in groupDict[groupName]])
    h0 = groupDict[groupName][0].ljust(l)
    msk = [True] * l
    for h in groupDict[groupName][1:]:
        msk = [m & (h0i == hi) for h0i,hi,m in zip(h0,h.ljust(l),msk)]
    # place rectangular box in the loc of non-common letter(s)
    commonltrs = ''.join([h0i if m else '\u25A1' for (h0i,m) in zip(h0,msk)])
    # make list with histogram name unique letters   
    histlbls = [''.join([hi for (hi,m) in zip(h,msk) if not m])
                   for h in groupDict[groupName]]
    return commonltrs,histlbls

def onRefineAll(event):
    '''Respond to the Refine All button. On the first press, all 
    refine check buttons are set as "on" and the button is relabeled 
    as 'C' (for clear). On the second press,  all refine check 
    buttons are set as "off" and the button is relabeled as 'S' (for 
    set).
    '''
    but = event.GetEventObject()
    refList = but.refList
    checkButList = but.checkButList
    
    #print('before',[item[0][item[1]] for item in refList])
    if but.GetLabelText() == 'S':
        setting = True
        but.SetLabelText('C')
    else:
        setting = False
        but.SetLabelText('S')
    for c in checkButList:
        c.SetValue(setting)
    for item in refList:
        item[0][item[1]] = setting
    #print('after ',[item[0][item[1]] for item in refList])

def onSetAll(event):
    '''Respond to the copy right button. Copies the first value to 
    all edit widgets
    '''
    but = event.GetEventObject()
    valList = but.valList
    valEditList = but.valEditList
    #print('before',[item[0][item[1]] for item in valList])
    firstVal = valList[0][0][valList[0][1]]
    for c in valEditList:
        c.ChangeValue(firstVal)
    #print('after',[item[0][item[1]] for item in valList])

def displayDataTable(rowLabels,Table,Sizer,Panel,lblRow=False,deltaMode=False,
                     lblSizer=None,lblPanel=None,CopyCtrl=True):
    '''Displays the data table in `Table` in Scrolledpanel `Panel`
    with wx.FlexGridSizer `Sizer`.
    '''
    firstentry = None
    #lblRow = True
    if lblSizer is None: lblSizer = Sizer
    if lblPanel is None: lblPanel = Panel
    checkButList = {}
    valEditList = {}
    lblDict = {}
    for row in rowLabels:
        checkButList[row] = []
        valEditList[row] = []
        # show the row labels, when not in a separate sizer
        if lblRow:
            # is a copy across and/or a refine all button needed?
            refList = []
            valList = []
            for hist in Table:
                if row not in Table[hist]: continue
                if 'val' in Table[hist][row]:
                    valList.append(Table[hist][row]['val'])
                if 'ref' in Table[hist][row]:
                    refList.append(Table[hist][row]['ref'])

            arr = None
            for hist in Table:
                if row not in Table[hist]: continue
                if 'rowlbl' in Table[hist][row]:
                    arr,key = Table[hist][row]['rowlbl']
                    break
            if arr is None: # text row labels
                w = wx.StaticText(lblPanel,label=row)
            else: # used for "renameable" sample vars (str)
                w = G2G.ValidatedTxtCtrl(lblPanel,arr,key,size=(125,-1))
            lblSizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
            lblDict[row] = w

            if len(refList) > 2:
                lbl = 'S'
                if all([l[i] for l,i in refList]): lbl = 'C'
                refAll = wx.Button(lblPanel,label=lbl,style=wx.BU_EXACTFIT)
                refAll.refList = refList
                refAll.checkButList = checkButList[row]
                lblSizer.Add(refAll,0,wx.ALIGN_CENTER_VERTICAL)               
                refAll.Bind(wx.EVT_BUTTON,onRefineAll)
            else:
                lblSizer.Add((-1,-1))
            # if len(valList) > 2:
            #     but = wx.Button(lblPanel,wx.ID_ANY,'\u2192',style=wx.BU_EXACTFIT)
            #     but.valList = valList
            #     but.valEditList = valEditList[row] 
            #     lblSizer.Add(but,0,wx.ALIGN_CENTER_VERTICAL)
            #     but.Bind(wx.EVT_BUTTON,onSetAll)
            # else:
            #     lblSizer.Add((-1,-1))
 
        for i,hist in enumerate(Table):
            if i == 1 and len(valList) > 2 and not deltaMode and CopyCtrl:
                but = wx.Button(Panel,wx.ID_ANY,'\u2192',style=wx.BU_EXACTFIT)
                but.valList = valList
                but.valEditList = valEditList[row]
                Sizer.Add(but,0,wx.ALIGN_CENTER_VERTICAL)
                but.Bind(wx.EVT_BUTTON,onSetAll)
            elif i == 1 and CopyCtrl:
                Sizer.Add((-1,-1))
            minval = None
            maxval = None
            # format the entry depending on what is defined
            if row not in Table[hist]:
                Sizer.Add((-1,-1))
                continue
            elif 'range' in Table[hist][row]:
                minval, maxval = Table[hist][row]['range']
            if ('init' in Table[hist][row] and
                      deltaMode and 'ref' in Table[hist][row]):
                arr,indx = Table[hist][row]['val']
                delta = arr[indx]
                arr,indx = Table[hist][row]['init']
                delta -= arr[indx]
                if abs(delta) < 9e-6: delta = 0.
                if delta == 0:
                    deltaS = ""
                else:
                    deltaS = f"\u0394 {delta:.4g} "
                valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                valrefsiz.Add(wx.StaticText(Panel,label=deltaS),0)
                arr,indx = Table[hist][row]['ref']
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                valrefsiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL)
                checkButList[row].append(w)
                Sizer.Add(valrefsiz,0,
                                 wx.EXPAND|wx.ALIGN_RIGHT)
            elif 'init' in Table[hist][row] and deltaMode:
                # does this ever happen?
                arr,indx = Table[hist][row]['val']
                delta = arr[indx]
                arr,indx = Table[hist][row]['init']
                delta -= arr[indx]
                if delta == 0:
                    deltaS = ""
                else:
                    deltaS = f"\u0394 {delta:.4g} "
                Sizer.Add(wx.StaticText(Panel,label=deltaS),0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
            elif 'val' in Table[hist][row] and 'ref' in Table[hist][row]:
                valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                arr,indx = Table[hist][row]['val']
                w = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                                         nDig=[9,7,'g'],
                                    xmin=minval,xmax=maxval)
                valEditList[row].append(w)
                valrefsiz.Add(w,0,WACV)
                if firstentry is None: firstentry = w
                arr,indx = Table[hist][row]['ref']
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                valrefsiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL)
                checkButList[row].append(w)
                Sizer.Add(valrefsiz,0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
            elif 'val' in Table[hist][row]:
                arr,indx = Table[hist][row]['val']
                nDig = [9,7,'g']
                if type(arr[indx]) is str: nDig = None
                w = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                            nDig=nDig,
                            xmin=minval,xmax=maxval,notBlank=False)
                valEditList[row].append(w)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
                if firstentry is None: firstentry = w

            elif 'ref' in Table[hist][row]:
                arr,indx = Table[hist][row]['ref']
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
                checkButList[row].append(w)
            elif 'str' in Table[hist][row]:
                Sizer.Add(wx.StaticText(Panel,label=Table[hist][row]['str']),0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
            else:
                print('Should not happen',Table[hist][row],hist,row)
    return firstentry,lblDict


def HistFrame(G2frame):
    '''Give up on side-by-side scrolled panels. Put everything 
    in a single FlexGridSizer.
    '''
    #---------------------------------------------------------------------
    # generate a dict with values for each histogram
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    CopyCtrl = True
    if G2frame.GroupInfo['displayMode'].startswith('Sample'):
        prmTable = getSampleVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Instrument'):
        prmTable = getInstVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Limits'):
        CopyCtrl = False
        prmTable = getLimitVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Background'):
        prmTable = getBkgVals(G2frame,Histograms)
        CopyCtrl = False
    else:
        print('Unexpected', G2frame.GroupInfo['displayMode'])
        return
    #debug# for hist in prmTable: printTable(phaseList[page],hist,prmTable[hist]) # see the dict
    # construct a list of row labels, attempting to keep the
    # order they appear in the original array
    rowLabels = []
    lpos = 0
    for hist in prmTable:
        prevkey = None
        for key in prmTable[hist]:
            if key not in rowLabels:
                if prevkey is None:
                    rowLabels.insert(lpos,key)
                    lpos += 1
                else:
                    rowLabels.insert(rowLabels.index(prevkey)+1,key)
            prevkey = key
    #=======  Generate GUI ===============================================
    # layout the window
    panel = midPanel = G2frame.dataWindow
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    G2G.HorizontalLine(mainSizer,panel)
    panel.SetSizer(mainSizer)
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    deltaMode = "\u0394" in G2frame.GroupInfo['displayMode']
    n = 2
    if CopyCtrl: n += 1 # add column for copy
    valSizer = wx.FlexGridSizer(0,len(prmTable)+n,3,10)
    mainSizer.Add(valSizer,1,wx.EXPAND)
    valSizer.Add(wx.StaticText(midPanel,label=' '))
    valSizer.Add(wx.StaticText(midPanel,label=' Ref '))
    #valSizer.Add(wx.StaticText(midPanel,label=' Copy '))
    for i,hist in enumerate(histLabels(G2frame)[1]):
            if i == 1 and CopyCtrl:
                if deltaMode:
                    valSizer.Add((-1,-1))
                elif CopyCtrl:
                    valSizer.Add(wx.StaticText(midPanel,label=' Copy '))
            valSizer.Add(wx.StaticText(midPanel,
                        label=f"\u25A1 = {hist}"),
                             0,wx.ALIGN_CENTER)
    firstentry,lblDict = displayDataTable(rowLabels,prmTable,valSizer,midPanel,
                                      lblRow=True,
                                      deltaMode=deltaMode,CopyCtrl=CopyCtrl)
    #G2frame.dataWindow.SetDataSize()
    if firstentry is not None:    # prevent scroll to show last entry
        wx.Window.SetFocus(firstentry)
        firstentry.SetInsertionPoint(0) # prevent selection of text in widget

def getSampleVals(G2frame,Histograms):
    '''Generate the Parameter Data Table (a dict of dicts) with 
    all Sample values for all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for Sample values.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.GroupInfo['groupName']
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # parameters to include in table
    parms = []
    parmDict = {}
    #indexDict = {}
   # loop over histograms in group
    for hist in groupDict[groupName]:
        #indexDict[hist] = {}
        histdata = Histograms[hist]
        hpD = {}
        hpD['Inst. name'] = {
            'val' : (histdata['Sample Parameters'],'InstrName')}
        #indexDict[hist]['Inst. name'] = {
        #    'val' : (hist,'Sample Parameters','InstrName')}
        hpD['Diff type'] = {
            'str' : histdata['Sample Parameters']['Type']}
        #indexDict[hist]['Diff type'] = {
        #    'str' : (hist,'Sample Parameters','Type')}
        arr = histdata['Sample Parameters']['Scale']
        hpD['Scale factor'] = {
            'val' : (arr,0),
            'ref' : (arr,1),}
        #indexDict[hist]['Scale factor'] = {
        #    'val' : (hist,'Sample Parameters','Scale',0),
        #    'ref' : (hist,1),}
        #breakpoint()
        #return
        # make a list of parameters to show
        histType = histdata['Instrument Parameters'][0]['Type'][0]
        dataType = histdata['Sample Parameters']['Type']
        if histType[2] in ['A','B','C']:
            parms.append(['Gonio. radius','Gonio radius','.3f'])
        #if 'PWDR' in histName:
        if dataType == 'Debye-Scherrer':
            if 'T' in histType:
                parms += [['Absorption','Sample abs, \xb5r/\u03bb',None,]]
            else:
                parms += [['DisplaceX',u'Sample X displ',None,],
                    ['DisplaceY','Sample Y displ',None,],
                    ['Absorption','Sample abs,\xb5\xb7r',None,]]
        elif dataType == 'Bragg-Brentano':
            parms += [['Shift','Sample displ',None,],
                ['Transparency','Sample transp',None],
                ['SurfRoughA','Surf rough A',None],
                ['SurfRoughB','Surf rough B',None]]
        #elif 'SASD' in histName:
        #    parms.append(['Thick','Sample thickness (mm)',[10,3]])
        #    parms.append(['Trans','Transmission (meas)',[10,3]])
        #    parms.append(['SlitLen',u'Slit length (Q,\xc5'+Pwrm1+')',[10,3]])
        parms.append(['Omega','Gonio omega',None])
        parms.append(['Chi','Gonio chi',None])
        parms.append(['Phi','Gonio phi',None])
        parms.append(['Azimuth','Detect azimuth',None])
        parms.append(['Time','time',None])
        parms.append(['Temperature','Sample T',None])
        parms.append(['Pressure','Sample P',None])

        # and loop over them
        for key,lbl,fmt in parms:
            if fmt is None and type(histdata['Sample Parameters'][key]) is list:
                arr = histdata['Sample Parameters'][key]
                hpD[lbl] = {
                    'val' : (arr,0),
                    'ref' : (arr,1),}
            elif fmt is None:
                hpD[lbl] = {
                     'val' : (histdata['Sample Parameters'],key)}
            elif type(fmt) is str:
                 hpD[lbl] = {
                     'str' : f"{histdata['Sample Parameters'][key]:{fmt}}"}

        for key in ('FreePrm1','FreePrm2','FreePrm3'):
            lbl = Controls[key]
            hpD[lbl] = {
                     'val' : (histdata['Sample Parameters'],key),
                     'rowlbl' : (Controls,key)
                }
        parmDict[hist] = hpD
    return parmDict

def getInstVals(G2frame,Histograms):
    '''Generate the Parameter Data Table (a dict of dicts) with 
    all Instrument Parameter values for all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI values.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.GroupInfo['groupName']
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # parameters to include in table
    parms = []
    parmDict = {}
   # loop over histograms in group
    for hist in groupDict[groupName]:
        histdata = Histograms[hist]
        insVal = Histograms[hist]['Instrument Parameters'][0]
        hpD = {}
        insType = insVal['Type'][1]
        try:
            hpD['Bank'] = {
                'str' : str(int(insVal['Bank'][1]))}
        except:
            pass
        hpD['Hist type'] = {
            'str' : insType}
        if insType[2] in ['A','B','C']:               #constant wavelength
            keylist = [('Azimuth','Azimuth','.3f'),]
            if 'Lam1' in insVal:
                keylist += [('Lam1','Lambda 1','.6f'),
                            ('Lam2','Lambda 2','.6f'),
                            (['Source',1],'Source','s'),
                            ('I(L2)/I(L1)','I(L2)/I(L1)',None)]
            else:
                keylist += [('Lam','Lambda',None),]
            itemList = ['Zero','Polariz.']
            if 'C' in insType:
                itemList += ['U','V','W','X','Y','Z','SH/L']
            elif 'B' in insType:
                itemList += ['U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1']
            else: #'A'
                itemList += ['U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1','SH/L']
            for lbl in itemList:
                keylist += [(lbl,lbl,None),]
        elif 'E' in insType:
            for lbl in ['XE','YE','ZE','WE']:
                keylist += [(lbl,lbl,'.6f'),]
            for lbl in ['A','B','C','X','Y','Z']:
                keylist += [(lbl,lbl,None),]
        elif 'T' in insType:
            keylist = [('fltPath','Flight path','.3f'),
                       ('2-theta','2\u03B8','.2f'),]
            for lbl in ['difC','difA','difB','Zero','alpha',
                            'beta-0','beta-1','beta-q',
                            'sig-0','sig-1','sig-2','sig-q','X','Y','Z']:
                keylist += [(lbl,lbl,None),]
        else:
            return {}
        for key,lbl,fmt in keylist:
                arr = insVal[key]
                if fmt is None:
                    hpD[lbl] = {
                        'init' : (arr,0),
                        'val' : (arr,1),
                        'ref' : (arr,2),}
                else:
                    hpD[lbl] = {
                        'str' : f'{arr[1]:{fmt}}'}
        parmDict[hist] = hpD
    return parmDict

def getLimitVals(G2frame,Histograms):
    '''Generate the Limits Data Table (a dict of dicts) with 
    all limits values for all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for limits values.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.GroupInfo['groupName']
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # parameters to include in table
    parmDict = {}
   # loop over histograms in group
    for hist in groupDict[groupName]:
        histdata = Histograms[hist]
        hpD = {}
        #breakpoint()
        hpD['Tmin'] = {
            'val' : (histdata['Limits'][1],0),
            'range': [histdata['Limits'][0][0],histdata['Limits'][0][1]]
            }
        hpD['Tmax'] = {
            'val' : (histdata['Limits'][1],1),
            'range': [histdata['Limits'][0][0],histdata['Limits'][0][1]]
            }
        for i,item in enumerate(histdata['Limits'][2:]):
            hpD[f'excl Low {i+1}'] = {
                'val' : (item,0),
                'range': [histdata['Limits'][0][0],histdata['Limits'][0][1]]}
            hpD[f'excl High {i+1}'] = {
                'val' : (item,1),
                'range': [histdata['Limits'][0][0],histdata['Limits'][0][1]]}
        parmDict[hist] = hpD
    # for i in parmDict:
    #     print(i)
    #     for j in  parmDict[i]:
    #         print('\t',j)
    #         for k in parmDict[i][j]:
    #             print('\t\t',k,parmDict[i][j][k])
    return parmDict


def getBkgVals(G2frame,Histograms):
    '''Generate the Background Data Table (a dict of dicts) with 
    all Background values for all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for 
    Background values.
    '''
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.GroupInfo['groupName']
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # parameters to include in table
    parms = []
    parmDict = {}
    # loop over histograms in group
    for hist in groupDict[groupName]:
        histdata = Histograms[hist]
        hpD = {}
        hpD['Function'] = {
            'str' : histdata['Background'][0][0]}
        hpD['ref flag'] = {
            'ref' : (histdata['Background'][0],1)}
        hpD['# Bkg terms'] = {
            'str' : str(int(histdata['Background'][0][2]))}
        hpD['# Debye terms'] = {
            'str' : str(int(histdata['Background'][1]['nDebye']))}
        for i,term in enumerate(histdata['Background'][1]['debyeTerms']):
            hpD[f'A #{i+1}'] = {
                'val' : (histdata['Background'][1]['debyeTerms'][i],0),
                'ref' : (histdata['Background'][1]['debyeTerms'][i],1),
                }
            hpD[f'R #{i+1}'] = {
                'val' : (histdata['Background'][1]['debyeTerms'][i],2),
                'ref' : (histdata['Background'][1]['debyeTerms'][i],3),
                }
            hpD[f'U #{i+1}'] = {
                'val' : (histdata['Background'][1]['debyeTerms'][i],4),
                'ref' : (histdata['Background'][1]['debyeTerms'][i],5),
                }
        hpD['# Bkg Peaks'] = {
            'str' : str(int(histdata['Background'][1]['nPeaks']))}
        for i,term in enumerate(histdata['Background'][1]['peaksList']):
            hpD[f'pos #{i+1}'] = {
                'val' : (histdata['Background'][1]['peaksList'][i],0),
                'ref' : (histdata['Background'][1]['peaksList'][i],1),
                }
            hpD[f'int #{i+1}'] = {
                'val' : (histdata['Background'][1]['peaksList'][i],2),
                'ref' : (histdata['Background'][1]['peaksList'][i],3),
                }
            hpD[f'sig #{i+1}'] = {
                'val' : (histdata['Background'][1]['peaksList'][i],4),
                'ref' : (histdata['Background'][1]['peaksList'][i],5),
                }
            hpD[f'gam #{i+1}'] = {
                'val' : (histdata['Background'][1]['peaksList'][i],6),
                'ref' : (histdata['Background'][1]['peaksList'][i],7),
                }
        if histdata['Background'][1]['background PWDR'][0]:
            val = 'yes'
        else:
            val = 'no'
        hpD['Fixed bkg file'] = {
            'str' : val}
        parmDict[hist] = hpD
    return parmDict


def HAPframe(G2frame):
    '''This creates two side-by-side scrolled panels, each containing 
    a FlexGridSizer.
    The panel to the left contains the labels for the sizer to the right.
    This way the labels are not scrolled horizontally and are always seen.
    The two vertical scroll bars are linked together so that the labels 
    are synced to the table of values.
    '''
    def selectPhase(event):
        'Display the selected phase'
        def OnScroll(event):
            'Synchronize vertical scrolling between the two scrolled windows'
            obj = event.GetEventObject()
            pos = obj.GetViewStart()[1]
            if obj == lblScroll:
                HAPScroll.Scroll(-1, pos)
            else:
                lblScroll.Scroll(-1, pos)
            event.Skip()
        #---------------------------------------------------------------------
        # selectPhase starts here. Find which phase is selected.
        if event:
            page = event.GetSelection()
            #print('page selected',page,phaseList[page])
        else:  # initial call when window is created
            page = 0
            #print('no page selected',phaseList[page])
        # generate a dict with HAP values for each phase (may not be the same)
        HAPtable = getHAPvals(G2frame,phaseList[page])
        #debug# for hist in HAPtable: printTable(phaseList[page],hist,HAPtable[hist]) # see the dict
        # construct a list of row labels, attempting to keep the
        # order they appear in the original array
        rowLabels = []
        lpos = 0
        for hist in HAPtable:
            prevkey = None
            for key in HAPtable[hist]:
                if key not in rowLabels:
                    if prevkey is None:
                        rowLabels.insert(lpos,key)
                        lpos += 1
                    else:
                        rowLabels.insert(rowLabels.index(prevkey)+1,key)
                prevkey = key
        #=======  Generate GUI ===============================================
        for panel in HAPtabs:
            if panel.GetSizer():
                panel.GetSizer().Destroy()          # clear out old widgets
        panel = HAPtabs[page]
        bigSizer = wx.BoxSizer(wx.HORIZONTAL)
        panel.SetSizer(bigSizer)
        # panel for labels; show scroll bars to hold the space
        lblScroll = wx.lib.scrolledpanel.ScrolledPanel(panel,
                        style=wx.VSCROLL|wx.HSCROLL|wx.ALWAYS_SHOW_SB)
        hpad = 3  # space between rows
        #lblSizer = wx.FlexGridSizer(0,3,hpad,2)
        lblSizer = wx.FlexGridSizer(0,2,hpad,2)
        lblScroll.SetSizer(lblSizer)
        bigSizer.Add(lblScroll,0,wx.EXPAND)

        # Create scrolled panel to display HAP data
        HAPScroll = wx.lib.scrolledpanel.ScrolledPanel(panel,
                        style=wx.VSCROLL|wx.HSCROLL|wx.ALWAYS_SHOW_SB)
        #HAPSizer = wx.FlexGridSizer(0,len(HAPtable),hpad,10)
        HAPSizer = wx.FlexGridSizer(0,len(HAPtable)+1,hpad,10)
        HAPScroll.SetSizer(HAPSizer)
        bigSizer.Add(HAPScroll,1,wx.EXPAND)
        
        # Bind scroll events to synchronize scrolling
        lblScroll.Bind(wx.EVT_SCROLLWIN, OnScroll)
        HAPScroll.Bind(wx.EVT_SCROLLWIN, OnScroll)
        # label columns with unique part of histogram names
        for i,hist in enumerate(histLabels(G2frame)[1]):
            if i == 1:
                HAPSizer.Add(wx.StaticText(HAPScroll,label='Copy'),
                             0,wx.ALIGN_CENTER)
            HAPSizer.Add(wx.StaticText(HAPScroll,label=f"\u25A1 = {hist}"),
                             0,wx.ALIGN_CENTER)
        w0 = wx.StaticText(lblScroll,label=' ')
        lblSizer.Add(w0)
        lblSizer.Add(wx.StaticText(lblScroll,label=' Ref '))
        #lblSizer.Add(wx.StaticText(lblScroll,label=' Copy '))
        firstentry,lblDict = displayDataTable(rowLabels,HAPtable,HAPSizer,HAPScroll,
                    lblRow=True,lblSizer=lblSizer,lblPanel=lblScroll)
        # get row sizes in data table
        HAPSizer.Layout()
        rowHeights = HAPSizer.GetRowHeights()
        # set row sizes in Labels
        # (must be done after HAPSizer row heights are defined)
        s = wx.Size(-1,rowHeights[0])
        w0.SetMinSize(s)
        # for i,row in enumerate(rowLabels):
        #     s = wx.Size(-1,rowHeights[i+1])
        #     w = wx.StaticText(lblScroll,label=row,size=s)
        #     lblSizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
        # lblDict = {}
        # for i,row in enumerate(rowLabels):
        #     w = wx.StaticText(lblScroll,label=row)
        #     lblDict[row] = w
        #     lblSizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
        for i,row in enumerate(rowLabels):
            s = wx.Size(-1,rowHeights[i+1])
            lblDict[row].SetMinSize(s)
        # Fit the scrolled windows to their content
        lblSizer.Layout()
        xLbl,_ = lblSizer.GetMinSize()
        xTab,yTab = HAPSizer.GetMinSize()
        lblScroll.SetSize((xLbl,yTab))
        lblScroll.SetMinSize((xLbl+15,yTab)) # add room for scroll bar
        lblScroll.SetVirtualSize(lblSizer.GetMinSize())
        HAPScroll.SetVirtualSize(HAPSizer.GetMinSize())
        lblScroll.SetupScrolling(scroll_x=True, scroll_y=True, rate_x=20, rate_y=20)
        HAPScroll.SetupScrolling(scroll_x=True, scroll_y=True, rate_x=20, rate_y=20)
        if firstentry is not None:    # prevent scroll to show last entry
            wx.Window.SetFocus(firstentry)
            firstentry.SetInsertionPoint(0) # prevent selection of text in widget

    #G2frame.dataWindow.ClearData()

    # layout the HAP window. This has histogram and phase info, so a
    # notebook is needed for phase name selection. (That could
    # be omitted for single-phase refinements, but better to remind the
    # user of the phase
    # topSizer = G2frame.dataWindow.topBox
    # topParent = G2frame.dataWindow.topPanel
    midPanel = G2frame.dataWindow
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    #botSizer = G2frame.dataWindow.bottomBox
    #botParent = G2frame.dataWindow.bottomPanel
    # label with shared portion of histogram name
    # topSizer.Add(wx.StaticText(topParent,
    #         label=f'HAP parameters for group "{histLabels(G2frame)[0]}"'),
    #                  0,WACV)
    # topSizer.Add((-1,-1),1,wx.EXPAND)
    # topSizer.Add(G2G.HelpButton(topParent,helpIndex=G2frame.dataWindow.helpKey))

    G2G.HorizontalLine(mainSizer,midPanel)
    midPanel.SetSizer(mainSizer)
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    if not Phases:
        mainSizer.Add(wx.StaticText(midPanel,
                           label='There are no phases in use'))
        G2frame.dataWindow.SetDataSize()
        return
    # notebook for phases
    HAPBook = G2G.GSNoteBook(parent=midPanel)
    mainSizer.Add(HAPBook,1,wx.ALL|wx.EXPAND,1)
    HAPtabs = []
    phaseList = []
    for phaseName in Phases:
        phaseList.append(phaseName)
        HAPtabs.append(wx.Panel(HAPBook))
        HAPBook.AddPage(HAPtabs[-1],phaseName)
    HAPBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, selectPhase)

    page = 0
    HAPBook.SetSelection(page)
    selectPhase(None)
    #G2frame.dataWindow.SetDataSize()

def getHAPvals(G2frame,phase):
    '''Generate the Parameter Data Table (a dict of dicts) with 
    all HAP values for the selected phase and all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for HAP values.
    '''
    sub = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
    item, cookie = G2frame.GPXtree.GetFirstChild(sub)
    PhaseData = None
    while item: # loop over phases
        phaseName = G2frame.GPXtree.GetItemText(item)
        if phase is None: phase = phaseName
        if phase == phaseName:
            PhaseData = G2frame.GPXtree.GetItemPyData(item)
            break
        item, cookie = G2frame.GPXtree.GetNextChild(sub, cookie)
    if PhaseData is None:
        print(f'Unexpected: Phase {phase!r} not found')
        return {}

    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    groupName = G2frame.GroupInfo['groupName']
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # loop over histograms in group
    parmDict = {}
    for hist in groupDict[groupName]:
        parmDict[hist] = makeHAPtbl(G2frame,phase,PhaseData,hist)
        #printTable(phase,hist,parmDict[hist])
    return parmDict

def makeHAPtbl(G2frame,phase,PhaseData,hist):
    '''Construct a Parameter Data Table dict, providing access 
    to the HAP variables for one phase and histogram.

    :return: the Parameter Data Table dict, as described above.
    '''
    SGData = PhaseData['General']['SGData']
    cell = PhaseData['General']['Cell'][1:]
    Amat,Bmat = G2lat.cell2AB(cell[:6])
    #G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, hist)
    #data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
    HAPdict = PhaseData['Histograms'][hist]

    parmDict = {}
    # phase fraction
    parmDict['Phase frac'] = {
        'val' : (HAPdict['Scale'],0),
        'ref' : (HAPdict['Scale'],1),}

    parmDict['LeBail extract'] = {
            'str' : "Yes" if HAPdict.get('LeBail') else '(off)'
        }

    # size values
    arr = HAPdict['Size']
    if arr[0] == 'isotropic':
        parmDict['Size'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
    elif arr[0] == 'uniaxial':
        parmDict['Size/Eq'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
        parmDict['Size/Ax'] = {
            'val' : (arr[1],1),
            'ref' : (arr[2],1),}
        parmDict['Size/dir'] = {
            'str' : ','.join([str(i) for i in arr[3]])}
    else:
        for i,lbl in enumerate(['S11','S22','S33','S12','S13','S23']):
            parmDict[f'Size/{lbl}'] = {
            'val' : (arr[4],i),
            'ref' : (arr[5],i),}
    parmDict['Size LGmix'] = {
        'val' : (arr[1],2),
        'ref' : (arr[2],2),}

    # microstrain values
    arr = HAPdict['Mustrain']
    if arr[0] == 'isotropic':
        parmDict['\u00B5Strain'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
    elif arr[0] == 'uniaxial':
        parmDict['\u00B5Strain/Eq'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
        parmDict['\u00B5Strain/Ax'] = {
            'val' : (arr[1],1),
            'ref' : (arr[2],1),}
        parmDict['\u00B5Strain/dir'] = {
            'str' : ','.join([str(i) for i in arr[3]])}
    else:
        Snames = G2spc.MustrainNames(SGData)
        for i,lbl in enumerate(Snames):
            if i >= len(arr[4]): break
            parmDict[f'\u00B5Strain/{lbl}'] = {
                'val' : (arr[4],i),
                'ref' : (arr[5],i),}
        muMean = G2spc.MuShklMean(SGData,Amat,arr[4][:len(Snames)])
        parmDict['\u00B5Strain/mean'] = {
            'str' : f'{muMean:.2f}'}
    parmDict['\u00B5Strain LGmix'] = {
        'val' : (arr[1],2),
        'ref' : (arr[2],2),}

    # Hydrostatic terms
    Hsnames = G2spc.HStrainNames(SGData)
    arr = HAPdict['HStrain']
    for i,lbl in enumerate(Hsnames):
        if i >= len(arr[0]): break
        parmDict[f'Size/{lbl}'] = {
            'val' : (arr[0],i),
            'ref' : (arr[1],i),}

    # Preferred orientation terms
    arr = HAPdict['Pref.Ori.']
    if arr[0] == 'MD':
        parmDict['March-Dollase'] = {
            'val' : (arr,1),
            'ref' : (arr,2),}
        parmDict['M-D/dir'] = {
            'str' : ','.join([str(i) for i in arr[3]])}
    else:
        parmDict['Spherical harmonics'] = {
            'ref' : (arr,2),}
        parmDict['SH order'] = {
            'str' : str(arr[4])}
        for lbl in arr[5]:
            parmDict[f'SP {lbl}']= {
                'val' : (arr[5],lbl),
                }
        parmDict['SH text indx'] = {
            'str' : f'{G2lat.textureIndex(arr[5]):.3f}'}

    # misc: Layer Disp, Extinction
    try:
        parmDict['Layer displ'] = {
        'val' : (HAPdict['Layer Disp'],0),
        'ref' : (HAPdict['Layer Disp'],1),}
    except KeyError:
        pass
    try:
        parmDict['Extinction'] = {
            'val' : (HAPdict['Extinction'],0),
            'ref' : (HAPdict['Extinction'],1),}
    except KeyError:
        pass
    try:
        parmDict['Babinet A'] = {
            'val' : (HAPdict['Babinet']['BabA'],0),
            'ref' : (HAPdict['Babinet']['BabA'],1),}
    except KeyError:
        pass
    try:
        parmDict['Babinet U'] = {
            'val' : (HAPdict['Babinet']['BabU'],0),
            'ref' : (HAPdict['Babinet']['BabU'],1),}
    except KeyError:
        pass
    return parmDict

def printTable(phase,hist,parmDict):
    # show the variables and values in data table array -- debug use only
    print(phase,hist)
    for sel in parmDict:
        arr = parmDict[sel]
        val = 'N/A'
        if 'val' in arr:
            val = arr['val'][0] [arr['val'][1]]
        if 'str' in arr:
            v = arr['str']
            val = f'"{v}"'
        ref = 'N/A'
        if 'ref' in arr:
            ref = arr['ref'][0] [arr['ref'][1]]
        print(f'{sel!r:20s}  {val}  {ref}')
    print('\n')

# code that did not work and is being abandoned for now
#
# def HistFrame(G2frame):
#     '''This creates two side-by-side scrolled panels, each containing 
#     a FlexGridSizer.
#     The panel to the left contains the labels for the sizer to the right.
#     This way the labels are not scrolled horizontally and are always seen.
#     The two vertical scroll bars are linked together so that the labels 
#     are synced to the table of values.
#     '''
#     def OnScroll(event):
#         'Synchronize vertical scrolling between the two scrolled windows'
#         obj = event.GetEventObject()
#         pos = obj.GetViewStart()[1]
#         if obj == lblScroll:
#             SampleScroll.Scroll(-1, pos)
#         else:
#             lblScroll.Scroll(-1, pos)
#         event.Skip()
#     #---------------------------------------------------------------------
#     # generate a dict with Sample values for each histogram
#     Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
#     HAPtable = getSampleVals(G2frame,Histograms)
#     #debug# for hist in HAPtable: printTable(phaseList[page],hist,HAPtable[hist]) # see the dict
#     # construct a list of row labels, attempting to keep the
#     # order they appear in the original array
#     rowLabels = []
#     lpos = 0
#     for hist in HAPtable:
#         prevkey = None
#         for key in HAPtable[hist]:
#             if key not in rowLabels:
#                 if prevkey is None:
#                     rowLabels.insert(lpos,key)
#                     lpos += 1
#                 else:
#                     rowLabels.insert(rowLabels.index(prevkey)+1,key)
#             prevkey = key
#     #=======  Generate GUI ===============================================
#     #G2frame.dataWindow.ClearData()

#     # layout the HAP window. This has histogram and phase info, so a
#     # notebook is needed for phase name selection. (That could
#     # be omitted for single-phase refinements, but better to remind the
#     # user of the phase
#     topSizer = G2frame.dataWindow.topBox
#     topParent = G2frame.dataWindow.topPanel
#     midPanel = G2frame.dataWindow
#     # this causes a crash, but somehow I need to clean up the old contents
#     #if midPanel.GetSizer():
#     #    for i in midPanel.GetSizer().GetChildren():
#     #        i.Destroy()          # clear out old widgets
#     #botSizer = G2frame.dataWindow.bottomBox
#     #botParent = G2frame.dataWindow.bottomPanel
#     # label with shared portion of histogram name
#     topSizer.Add(wx.StaticText(topParent,
#             label=f'Sample parameters for group "{histLabels(G2frame)[0]}"'),
#                      0,WACV)
#     topSizer.Add((-1,-1),1,wx.EXPAND)
#     topSizer.Add(G2G.HelpButton(topParent,helpIndex=G2frame.dataWindow.helpKey))

#     panel = wx.Panel(midPanel)
#     panel.SetSize(midPanel.GetSize())
#     mainSizer = wx.BoxSizer(wx.VERTICAL)
#     G2G.HorizontalLine(mainSizer,panel)
#     panel.SetSizer(mainSizer)
#     Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()

    
#     bigSizer = wx.BoxSizer(wx.HORIZONTAL)
#     mainSizer.Add(bigSizer,1,wx.EXPAND)
#     if True:
#         # panel for labels; show scroll bars to hold the space
#         lblScroll = wx.lib.scrolledpanel.ScrolledPanel(panel,
#                         style=wx.VSCROLL|wx.HSCROLL|wx.ALWAYS_SHOW_SB)
#         hpad = 3  # space between rows
#         lblSizer = wx.FlexGridSizer(0,1,hpad,10)
#         lblScroll.SetSizer(lblSizer)
#         bigSizer.Add(lblScroll,0,wx.EXPAND)

#         # Create scrolled panel to display HAP data
#         SampleScroll = wx.lib.scrolledpanel.ScrolledPanel(panel,
#                         style=wx.VSCROLL|wx.HSCROLL|wx.ALWAYS_SHOW_SB)
#         SampleSizer = wx.FlexGridSizer(0,len(HAPtable),hpad,10)
#         SampleScroll.SetSizer(SampleSizer)
#         bigSizer.Add(SampleScroll,1,wx.EXPAND)
        
#         # Bind scroll events to synchronize scrolling
#         lblScroll.Bind(wx.EVT_SCROLLWIN, OnScroll)
#         SampleScroll.Bind(wx.EVT_SCROLLWIN, OnScroll)
#         # label columns with unique part of histogram names
#         for hist in histLabels(G2frame)[1]:
#             SampleSizer.Add(wx.StaticText(SampleScroll,label=f"\u25A1 = {hist}"),
#                              0,wx.ALIGN_CENTER)
    
#         displayDataTable(rowLabels,HAPtable,SampleSizer,SampleScroll)
#         # get row sizes in data table
#         SampleSizer.Layout()
#         rowHeights = SampleSizer.GetRowHeights()
#         # match rose sizes in Labels
#         # (must be done after SampleSizer row heights are defined)
#         s = wx.Size(-1,rowHeights[0])
#         lblSizer.Add(wx.StaticText(lblScroll,label=' ',size=s))
#         for i,row in enumerate(rowLabels):
#             s = wx.Size(-1,rowHeights[i+1])
#             lblSizer.Add(wx.StaticText(lblScroll,label=row,size=s),0,
#                              wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
#         # Fit the scrolled windows to their content
#         lblSizer.Layout()
#         xLbl,_ = lblSizer.GetMinSize()
#         xTab,yTab = SampleSizer.GetMinSize()
#         lblScroll.SetSize((xLbl,yTab))
#         lblScroll.SetMinSize((xLbl+15,yTab)) # add room for scroll bar
#         lblScroll.SetVirtualSize(lblSizer.GetMinSize())
#         SampleScroll.SetVirtualSize(SampleSizer.GetMinSize())
#         lblScroll.SetupScrolling(scroll_x=True, scroll_y=True, rate_x=20, rate_y=20)
#         SampleScroll.SetupScrolling(scroll_x=True, scroll_y=True, rate_x=20, rate_y=20)
#         breakpoint()
#         wx.CallLater(100,G2frame.SendSizeEvent)
#     G2frame.dataWindow.SetDataSize()
