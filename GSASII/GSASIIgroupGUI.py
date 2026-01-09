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

** Parameter Data Table **

For use to create GUI tables and to copy values between histograms, parameters 
are organized in a dict where each dict entry has contents of form:

 *  dict['__dataSource'] : SourceArray

 *  dict[histogram]['label'] : `innerdict`

where `label` is the text shown on the row label and `innerdict` can contain
one or more of the following elements:

 *  'val'    : (key1,key2,...)
 *  'ref'    : (key1, key2,...)
 *  'range'  : (float,float)
 *  'str'    : (key1,key2,...)
 *  'fmt'    : str
 *  'txt'    : str
 *  'init'   : float
 *  'rowlbl' : (array, key)

One of 'val', 'ref' or 'str' elements will be present. 

 *  The 'val' tuple provides a reference to the float value for the 
    defined quantity, where SourceArray[histogram][key1][key2][...] 
    provides r/w access to the parameter.

 *  The 'ref' tuple provides a reference to the bool value, where 
    SourceArray[histogram][key1][key2][...] provides r/w access to the 
    refine flag for the labeled quantity

    Both 'ref' and 'val' are usually defined together, but either may 
    occur alone. These exceptions will be for parameters where a single
    refine flag is used for a group of parameters or for non-refined 
    parameters.

 *  The 'str' value is something that cannot be edited from the GUI; If 'str' is
    present, the only other possible entries that can be present is either 'fmt' 
    or 'txt. 
    'str' is used for a parameter value that is typically computed or must be 
    edited in the histogram section.

 *  The 'fmt' value is a string used to format the 'str' value to 
    convert it to a string, if it is a float or int value.

 *  The 'txt' value is a string that replaces the value in 'str'.

 *  The 'init' value is also something that cannot be edited. 
    These 'init' values are used for Instrument Parameters 
    where there is both a current value for the parameter as 
    well as an initial value (usually read from the instrument 
    parameters file when the histogram is read. If 'init' is 
    present in `innerdict`, there will also be a 'val' entry 
    in `innerdict` and likely a 'ref' entry as well.

 *  The 'rowlbl' value provides a reference to a str value that 
    will be an editable row label (FreePrmX sample parametric 
    values). 

 *  The 'range' list/tuple provides min and max float value for the 
    defined quantity to be defined. Use None for any value that
    should not be enforced. The 'range' values will be used as limits
    for the entry widget.

'''

# import math
# import os
import re
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
# from . import GSASIIpwdGUI as G2pdG
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

def SearchGroups(G2frame,Histograms,hist):
    '''Determine which histograms are in groups, called by SearchGroups in 
    :func:`GSASIIdataGUI.UpdateControls`.
    '''
    repeat = True
    srchStr = hist[5:]
    msg = ('Edit the histogram name below placing a question mark (?) '
            'at the location '
            'of characters that change between groups of '
            'histograms. Use backspace or delete to remove '
            'characters that should be ignored as they will '
            'vary within a histogram group (e.g. Bank 1). '
            'Be sure to leave enough characters so the string '
            'can be uniquely matched.')
    while repeat:
        srchStr = G2G.StringSearchTemplate(G2frame,'Set match template',
                                               msg,srchStr)
        if srchStr is None: return {} # cancel pressed
        reSrch = re.compile(srchStr.replace('.',r'\.').replace('?','.'))
        setDict = {}
        keyList = []
        noMatchCount = 0
        for hist in Histograms:
            if hist.startswith('PWDR '):
                m = reSrch.search(hist)
                if m:
                    key = hist[m.start():m.end()]
                    setDict[hist] = key
                    if key not in keyList: keyList.append(key)
                else:
                    noMatchCount += 1
        groupDict = {}
        groupCount = {}
        for k in keyList:
            groupDict[k] = [hist for hist,key in setDict.items() if k == key]
            groupCount[k] = len(groupDict[k])

        msg1 = f'With search template "{srchStr}", '

        buttons = []
        OK = True
        if len(groupCount) == 0:
            msg1 += f'there are ho histograms in any groups.'
        elif min(groupCount.values()) == max(groupCount.values()):
            msg1 += f'there are {len(groupDict)} groups with {min(groupCount.values())} histograms each.'
        else:
            msg1 += (f'there are {len(groupDict)} groups with between {min(groupCount.values())}'
                       f' and {min(groupCount.values())} histograms in each.')
        if noMatchCount:
            msg1 += f"\n\nNote that {noMatchCount} PWDR histograms were not included in any groups."
        # place a sanity check limit on the number of histograms in a group
        if len(groupCount) == 0:
            OK = False
        elif max(groupCount.values()) >= 150:
            OK = False
            msg1 += '\n\nThis exceeds the maximum group length of 150 histograms'
        elif min(groupCount.values()) == max(groupCount.values()) == 1:
            OK = False
            msg1 += '\n\nEach histogram is in a separate group. Grouping histograms only makes sense with multiple histograms in at least some groups.'
        if OK:
            buttons += [('OK', lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_OK))]
        buttons += [('try again', lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_CANCEL))]
        res = G2G.ShowScrolledInfo(G2frame,msg1,header='Grouping result',
                    buttonlist=buttons,height=150)
        if res == wx.ID_OK:
            repeat = False

    return {'groupDict':groupDict,'notGrouped':noMatchCount,'template':srchStr}

def UpdateGroup(G2frame,item,plot=True):
    def onDisplaySel(event):
        G2frame.GroupInfo['displayMode'] = dsplType.GetValue()
        wx.CallAfter(UpdateGroup,G2frame,item,False)

    def copyPrep():
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        groupDict = Controls.get('Groups',{}).get('groupDict',{})
        groupName = G2frame.GroupInfo['groupName']
        # make a list of groups of the same length as the current
        curLen = len(groupDict[groupName])
        matchGrps = []
        selList = []
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
        return selList,groupDict,groupName

    def OnCopyAll(event):
        res = copyPrep()
        if res is None: return
        selList,groupDict,groupName = res
        dataSource = prmArray['_dataSource']
        for h in selList: # selected groups
            for src,dst in zip(groupDict[groupName],groupDict[h]): # histograms in groups (same length enforced)
                for i in prmArray[src]:
                    for j in ('ref','val','str'):
                        if j in prmArray[src][i]:
                            if prmArray[src][i][j] is None: continue
                            try:
                                arr,indx = indexArrayRef(dataSource,dst,prmArray[src][i][j])
                                arr[indx] = indexArrayVal(dataSource,src,prmArray[src][i][j])
                            except Exception as msg: # could hit an error if an array element is not defined
                                pass
                                #print(msg)
                                #print('error with',i,dst)
    def OnCopySel(event):
        res = copyPrep()
        if res is None: return
        selList,groupDict,groupName = res
        dataSource = prmArray['_dataSource']
        choices = []
        for src in groupDict[groupName]:
            for i in prmArray[src]:
                for j in ('ref','val','str'):
                    if prmArray[src][i].get(j) is None: 
                        continue
                    if i not in choices:
                        choices.append(i)
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Copy which items?', 'Copy what?', choices)
        itemList = []
        try:
            if dlg.ShowModal() == wx.ID_OK:
                itemList = [choices[i] for i in dlg.GetSelections()]
        finally:
            dlg.Destroy()
        if len(itemList) == 0: return
        for h in selList: # selected groups
            for src,dst in zip(groupDict[groupName],groupDict[h]): # histograms in groups (same length enforced)
                for i in prmArray[src]:
                    if i not in itemList: continue
                    for j in ('ref','val','str'):
                        if j in prmArray[src][i]:
                            if prmArray[src][i][j] is None: continue
                            try:
                                arr,indx = indexArrayRef(dataSource,dst,prmArray[src][i][j])
                                arr[indx] = indexArrayVal(dataSource,src,prmArray[src][i][j])
                            except Exception as msg: # could hit an error if an array element is not defined
                                pass
                                #print(msg)
                                #print('error with',i,dst)

    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
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
        HAPframe(G2frame,Histograms,Phases)
    else:
        HistFrame(G2frame,Histograms)
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

def indexArrayRef(dataSource,hist,arrayIndices):
    indx = arrayIndices[-1]
    arr = dataSource[hist]
    for i in  arrayIndices[:-1]:
        arr = arr[i]
    return arr,indx

def indexArrayVal(dataSource,hist,arrayIndices):
    if arrayIndices is None: return None
    arr = dataSource[hist]
    for i in  arrayIndices:
        arr = arr[i]
    return arr

ClearAllLbl = '☐' # previously 'C'
SetAllLbl = '☑'   # previously 'S'
def onRefineAll(event):
    '''Respond to the Refine All button. On the first press, all 
    refine check buttons are set as "on" and the button is relabeled 
    as ClearAllLbl (for clear). On the second press,  all refine check 
    buttons are set as "off" and the button is relabeled as SetAllLbl (for 
    set).
    '''
    but = event.GetEventObject()
    dataSource = but.refDict['dataSource']
    checkButList = but.checkButList
    
    if but.GetLabelText() == SetAllLbl:
        setting = True
        but.SetLabelText(ClearAllLbl)
    else:
        setting = False
        but.SetLabelText(SetAllLbl)
    for c in checkButList:
        c.SetValue(setting)
    for item,hist in zip(but.refDict['arrays'],but.refDict['hists']):
        arr,indx = indexArrayRef(dataSource,hist,item)
        arr[indx] = setting

def onSetAll(event):
    '''Respond to the copy right button. Copies the first value to 
    all edit widgets
    '''
    but = event.GetEventObject()
    dataSource = but.valDict['dataSource']
    valList = but.valDict['arrays']
    histList = but.valDict['hists']
    valEditList = but.valDict['valEditList']
    firstVal = indexArrayVal(dataSource,histList[0],valList[0])
    for c in valEditList:
        c.ChangeValue(firstVal)

def displayDataArray(rowLabels,DataArray,Sizer,Panel,lblRow=False,deltaMode=False,
                     lblSizer=None,lblPanel=None,CopyCtrl=True):
    '''Displays the data table in `DataArray` in Scrolledpanel `Panel`
    with wx.FlexGridSizer `Sizer`.
    '''
    firstentry = None
    #lblRow = True
    if lblSizer is None: lblSizer = Sizer
    if lblPanel is None: lblPanel = Panel
    checkButList = {}
    valEditList = {}
    lblDict = {}
    dataSource = DataArray['_dataSource']
    for row in rowLabels:
        checkButList[row] = []
        valEditList[row] = []
        # show the row labels, when not in a separate sizer
        if lblRow:
            # is a copy across and/or a refine all button needed?
            refList = []
            valList = []
            for hist in DataArray:
                if row not in DataArray[hist]: continue
                if 'val' in DataArray[hist][row]:
                    valList.append(DataArray[hist][row]['val'])
                if 'ref' in DataArray[hist][row]:
                    refList.append(DataArray[hist][row]['ref'])

            arr = None
            histList = []
            for hist in DataArray:
                if row not in DataArray[hist]: continue
                histList.append(hist)
                if 'rowlbl' in DataArray[hist][row]:
                    arr,key = DataArray[hist][row]['rowlbl']
                    break
            if arr is None: # text row labels
                w = wx.StaticText(lblPanel,label=row)
            else: # used for "renameable" sample vars (str)
                w = G2G.ValidatedTxtCtrl(lblPanel,arr,key,size=(125,-1))
            lblSizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
            lblDict[row] = w

            if len(refList) > 2:
                lbl = SetAllLbl
                if all([indexArrayVal(dataSource,hist,i) for i in refList]): lbl = ClearAllLbl
                refAll = wx.Button(lblPanel,label=lbl,style=wx.BU_EXACTFIT)
                font = refAll.GetFont()
                font.PointSize += 5
                refAll.SetFont(font)
                refAll.refDict = {'arrays': refList, 'hists': histList,
                                  'dataSource':dataSource}
                refAll.checkButList = checkButList[row]
                lblSizer.Add(refAll,0,wx.ALIGN_CENTER_VERTICAL)               
                refAll.Bind(wx.EVT_BUTTON,onRefineAll)
            else:
                lblSizer.Add((-1,-1))

        i = -1
        for hist in DataArray:
            if hist == '_dataSource': continue
            i += 1
            if i == 1 and len(valList) > 2 and not deltaMode and CopyCtrl:
                but = wx.Button(Panel,wx.ID_ANY,'\u2192',style=wx.BU_EXACTFIT)
                but.valDict = {'arrays': valList, 'hists': histList,
                                  'dataSource':dataSource,
                                   'valEditList' :valEditList[row]}
                Sizer.Add(but,0,wx.ALIGN_CENTER_VERTICAL)
                but.Bind(wx.EVT_BUTTON,onSetAll)
            elif i == 1 and CopyCtrl:
                Sizer.Add((-1,-1))
            minval = None
            maxval = None
            # format the entry depending on what is defined
            if row not in DataArray[hist]:
                Sizer.Add((-1,-1))
                continue
            elif 'range' in DataArray[hist][row]:
                minval, maxval = DataArray[hist][row]['range']
            if ('init' in DataArray[hist][row] and
                      deltaMode and 'ref' in DataArray[hist][row]):
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                delta = arr[indx]
                arr,indx = DataArray[hist][row]['init']
                delta -= arr[indx]
                if abs(delta) < 9e-6: delta = 0.
                if delta == 0:
                    deltaS = ""
                else:
                    deltaS = f"\u0394 {delta:.4g} "
                valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                valrefsiz.Add(wx.StaticText(Panel,label=deltaS),0)
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['ref'])
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                valrefsiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL)
                checkButList[row].append(w)
                Sizer.Add(valrefsiz,0,
                                 wx.EXPAND|wx.ALIGN_RIGHT)
            elif 'init' in DataArray[hist][row] and deltaMode:
                # does this ever happen?
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                delta = arr[indx]
                arr,indx = DataArray[hist][row]['init']
                delta -= arr[indx]
                if delta == 0:
                    deltaS = ""
                else:
                    deltaS = f"\u0394 {delta:.4g} "
                Sizer.Add(wx.StaticText(Panel,label=deltaS),0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
            elif 'val' in DataArray[hist][row] and 'ref' in DataArray[hist][row]:
                valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                w = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                                         nDig=[9,7,'g'],
                                    xmin=minval,xmax=maxval)
                valEditList[row].append(w)
                valrefsiz.Add(w,0,WACV)
                if firstentry is None: firstentry = w
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['ref'])
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                valrefsiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL)
                checkButList[row].append(w)
                Sizer.Add(valrefsiz,0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
            elif 'val' in DataArray[hist][row]:
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                nDig = [9,7,'g']
                if type(arr[indx]) is str: nDig = None
                w = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                            nDig=nDig,
                            xmin=minval,xmax=maxval,notBlank=False)
                valEditList[row].append(w)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
                if firstentry is None: firstentry = w

            elif 'ref' in DataArray[hist][row]:
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['ref'])
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
                checkButList[row].append(w)
            elif 'str' in DataArray[hist][row]:
                val = indexArrayVal(dataSource,hist,DataArray[hist][row]['str'])
                if 'txt' in DataArray[hist][row]:
                    val = DataArray[hist][row]['txt']
                elif 'fmt' in DataArray[hist][row]:
                    f = DataArray[hist][row]['fmt']
                    val = f'{val:{f}}'
                Sizer.Add(wx.StaticText(Panel,label=val),0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
            else:
                print('Should not happen',DataArray[hist][row],hist,row)
    return firstentry,lblDict

def HistFrame(G2frame,Histograms):
    '''Put everything in a single FlexGridSizer.
    '''
    #---------------------------------------------------------------------
    # generate a dict with values for each histogram
    CopyCtrl = True
    global prmArray
    if G2frame.GroupInfo['displayMode'].startswith('Sample'):
        prmArray = getSampleVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Instrument'):
        prmArray = getInstVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Limits'):
        CopyCtrl = False
        prmArray = getLimitVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Background'):
        prmArray = getBkgVals(G2frame,Histograms)
        CopyCtrl = False
    else:
        prmArray = None
        print('Unexpected', G2frame.GroupInfo['displayMode'])
        return
    rowLabels = []
    lpos = 0
    nonZeroRows = []
    dataSource = prmArray['_dataSource']
    for hist in prmArray:
        if hist == '_dataSource': continue
        cols = len(prmArray)
        prevkey = None
        for key in prmArray[hist]:
            # find delta-terms that are non-zero
            if '\u0394' in G2frame.GroupInfo['displayMode']:
                if 'val' in prmArray[hist][key] and 'init' in prmArray[hist][key]:
                    arr,indx = prmArray[hist][key]['init']
                    val = indexArrayVal(dataSource,hist,prmArray[hist][key]['val'])
                    if abs(val-arr[indx]) > 1e-5: nonZeroRows.append(key)
            if key not in rowLabels:
                if prevkey is None:
                    rowLabels.insert(lpos,key)
                    lpos += 1
                else:
                    rowLabels.insert(rowLabels.index(prevkey)+1,key)
            prevkey = key
    # remove rows where delta-terms are all zeros
    if '\u0394' in G2frame.GroupInfo['displayMode']:
        rowLabels = [i for i in rowLabels if i in nonZeroRows]
    #=======  Generate GUI ===============================================
    # layout the window
    panel = midPanel = G2frame.dataWindow
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    G2G.HorizontalLine(mainSizer,panel)
    panel.SetSizer(mainSizer)
    deltaMode = "\u0394" in G2frame.GroupInfo['displayMode']
    n = 2
    if CopyCtrl and len(prmArray) > 2: n += 1 # add column for copy (when more than one histogram)
    valSizer = wx.FlexGridSizer(0,len(prmArray)+n-1,3,10)
    mainSizer.Add(valSizer,1,wx.EXPAND)
    valSizer.Add(wx.StaticText(midPanel,label=' '))
    valSizer.Add(wx.StaticText(midPanel,label=' Ref '))
    for i,hist in enumerate(histLabels(G2frame)[1]):
            if i == 1 and CopyCtrl:
                if deltaMode:
                    valSizer.Add((-1,-1))
                elif CopyCtrl:
                    valSizer.Add(wx.StaticText(midPanel,label=' Copy '))
            valSizer.Add(wx.StaticText(midPanel,
                        label=f"\u25A1 = {hist}"),
                             0,wx.ALIGN_CENTER)
    firstentry,lblDict = displayDataArray(rowLabels,prmArray,valSizer,midPanel,
                                      lblRow=True,
                                      deltaMode=deltaMode,CopyCtrl=CopyCtrl)
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
    indexDict = {'_dataSource':Histograms}
    def histderef(hist,l):
        a = Histograms[hist]
        for i in l:
            a = a[i]
        return a
   # loop over histograms in group
    for hist in groupDict[groupName]:
        indexDict[hist] = {}
        indexDict[hist]['Inst. name'] = {
            'val' : ('Sample Parameters','InstrName')}
        indexDict[hist]['Diff type'] = {
            'str' : ('Sample Parameters','Type')}
        indexDict[hist]['Scale factor'] = {
            'val' : ('Sample Parameters','Scale',0),
            'ref' : ('Sample Parameters','Scale',1),}
        # make a list of parameters to show
        histType = Histograms[hist]['Instrument Parameters'][0]['Type'][0]
        dataType = Histograms[hist]['Sample Parameters']['Type']
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
            if fmt is None and type(Histograms[hist]['Sample Parameters'][key]) is list:
                indexDict[hist][lbl] = {
                    'val' : ('Sample Parameters',key,0),
                    'ref' : ('Sample Parameters',key,1),}

            elif fmt is None:
                indexDict[hist][lbl] = {
                    'val' : ('Sample Parameters',key)}
            elif type(fmt) is str:
                indexDict[hist][lbl] = {
                    'str' : ('Sample Parameters',key),
                    'fmt' : fmt}
            
        for key in ('FreePrm1','FreePrm2','FreePrm3'):
            lbl = Controls[key]
            indexDict[hist][lbl] = {
                    'val' : ('Sample Parameters',key),
                     'rowlbl' : (Controls,key)
                }
    return indexDict

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
    indexDict = {'_dataSource':Histograms}
   # loop over histograms in group
    for hist in groupDict[groupName]:
        indexDict[hist] = {}
        insVal = Histograms[hist]['Instrument Parameters'][0]
        insType = insVal['Type'][1]
        if 'Bank' in Histograms[hist]['Instrument Parameters'][0]:
            indexDict[hist]['Bank'] = {
                'str' : ('Instrument Parameters',0,'Bank',1),
                'fmt' : '.0f'
                }
        indexDict[hist]['Hist type'] = {
                'str' : ('Instrument Parameters',0,'Type',1),
                }
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
                if fmt is None:
                    indexDict[hist][lbl] = {
                        'init' : (insVal[key],0),
                        'val' : ('Instrument Parameters',0,key,1),
                        'ref' : ('Instrument Parameters',0,key,2),}
                else:
                    indexDict[hist][lbl] = {
                        'str' : ('Instrument Parameters',0,key,1),
                        'fmt' : fmt
                }
    return indexDict

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
    indexDict = {'_dataSource':Histograms}
   # loop over histograms in group
    for hist in groupDict[groupName]:
        indexDict[hist] = {}
        for lbl,indx in [('Tmin',0),('Tmax',1)]:
            indexDict[hist][lbl] = {
                'val' : ('Limits',1,indx),
                'range': [Histograms[hist]['Limits'][0][0],
                          Histograms[hist]['Limits'][0][1]]
                }
        for i,item in enumerate(Histograms[hist]['Limits'][2:]):
            for l,indx in [('Low',0),('High',1)]:
                lbl = f'excl {l} {i+1}'
                indexDict[hist][lbl] = {
                    'val' : ('Limits',2+i,indx),
                    'range': [Histograms[hist]['Limits'][0][0],
                              Histograms[hist]['Limits'][0][1]]}
    return indexDict

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
    indexDict = {'_dataSource':Histograms}
    # loop over histograms in group
    for hist in groupDict[groupName]:
        indexDict[hist] = {}
        for lbl,indx,typ in [('Function',0,'str'),
                             ('ref flag',1,'ref'),
                             ('# Bkg terms',2,'str')]:
            indexDict[hist][lbl] = {
                typ : ('Background',0,indx)
                }
            if indx == 2:
                indexDict[hist][lbl]['fmt'] = '.0f'
        indexDict[hist]['# Debye terms'] = {
            'str' : ('Background',1,'nDebye'),
            'fmt' : '.0f'}
        for i,term in enumerate(Histograms[hist]['Background'][1]['debyeTerms']):
            for indx,l in enumerate(['A','R','U']):
                lbl = f'{l} #{i+1}'
                indexDict[hist][lbl] = {
                    'val' : ('Background',1,'debyeTerms',i,2*indx),
                    'ref' : ('Background',1,'debyeTerms',i,2*indx+1)}
        indexDict[hist]['# Bkg Peaks'] = {
            'str' : ('Background',1,'nPeaks'),
            'fmt' : '.0f'}
        for i,term in enumerate(Histograms[hist]['Background'][1]['peaksList']):
            for indx,l in enumerate(['pos','int','sig','gam']):
                lbl = f'{l} #{i+1}'
                indexDict[hist][lbl] = {
                    'val' : ('Background',1,'peaksList',i,2*indx),
                    'ref' : ('Background',1,'peaksList',i,2*indx+1)}
        if Histograms[hist]['Background'][1]['background PWDR'][0]:
            val = 'yes'
        else:
            val = 'no'
        indexDict[hist]['Fixed bkg file'] = {
            'str' : ('Background',1,'background PWDR',0),
            'txt' : val}
    return indexDict

def HAPframe(G2frame,Histograms,Phases):
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
        global prmArray
        prmArray = getHAPvals(G2frame,phaseList[page],Histograms,Phases)
        # construct a list of row labels, attempting to keep the
        # order they appear in the original array
        rowLabels = []
        lpos = 0
        for hist in prmArray:
            if hist == '_dataSource': continue
            prevkey = None
            for key in prmArray[hist]:
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
        lblSizer = wx.FlexGridSizer(0,2,hpad,2)
        lblScroll.SetSizer(lblSizer)
        bigSizer.Add(lblScroll,0,wx.EXPAND)

        # Create scrolled panel to display HAP data
        HAPScroll = wx.lib.scrolledpanel.ScrolledPanel(panel,
                        style=wx.VSCROLL|wx.HSCROLL|wx.ALWAYS_SHOW_SB)
        HAPSizer = wx.FlexGridSizer(0,len(prmArray),hpad,10)
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
        firstentry,lblDict = displayDataArray(rowLabels,prmArray,HAPSizer,HAPScroll,
                    lblRow=True,lblSizer=lblSizer,lblPanel=lblScroll)
        # get row sizes in data table
        HAPSizer.Layout()
        rowHeights = HAPSizer.GetRowHeights()
        # set row sizes in Labels
        # (must be done after HAPSizer row heights are defined)
        s = wx.Size(-1,rowHeights[0])
        w0.SetMinSize(s)
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

    G2G.HorizontalLine(mainSizer,midPanel)
    midPanel.SetSizer(mainSizer)
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

def getHAPvals(G2frame,phase,Histograms,Phases):
    '''Generate the Parameter Data Table (a dict of dicts) with 
    all HAP values for the selected phase and all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for HAP values.
    '''
    PhaseData = Phases[phase]
    SGData = PhaseData['General']['SGData']
    cell = PhaseData['General']['Cell'][1:]
    Amat,Bmat = G2lat.cell2AB(cell[:6])
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    groupName = G2frame.GroupInfo['groupName']
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    indexDict = {'_dataSource':PhaseData['Histograms']}
    for hist in groupDict[groupName]:
        indexDict[hist] = {}
        # phase fraction
        indexDict[hist]['Phase frac'] = {
            'val' : ('Scale',0),
            'ref' : ('Scale',1),}
        PhaseData['Histograms'][hist]['LeBail'] = PhaseData['Histograms'][hist].get('LeBail',False)
        indexDict[hist]['LeBail extract'] = {
            'str' : ('LeBail',),
            'txt' : "Yes" if PhaseData['Histograms'][hist]['LeBail'] else '(off)'}
        # size values
        if PhaseData['Histograms'][hist]['Size'][0] == 'isotropic':
            indexDict[hist]['Size'] = {
                'val' : ('Size',1,0),
                'ref' : ('Size',2,0),}
        elif PhaseData['Histograms'][hist]['Size'][0] == 'uniaxial':
            indexDict[hist]['Size/Eq'] = {
                'val' : ('Size',1,0),
                'ref' : ('Size',2,0),}
            indexDict[hist]['Size/Ax'] = {
                'val' : ('Size',1,1),
                'ref' : ('Size',2,1),}
            indexDict[hist]['Size/dir'] = {
                'str' : ('Size',3),
                'txt' : ','.join([str(i) for i in PhaseData['Histograms'][hist]['Size'][3]])}
        else:
            for i,lbl in enumerate(['S11','S22','S33','S12','S13','S23']):
                indexDict[hist][f'Size/{lbl}'] = {
                    'val' : ('Size',4,i),
                    'ref' : ('Size',5,i),}
        indexDict[hist]['Size LGmix'] = {
            'val' : ('Size',1,2),
            'ref' : ('Size',2,2),}
        # microstrain values
        if PhaseData['Histograms'][hist]['Mustrain'][0] == 'isotropic':
            indexDict[hist]['\u00B5Strain'] = {
                'val' : ('Mustrain',1,0),
                'ref' : ('Mustrain',2,0),}
        elif PhaseData['Histograms'][hist]['Mustrain'][0] == 'uniaxial':
            indexDict[hist]['\u00B5Strain/Eq'] = {
                'val' : ('Mustrain',1,0),
                'ref' : ('Mustrain',2,0),}
            indexDict[hist]['\u00B5Strain/Ax'] = {
                'val' : ('Mustrain',1,1),
                'ref' : ('Mustrain',2,1),}
            indexDict[hist]['\u00B5Strain/dir'] = {
                'str' : ('Mustrain',3),
                'txt' : ','.join([str(i) for i in PhaseData['Histograms'][hist]['Mustrain'][3]])}
        else:
            Snames = G2spc.MustrainNames(SGData)
            for i,lbl in enumerate(Snames):
                if i >= len(PhaseData['Histograms'][hist]['Mustrain'][4]): break
                indexDict[hist][f'\u00B5Strain/{lbl}'] = {
                    'val' : ('Mustrain',4,i),
                    'ref' : ('Mustrain',5,i),}
            muMean = G2spc.MuShklMean(SGData,Amat,PhaseData['Histograms'][hist]['Mustrain'][4][:len(Snames)])
            indexDict[hist]['\u00B5Strain/mean'] = {
                'str' : None,
                'txt' : f'{muMean:.2f}'}
        indexDict[hist]['\u00B5Strain LGmix'] = {
            'val' : ('Mustrain',1,2),
            'ref' : ('Mustrain',2,2),}

        # Hydrostatic terms
        Hsnames = G2spc.HStrainNames(SGData)
        for i,lbl in enumerate(Hsnames):
            if i >= len(PhaseData['Histograms'][hist]['HStrain'][0]): break
            indexDict[hist][f'Size/{lbl}'] = {
                'val' : ('HStrain',0,i),
                'ref' : ('HStrain',1,i),}

        # Preferred orientation terms
        if PhaseData['Histograms'][hist]['Pref.Ori.'][0] == 'MD':
            indexDict[hist]['March-Dollase'] = {
                'val' : ('Pref.Ori.',1),
                'ref' : ('Pref.Ori.',2),}
            indexDict[hist]['M-D/dir'] = {
                'str' : ('Pref.Ori.',3),
                'txt' : ','.join([str(i) for i in PhaseData['Histograms'][hist]['Pref.Ori.'][3]])}
        else:
            indexDict[hist]['Spherical harmonics'] = {
                'ref' : ('Pref.Ori.',2),}
            indexDict[hist]['SH order'] = {
                'str' : ('Pref.Ori.',4),
                'fmt' : '.0f'}
            for lbl in PhaseData['Histograms'][hist]['Pref.Ori.'][5]:
                indexDict[hist][f'SP {lbl}']= {
                    'val' : ('Pref.Ori.',5,lbl),
                    }
            indexDict[hist]['SH txtr indx'] = {
                'str' : None,
                'txt' : f'{G2lat.textureIndex(PhaseData['Histograms'][hist]['Pref.Ori.'][5]):.3f}'}
        # misc: Layer Disp, Extinction
        if 'Layer Disp' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Layer displ'] = {
                'val' : ('Layer Disp',0),
                'ref' : ('Layer Disp',1),}                
        if 'Extinction' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Extinction'] = {
                'val' : ('Extinction',0),
                'ref' : ('Extinction',1),}
        if 'Babinet' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Babinet A'] = {
                'val' : ('Babinet','BabA',0),
                'ref' : ('Babinet','BabA',1),}
        if 'Babinet' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Babinet U'] = {
                'val' : ('Babinet','BabU',0),
                'ref' : ('Babinet','BabU',1),}
    return indexDict
