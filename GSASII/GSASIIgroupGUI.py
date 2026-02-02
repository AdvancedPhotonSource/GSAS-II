# -*- coding: utf-8 -*-
'''
GSASIIgroupGUI: Routines for working with groups of histograms.

Outline of code in GSASIIgroupGUI:

* :func:`SearchGroups` is used to define which histograms are in each 
  group. Defines variables in Controls['Groups'] (see below)

* :func:`UpdateGroup` is used to create the contents of the DataWindow

* :func:`HistFrame` is called by :func:`UpdateGroup` to display the 
  main contents of the window 

* :func:`displayDataArray` is called by :func:`HistFrame` to create the 
  parameter table shown in the window, from the **Parameter Data Table** 
  entries, as described below.

A few routines are optionally used to color parameters by the shift
applied in the last refinement. Routine :func:`computeShiftInfo` computes
the shifts for all parameters, creates a set of steps for them and a color
scale. All these are placed in global array ``shiftInfo``. The 
function :func:`colorCaption` creates a caption showing the levels and colors
developed in :func:`computeShiftInfo`, while :func:`ColorTxtbox` applies
those colors to ValidatedTextCtrl entry widgets.

The settings in Controls entry ['Groups'] that define how histograms are 
divided into groups are:

* Controls['Groups']['groupDict'] 
   a dict where each key is the name of the group and the value is a list of 
   histograms in the group
* Controls['Groups']['notGrouped']
   a count of the number of histograms that are not in any group
* Controls['Groups']['template'] 
   the string used to set the grouping

** Parameter Data Table **

For use in creating GUI tables and to copy values between histograms, 
a dict is created in :func:`getSampleVals`, :func:`getInstVals`, 
:func:`getLimitVals`, :func:`getBkgVals` or :func:`getHAPvals` with 
the parameters used in that data tree item. The dict is organized 
with the following entries: 

 *  dict['__dataSource'] : SourceArray

 *  dict[histogram]['label'] : `innerdict`

where `label` is the text shown on the row label and `innerdict` is a dict 
containing one or more of the following elements:

 *  'val'    : (key1,key2,...)
 *  'ref'    : (key1, key2,...)
 *  'range'  : (float,float)
 *  'str'    : (key1,key2,...)
 *  'fmt'    : str
 *  'txt'    : str
 *  'init'   : float
 *  'rowlbl' : (array, key)
 *  'setintfunc' : function
 *  'prmname': str

One of 'val', 'ref' or 'str' elements must be present. 

 *  The 'val' tuple provides a reference to the float value for the 
    defined quantity, where SourceArray[histogram][key1][key2][...] 
    provides r/w access to the parameter. 

 *  The 'ref' tuple provides a reference to the bool value, where 
    SourceArray[histogram][key1][key2][...] provides r/w access to the 
    refine flag for the labeled quantity. 

    Both 'ref' and 'val' are usually defined together, but either may 
    occur alone. These exceptions will be for parameters where a single
    refine flag is used for a group of parameters or for non-refined 
    parameters.

    Note that 
    :func:`indexArrayVal` is used to obtain the value for the 
    quantity from the keys in the tuple for the 'ref' and 'val' entries,
    while :func:`indexArrayRef` 
    is used to obtain the dict/list reference and the key/index for that
    object from the saved info (so that the value can be changed).

 *  The 'str' value is something that cannot be edited from the GUI; If 'str' is
    present, the only other possible entries that can be present is either 'fmt' 
    or 'txt. 
    'str' is used for a parameter value that is typically computed or 
    cannot be edited in the group tables (edit by selecting the 
    appropriate in the histogram tree entry.

Optional elements in the array:

 *  The 'fmt' value is a string used to format a 'str' value to 
    convert it to a string. Used for float or int values.

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

 *  The 'setintfunc' value defines a function that is executed when the 
    'ref' tuple is changed. When this is supplied, a spin box is used 
    to set the value rather than a TextCtrl.

 *  The 'prmname' name entry lists the parameter name, as used in 
    varyList, parmDict etc. 

The full list of routines used in mod:`GSASIIgroupGUI` are:
'''

import math
import re
import wx

from . import GSASIIdataGUI as G2gd
from . import GSASIIspc as G2spc
from . import GSASIIlattice as G2lat
from . import GSASIIctrlGUI as G2G
from . import GSASIIpwdplot as G2pwpl

WACV = wx.ALIGN_CENTER_VERTICAL

shiftInfo = {'colorMode':False}  # home for into to color TextCtrls

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

def UpdateGroup(G2frame,selectedGrp,plot=True):
    '''Code to display parameters from all sets of grouped histograms
    '''
    def OnColorMode(event):
        if event.GetId() == G2G.wxID_GRPLOG:
            shiftInfo['colorMode'] = 'log'
        elif event.GetId() == G2G.wxID_GRPLIN:
            shiftInfo['colorMode'] = 'linear'
        else:
            shiftInfo['colorMode'] = False
        wx.CallAfter(UpdateGroup,G2frame,selectedGrp,False)
        
    def onDisplaySel(event):
        '''Change the type of parameters that are shown in response to the 
        pulldown in the upper left
        '''
        G2frame.GroupInfo['displayMode'] = dsplType.GetValue()
        wx.CallAfter(UpdateGroup,G2frame,selectedGrp,False)

    def onPhaseSel(event):
        '''Change the selected phase (in Hist/Phase mode only)
        '''
        G2frame.GroupInfo['selectedPhase'] = G2frame.GroupInfo['phaseSel'].GetValue()
        wx.CallAfter(UpdateGroup,G2frame,selectedGrp,False)

    def copyPrep():
        '''Prepare to copy parameters by locating groups that have the
        same number of histograms as the currently selected group and 
        then select from those groups. Initiated by copy all/copy selected 
        menu command.

        Could check even more closely that histograms match up properly,
        but for now I don't think that is needed.
        '''
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
        '''Copy all displayed parameters to the selected histograms.

        Called by menu command.
        '''
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

    def OnCopySel(event):
        '''Copy selected parameters from a list based on what is displayed 
        to the selected histograms.

        Called by menu command.
        '''
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
                                
    ##### beginning of code to display Groups of Histograms ###
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    if not hasattr(G2frame,'GroupInfo'):
        G2frame.GroupInfo = {}
    computeShiftInfo(G2frame) # Setup for color by shift/sigma
    # start with Sample but reuse the last displayed item
    G2frame.GroupInfo['displayMode'] = G2frame.GroupInfo.get('displayMode','Sample')
    #G2frame.GroupInfo['displayMode'] = G2frame.GroupInfo.get('displayMode','Hist/Phase') # debug code
    G2frame.GroupInfo['groupName'] = G2frame.GPXtree.GetItemText(selectedGrp)
    G2frame.GroupInfo['selectedPhase'] = G2frame.GroupInfo.get('selectedPhase',
                                                list(Phases.keys())[0])
    G2frame.GroupInfo['Redraw'] = (UpdateGroup,G2frame,selectedGrp,False)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.GroupMenu)
    G2frame.Bind(wx.EVT_MENU, OnCopyAll, id=G2G.wxID_GRPALL)
    G2frame.Bind(wx.EVT_MENU, OnCopySel, id=G2G.wxID_GRPSEL)
    for i in (G2G.wxID_GRPNONE, G2G.wxID_GRPLOG, G2G.wxID_GRPLIN):
        G2frame.Bind(wx.EVT_MENU, OnColorMode, id=i)
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
    G2frame.GroupInfo['phaseSel'] = wx.ComboBox(topParent,wx.ID_ANY,
                                value=G2frame.GroupInfo['selectedPhase'],
                                choices=list(Phases),
                                style=wx.CB_READONLY|wx.CB_DROPDOWN)
    G2frame.GroupInfo['phaseSel'].Bind(wx.EVT_COMBOBOX, onPhaseSel)
    topSizer.Add(G2frame.GroupInfo['phaseSel'],0,wx.RIGHT|wx.LEFT|WACV,10)
    topSizer.Add(wx.StaticText(topParent,
            label=f' Group "{histLabels(G2frame)[0]}" parameters'),
                     0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(topParent,helpIndex=G2frame.dataWindow.helpKey))

    HistFrame(G2frame,Histograms,Phases)
    if plot: G2pwpl.PlotPatterns(G2frame,plotType='GROUP')
    G2frame.dataWindow.SetDataSize()
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
    '''Find the array and key for a refinement flag. (From a 'ref' or 
    'val' entry.)

    :returns: arr, indx where arr[indx] is the refinement 
       flag value (bool for 'ref') or the value (float for 'val')
    '''
    indx = arrayIndices[-1]
    arr = dataSource[hist]
    for i in  arrayIndices[:-1]:
        arr = arr[i]
    return arr,indx

def indexArrayVal(dataSource,hist,arrayIndices):
    '''Find the value for a GSAS-II parameter in the data tree 
    from a 'val' or the refinement flag from the entry

    :returns: the parameter value (usually float or bool)
    '''
    if arrayIndices is None: return None
    arr = dataSource[hist]
    for i in  arrayIndices:
        arr = arr[i]
    return arr

ClearAllLbl = '☐'
SetAllLbl = '☑'
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
    all other edit widgets in that row.
    '''
    but = event.GetEventObject()
    firstVal = indexArrayVal(but.valDict['dataSource'],
                    but.valDict['hists'][0],but.valDict['arrays'][0])
    for f in but.valDict['valSetFxnList']:
        f(firstVal)

def onResetAll(event):
    '''Respond to the Reset button. Copies the initial setting (which 
    should be index 0) to the current setting (which should be index 1) 
    for all histograms in the current row. This is only used for 
    instrument parameters. 
    '''
    but = event.GetEventObject()
    dataSource = but.valDict['dataSource']
    for item,hist in zip(but.valDict['arrays'],but.valDict['hists']):
        arr,indx = indexArrayRef(dataSource,hist,item)
        if indx == 1:
            arr[1] = arr[0]
    G2frame = but.GetTopLevelParent()
    if G2frame.GroupInfo.get('Redraw'):
        wx.CallAfter(*G2frame.GroupInfo.get('Redraw'))

def displayDataArray(rowLabels,DataArray,Sizer,Panel,deltaMode=False,
                     CopyCtrl=True):
    '''Displays the data table in `DataArray` in Scrolledpanel `Panel`
    with wx.FlexGridSizer `Sizer`.
    '''
    def OnSpin(evt):
        '''respond when a SpinButton entry is changed (currently used for 
        background terms only)
        '''
        spin = (evt.GetEventObject())
        spin.arr[spin.indx] = evt.GetEventObject().GetValue()
        spin.txt.SetLabel(str(spin.arr[spin.indx]))
        spin.setTermsFnx()
    firstentry = None
    checkButList = {}
    valSetFxnList = {}
    lblDict = {}
    dataSource = DataArray['_dataSource']
    for row in rowLabels:
        checkButList[row] = []
        valSetFxnList[row] = []
        # show the row labels
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
            w = wx.StaticText(Panel,label=row)
        else: # used for "renameable" sample vars (str)
            w = G2G.ValidatedTxtCtrl(Panel,arr,key,size=(125,-1))
        Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
        lblDict[row] = w

        if deltaMode:
            resetBut = wx.Button(Panel,label='↩',style=wx.BU_EXACTFIT)
            resetBut.valDict = {'arrays': valList, 'hists': histList,
                                  'dataSource':dataSource}
            Sizer.Add(resetBut,0,wx.ALIGN_CENTER_VERTICAL)               
            resetBut.Bind(wx.EVT_BUTTON,onResetAll)

        if len(refList) > 2: # >two refinement flags in this row. Include ste/clear all button
            lbl = SetAllLbl
            if all([indexArrayVal(dataSource,hist,i) for i in refList]): lbl = ClearAllLbl
            refAll = wx.Button(Panel,label=lbl,style=wx.BU_EXACTFIT)
            font = refAll.GetFont()
            font.PointSize += 5
            refAll.SetFont(font)
            refAll.refDict = {'arrays': refList, 'hists': histList,
                              'dataSource':dataSource}
            refAll.checkButList = checkButList[row]
            Sizer.Add(refAll,0,wx.ALIGN_CENTER_VERTICAL)               
            refAll.Bind(wx.EVT_BUTTON,onRefineAll)
        else:
            Sizer.Add((-1,-1))

        i = -1
        for hist in DataArray:
            if hist == '_dataSource': continue
            i += 1
            if i == 1 and len(valList) > 2 and not deltaMode and CopyCtrl:
                # copy button; place after 0th column 
                but = wx.Button(Panel,wx.ID_ANY,'\u2192',style=wx.BU_EXACTFIT)
                but.valDict = {'arrays': valList, 'hists': histList,
                                  'dataSource':dataSource,
                                   'valSetFxnList' :valSetFxnList[row]}
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
                # instrument parameter (has 'init' value) and in "delta"
                # mode: show difference
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
                # a non-refined instrument parameter that has an initial value
                # I don't think this ever happens
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
            elif 'setintfunc' in DataArray[hist][row]:
                # integer value that has a function to be called when changed
                # GUI: provide a Spin Box to increment or decrement
                # Currently used for background terms only.
                spinsiz = wx.BoxSizer(wx.HORIZONTAL)
                spin = wx.SpinButton(Panel, wx.ID_ANY)
                spin.Bind(wx.EVT_SPIN, OnSpin)
                spin.txt = wx.StaticText(Panel,label='?')                
                spinsiz.Add(spin.txt,0,wx.ALIGN_CENTER_VERTICAL)
                spinsiz.Add((5,-1))
                spinsiz.Add(spin,0,wx.ALIGN_CENTER_VERTICAL)
                Sizer.Add(spinsiz,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)

                spin.SetRange(1,36)
                spin.arr,spin.indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                spin.setTermsFnx = DataArray[hist][row]['setintfunc']
                spin.SetValue(spin.arr[spin.indx])
                spin.txt.SetLabel(str(spin.arr[spin.indx]))
                def SetVal(newval,spin=spin):
                    'Used to set a value for the current spinbutton & associated StaticText'
                    spin.arr[spin.indx] = newval
                    spin.SetValue(spin.arr[spin.indx])
                    spin.txt.SetLabel(str(newval))
                    spin.setTermsFnx()
                valSetFxnList[row].append(SetVal)
            elif 'val' in DataArray[hist][row] and 'ref' in DataArray[hist][row]:
                # a value that has a refine flag
                valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                txtctrl = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                                         nDig=[9,7,'g'],
                                    xmin=minval,xmax=maxval)
                valSetFxnList[row].append(txtctrl.ChangeValue)
                valrefsiz.Add(txtctrl,0,WACV)
                if firstentry is None: firstentry = txtctrl
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['ref'])
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                valrefsiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL)
                checkButList[row].append(w)
                Sizer.Add(valrefsiz,0,
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
                ColorTxtbox(txtctrl,DataArray[hist][row],row)
            elif 'val' in DataArray[hist][row]:
                # a value that is not refined, but can be edited
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['val'])
                nDig = [9,7,'g']
                if type(arr[indx]) is str: nDig = None
                w = G2G.ValidatedTxtCtrl(Panel,arr,indx,size=(80,-1),
                            nDig=nDig,
                            xmin=minval,xmax=maxval,notBlank=False)
                valSetFxnList[row].append(w.ChangeValue)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
                if firstentry is None: firstentry = w
                ColorTxtbox(w,DataArray[hist][row],row)

            elif 'ref' in DataArray[hist][row]:
                # a refinement flag for a group of parameters
                arr,indx = indexArrayRef(dataSource,hist,DataArray[hist][row]['ref'])
                w = G2G.G2CheckBox(Panel,'',arr,indx)
                Sizer.Add(w,0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
                checkButList[row].append(w)
            elif 'str' in DataArray[hist][row]:
                # a value that is not refined and cannot be edited
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

def HistFrame(G2frame,Histograms,Phases):
    '''Display all the parameters of the selected type. 
    This puts everything in a single FlexGridSizer.
    '''
    #---------------------------------------------------------------------
    # generate a dict with values for each histogram
    CopyCtrl = True
    global prmArray
    G2frame.GroupInfo['phaseSel'].Show(False)
    if G2frame.GroupInfo['displayMode'].startswith('Sample'):
        prmArray = getSampleVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Instrument'):
        prmArray = getInstVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Limits'):
        CopyCtrl = False
        prmArray = getLimitVals(G2frame,Histograms)
    elif G2frame.GroupInfo['displayMode'].startswith('Background'):
        prmArray = getBkgVals(G2frame,Histograms)
        CopyCtrl = True
    elif G2frame.GroupInfo['displayMode'].startswith('Hist/Phase') and Phases:
        phase = G2frame.GroupInfo.get('selectedPhase',list(Phases.keys())[0])
        prmArray = getHAPvals(G2frame,phase,Histograms,Phases)
        G2frame.GroupInfo['phaseSel'].Show(True)
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
        #cols = len(prmArray)
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
    n = 2 # extra columns
    # remove rows where delta-terms are all zeros
    if '\u0394' in G2frame.GroupInfo['displayMode']:
        rowLabels = [i for i in rowLabels if i in nonZeroRows]
        n += 1
    #=======  Generate GUI ===============================================
    # layout the window
    midPanel = G2frame.dataWindow
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    G2G.HorizontalLine(mainSizer,midPanel)
    midPanel.SetSizer(mainSizer)
    deltaMode = "\u0394" in G2frame.GroupInfo['displayMode']
    if CopyCtrl and len(prmArray) > 2: n += 1 # add column for copy (when more than one histogram)
    if deltaMode and not rowLabels:
        mainSizer.Add(wx.StaticText(midPanel,
                                    label='No parameters deviate from initial values'),
                      0,wx.ALL,10)
    elif not rowLabels:  # I don't think this can happen
        mainSizer.Add(wx.StaticText(midPanel,label='No parameters to display'))
    else:
        valSizer = wx.FlexGridSizer(0,len(prmArray)+n-1,3,10)
        mainSizer.Add(valSizer,1,wx.EXPAND)
        # place column headers into table
        valSizer.Add(wx.StaticText(midPanel,label=' '))
        if '\u0394' in G2frame.GroupInfo['displayMode']:
            valSizer.Add(wx.StaticText(midPanel,label='Reset'))
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
        # fill the remaining rows in the table
        firstentry,lblDict = displayDataArray(rowLabels,prmArray,valSizer,midPanel,
                                      deltaMode=deltaMode,CopyCtrl=CopyCtrl)
        
        if firstentry is not None:    # prevent scroll to show last entry
            wx.Window.SetFocus(firstentry)
            firstentry.SetInsertionPoint(0) # prevent selection of text in widget
        colorCaption(G2frame)

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
        prmPrefix = f':{Histograms[hist]['hId']}:'
        indexDict[hist] = {}
        indexDict[hist]['Inst. name'] = {
            'val' : ('Sample Parameters','InstrName')}
        indexDict[hist]['Diff type'] = {
            'str' : ('Sample Parameters','Type')}
        indexDict[hist]['Scale factor'] = {
            'val' : ('Sample Parameters','Scale',0),
            'ref' : ('Sample Parameters','Scale',1),
            'prmname': prmPrefix+'Scale',
            }
        # make a list of parameters to show
        histType = Histograms[hist]['Instrument Parameters'][0]['Type'][0]
        dataType = Histograms[hist]['Sample Parameters']['Type']
        if histType[2] in ['A','B','C']:
            parms.append(['Gonio. radius','Gonio radius','.3f'])
        #if 'PWDR' in histName:
        if dataType == 'Debye-Scherrer':
            if 'T' in histType:
                parms += [['Absorption','Sample abs, \xb5r/\u03bb',None,'Absorption']]
            else:
                parms += [
                    ['DisplaceX','Sample X displ',None,'DisplaceX'],
                    ['DisplaceY','Sample Y displ',None,'DisplaceY'],
                    ['Absorption','Sample abs,\xb5\xb7r',None,'Absorption']]
        elif dataType == 'Bragg-Brentano':
            parms += [
                ['Shift','Sample displ',None,'Shift'],
                ['Transparency','Sample transp',None,'Transparency'],
                ['SurfRoughA','Surf rough A',None,'SurfRoughA'],
                ['SurfRoughB','Surf rough B',None,'SurfRoughB']]
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
        for entry in parms:
            prmNam = None
            if len(entry) == 3:
                key,lbl,fmt = entry
            else:
                key,lbl,fmt,prmNam = entry
            if fmt is None and type(Histograms[hist]['Sample Parameters'][key]) is list:
                indexDict[hist][lbl] = {
                    'val' : ('Sample Parameters',key,0),
                    'ref' : ('Sample Parameters',key,1),}
                if prmNam:
                    indexDict[hist][lbl]['prmname'] = prmPrefix+prmNam

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
        prmPrefix = f':{Histograms[hist]['hId']}:'
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
                    # for instrument parameters assume prmname is based on key
                    indexDict[hist][lbl]['prmname'] = prmPrefix+key
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
        prmPrefix = f':{Histograms[hist]['hId']}:'
        indexDict[hist] = {}
        for lbl,indx,typ in [('Function type',0,'str'),
                             ('bkg refine flag',1,'ref'),
                             ('# Bkg terms',2,'val')]:
            indexDict[hist][lbl] = {
                typ : ('Background',0,indx)
                }
            if indx == 2:
                indexDict[hist][lbl]['fmt'] = '.0f'
                def OnChangeBkgTerms(Histograms=Histograms,hist=hist):
                    'set the number of terms to match the new number'
                    nterms = Histograms[hist]['Background'][0][2]
                    Histograms[hist]['Background'][0][3:] = (
                        Histograms[hist]['Background'][0][3:] +
                        36*[0.0])[:nterms]
                indexDict[hist][lbl]['setintfunc'] = OnChangeBkgTerms
        indexDict[hist]['# Debye terms'] = {
            'str' : ('Background',1,'nDebye'),
            'fmt' : '.0f'}
        for i,term in enumerate(Histograms[hist]['Background'][1]['debyeTerms']):
            for indx,l in enumerate(['A','R','U']):
                lbl = f'{l} #{i+1}'
                indexDict[hist][lbl] = {
                    'val' : ('Background',1,'debyeTerms',i,2*indx),
                    'ref' : ('Background',1,'debyeTerms',i,2*indx+1),
                    'prmname': f'{prmPrefix}Debye{l};{i}'}
        indexDict[hist]['# Bkg Peaks'] = {
            'str' : ('Background',1,'nPeaks'),
            'fmt' : '.0f'}
        for i,term in enumerate(Histograms[hist]['Background'][1]['peaksList']):
            for indx,l in enumerate(['pos','int','sig','gam']):
                lbl = f'{l} #{i+1}'
                indexDict[hist][lbl] = {
                    'val' : ('Background',1,'peaksList',i,2*indx),
                    'ref' : ('Background',1,'peaksList',i,2*indx+1),
                    'prmname': f'{prmPrefix}BkPk{l};{i}'}
        if Histograms[hist]['Background'][1]['background PWDR'][0]:
            val = 'yes'
        else:
            val = 'no'
        indexDict[hist]['Fixed bkg file'] = {
            'str' : ('Background',1,'background PWDR',0),
            'txt' : val}
        if Histograms[hist]['Background'][1]['background PWDR'][0]:
            indexDict[hist]['Fixed bkg mult'] = {
                'ref' : ('Background',1,'background PWDR',2),
                'val' : ('Background',1,'background PWDR',1),
                'prmname': f'{prmPrefix}BF mult'}
    return indexDict

def getHAPvals(G2frame,phase,Histograms,Phases):
    '''Generate the Parameter Data Table (a dict of dicts) with 
    all HAP values for the selected phase and all histograms in the 
    selected histogram group (from G2frame.GroupInfo['groupName']).
    This will be used to generate the contents of the GUI for HAP values.
    '''
    PhaseData = Phases[phase]
    pId = Phases[phase]['pId']
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
        prmPrefix = f'{pId}:{Histograms[hist]['hId']}:'
        indexDict[hist] = {}
        # phase fraction
        indexDict[hist]['Phase frac'] = {
            'val' : ('Scale',0),
            'ref' : ('Scale',1),
            'prmname': prmPrefix+'Scale',
            }
        PhaseData['Histograms'][hist]['LeBail'] = PhaseData['Histograms'][hist].get('LeBail',False)
        indexDict[hist]['LeBail extract'] = {
            'str' : ('LeBail',),
            'txt' : "Yes" if PhaseData['Histograms'][hist]['LeBail'] else '(off)'}
        # size values
        if PhaseData['Histograms'][hist]['Size'][0] == 'isotropic':
            indexDict[hist]['Size'] = {
                'val' : ('Size',1,0),
                'ref' : ('Size',2,0),
                'prmname': prmPrefix+'Size;i',
                }
        elif PhaseData['Histograms'][hist]['Size'][0] == 'uniaxial':
            indexDict[hist]['Size/Eq'] = {
                'val' : ('Size',1,0),
                'ref' : ('Size',2,0),
                'prmname': prmPrefix+'Size;i',}
            indexDict[hist]['Size/Ax'] = {
                'val' : ('Size',1,1),
                'ref' : ('Size',2,1),
                'prmname': prmPrefix+'Size;a',}
            indexDict[hist]['Size/dir'] = {
                'str' : ('Size',3),
                'txt' : ','.join([str(i) for i in PhaseData['Histograms'][hist]['Size'][3]])}
        else:
            for i,lbl in enumerate(['S11','S22','S33','S12','S13','S23']):
                indexDict[hist][f'Size/{lbl}'] = {
                    'val' : ('Size',4,i),
                    'ref' : ('Size',5,i),
                    'prmname': prmPrefix+f'Size;{i}',}
        indexDict[hist]['Size LGmix'] = {
            'val' : ('Size',1,2),
            'ref' : ('Size',2,2),
            'prmname': prmPrefix+'Size;mx',}
        # microstrain values
        if PhaseData['Histograms'][hist]['Mustrain'][0] == 'isotropic':
            indexDict[hist]['\u00B5Strain'] = {
                'val' : ('Mustrain',1,0),
                'ref' : ('Mustrain',2,0),
                'prmname': prmPrefix+'Mustrain;i',}
        elif PhaseData['Histograms'][hist]['Mustrain'][0] == 'uniaxial':
            indexDict[hist]['\u00B5Strain/Eq'] = {
                'val' : ('Mustrain',1,0),
                'ref' : ('Mustrain',2,0),
                'prmname': prmPrefix+'Mustrain;i',}
            indexDict[hist]['\u00B5Strain/Ax'] = {
                'val' : ('Mustrain',1,1),
                'ref' : ('Mustrain',2,1),
                'prmname': prmPrefix+'Mustrain;a',}
            indexDict[hist]['\u00B5Strain/dir'] = {
                'str' : ('Mustrain',3),
                'txt' : ','.join([str(i) for i in PhaseData['Histograms'][hist]['Mustrain'][3]])}
        else:
            Snames = G2spc.MustrainNames(SGData)
            for i,lbl in enumerate(Snames):
                if i >= len(PhaseData['Histograms'][hist]['Mustrain'][4]): break
                indexDict[hist][f'\u00B5Strain/{lbl}'] = {
                    'val' : ('Mustrain',4,i),
                    'ref' : ('Mustrain',5,i),
                    'prmname': prmPrefix+f'Mustrain;{i}',}
            muMean = G2spc.MuShklMean(SGData,Amat,PhaseData['Histograms'][hist]['Mustrain'][4][:len(Snames)])
            indexDict[hist]['\u00B5Strain/mean'] = {
                'str' : None,
                'txt' : f'{muMean:.2f}'}
        indexDict[hist]['\u00B5Strain LGmix'] = {
            'val' : ('Mustrain',1,2),
            'ref' : ('Mustrain',2,2),
            'prmname': prmPrefix+'Mustrain;mx',}

        # Hydrostatic terms
        Hsnames = G2spc.HStrainNames(SGData)
        for i,lbl in enumerate(Hsnames):
            if i >= len(PhaseData['Histograms'][hist]['HStrain'][0]): break
            indexDict[hist][f'Recip Cell/{lbl}'] = {
                'val' : ('HStrain',0,i),
                'ref' : ('HStrain',1,i),
                'prmname': prmPrefix+f'{lbl}',}
        # Preferred orientation terms
        if PhaseData['Histograms'][hist]['Pref.Ori.'][0] == 'MD':
            indexDict[hist]['March-Dollase'] = {
                'val' : ('Pref.Ori.',1),
                'ref' : ('Pref.Ori.',2),
                'prmname': prmPrefix+'MD',}
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
                    'prmname': prmPrefix+f'{lbl}',}
            indexDict[hist]['SH txtr indx'] = {
                'str' : None,
                'txt' : f'{G2lat.textureIndex(PhaseData["Histograms"][hist]["Pref.Ori."][5]):.3f}'}
        # misc: Layer Disp, Extinction
        if 'Layer Disp' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Layer displ'] = {
                'val' : ('Layer Disp',0),
                'ref' : ('Layer Disp',1),
                'prmname': prmPrefix+'LayerDisp'}
        if 'Extinction' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Extinction'] = {
                'val' : ('Extinction',0),
                'ref' : ('Extinction',1),
                'prmname': prmPrefix+'Extinction',}
        if 'Babinet' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Babinet A'] = {
                'val' : ('Babinet','BabA',0),
                'ref' : ('Babinet','BabA',1),
                'prmname': prmPrefix+'BabA',}
        if 'Babinet' in PhaseData['Histograms'][hist]:
            indexDict[hist]['Babinet U'] = {
                'val' : ('Babinet','BabU',0),
                'ref' : ('Babinet','BabU',1),
                'prmname': prmPrefix+'BabU',}
    return indexDict

def computeShiftInfo(G2frame):
    '''Define the contents of the array used to color 
    ValidatedTextCtrls by shift/sigma
    '''
    covdata = G2frame.GPXtree.GetItemPyData(
                  G2gd.GetGPXtreeItemId(G2frame,
                                        G2frame.root,'Covariance'))
    shiftOsig = {}
    if len(covdata['varyList']) == len(covdata['sig']) == len(covdata['Lastshft']):
        for var,sig,shift in zip(covdata['varyList'],covdata['sig'],covdata['Lastshft']):
            try:
                shiftOsig[var] = float(shift/sig)
            except:
                pass
    else:
        print('Warning: in Covariance, varyList, sig and Lastshft do not have same lengths!')
    vmin = None
    vmax = 1.0 # shifts less than 1 are small
    for v in shiftOsig.values():
        v = abs(v)
        if vmin is None:
            vmin = v
        if vmax is None:
            vmax = v
        vmin = min(vmin,v)
        vmax = max(vmax,v)
    shiftLinearLevels = []
    shiftLogLevels = []
    shiftColors = []
    if shiftOsig: # come up with levels and colors
        shiftSteps = 5
        for i in range(shiftSteps):
            shiftLinearLevels.append(vmin+i*(vmax-vmin)/shiftSteps)
        lmin = math.log(max(vmin,1e-8))
        lmax = math.log(max(vmax,1e-8))
        for i in range(shiftSteps):
            shiftLogLevels.append(math.exp(lmin + i*(lmax-lmin)/shiftSteps))
        for i in range(shiftSteps):
            #g = int(255.5-i*255/(shiftSteps-1))
            #shiftColors.append((255,g,0,255))  # yellow to red
            #shiftColors.append((255-g,g,0,255))  # green to red
            # generate light pink to dark pink
            g = int(230.5-i*130/(shiftSteps-1))
            shiftColors.append((255,g,g,255))  # white to red
    shiftInfo['shiftLinearLevels'] = shiftLinearLevels
    shiftInfo['shiftLogLevels'] = shiftLogLevels
    shiftInfo['shiftColors'] = shiftColors
    shiftInfo['colorCaption'] = colorCaption
    shiftInfo['shiftOsig'] = shiftOsig

def colorCaption(G2frame):
    'Create legend info showing the range of shifts'
    colorSizer = G2frame.dataWindow.bottomBox
    colorSizer.Clear(True)
    cp = G2frame.dataWindow.bottomPanel
    if shiftInfo['colorMode'] == 'log':
        levels = shiftInfo['shiftLogLevels']
    elif shiftInfo['colorMode'] == 'linear':
        levels = shiftInfo['shiftLinearLevels']
    else:
        return None
    colorSizer.Add(wx.StaticText(cp,label='Shift/\u03c3 by color:'),0,WACV)
    for i,(level,color) in enumerate(zip(levels,shiftInfo['shiftColors'])):
        colorSizer.Add((15,-1))
        #b = wx.StaticText(cp,label=f'{i}: >{level:.4g}')
        b = wx.StaticText(cp,label=f'>{level:.4g}')
        b.SetBackgroundColour(color)
        colorSizer.Add(b,0,WACV)

def ColorTxtbox(txtctrl,DataArray,row):
    '''Color a ValidatedTextCtrl textbox based on shift-over-sigma values
    '''
    if shiftInfo['colorMode'] == 'log':
        levels = shiftInfo['shiftLogLevels']
    elif shiftInfo['colorMode'] == 'linear':
        levels = shiftInfo['shiftLinearLevels']
    else:
        return
    if 'prmname' in DataArray:
        prmname = DataArray['prmname']
        if prmname in shiftInfo['shiftOsig']:
            val = abs(shiftInfo['shiftOsig'][prmname])            
            for (level,color) in reversed(list(
                zip(levels,shiftInfo['shiftColors']))):
                if val >= level:
                    txtctrl.SaveBackgroundColor(color)
                    break
            else: # unexpected
                txtctrl.SaveBackgroundColor(color)
