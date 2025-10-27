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
'''

# import math
# import os
# import re
# import copy
# import platform
# import pickle
# import sys
# import random as ran

import numpy as np
# import numpy.ma as ma
import wx

from . import GSASIIpath
from . import GSASIIdataGUI as G2gd
# from . import GSASIIobj as G2obj
# import GSASIIpwdGUI as G2pdG
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

def UpdateGroup(G2frame,item):
    G2gd.SetDataMenuBar(G2frame)
    G2frame.dataWindow.helpKey = "Groups/Powder"
    # topSizer = G2frame.dataWindow.topBox
    # parent = G2frame.dataWindow.topPanel
    # topSizer.Add(wx.StaticText(parent,label=' Group edit goes here someday'),0,WACV)
    # topSizer.Add((-1,-1),1,wx.EXPAND)
    # topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    # G2G.HorizontalLine(G2frame.dataWindow.GetSizer(),G2frame.dataWindow)
    # #G2frame.dataWindow.GetSizer().Add(text,1,wx.ALL|wx.EXPAND)
    G2frame.groupName = G2frame.GPXtree.GetItemText(item)
    HAPframe(G2frame)
    G2pwpl.PlotPatterns(G2frame,plotType='GROUP')
    
def histLabels(G2frame):
    # find portion of hist name that is the same and different
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.groupName
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    h0 = groupDict[groupName][0]
    msk = [True] * len(h0)
    for h in groupDict[groupName][1:]:
        msk = [m & (h0i == hi) for h0i,hi,m in zip(h0,h,msk)]
    # place rectangular box in the loc of non-common letter(s)
    commonltrs = ''.join([h0i if m else '\u25A1' for (h0i,m) in zip(h0,msk)])
    #for i,h in enumerate(groupDict[groupName]):
    # make list with histogram name unique letters   
    histlbls = [''.join([hi for (hi,m) in zip(h,msk) if not m])
                   for h in groupDict[groupName]]
    return commonltrs,histlbls

def HAPframe(G2frame):
    def OnPageChanged(event):
        'respond to a notebook tab'
        if event:
            page = event.GetSelection()
            print('page selected',page,phaseList[page])
        else:
            page = 0
            print('no page selected',phaseList[page])
        HAPtable = getHAPvals(G2frame,phaseList[page])
        #for hist in HAPtable:
        #    showTable(phaseList[page],hist,HAPtable[hist])
        # clear out old widgets
        for panel in HAPtabs:
            if panel.GetSizer():
                panel.GetSizer().Destroy()                
        panel = HAPtabs[page]
        HAPSizer = wx.FlexGridSizer(0,len(HAPtable)+1,2,10)
        # construct a list of row labels, attempting to keep the
        # order they appear in the original array
        rowsLbls = []
        lpos = 0
        for hist in HAPtable:
            prevkey = None
            for key in HAPtable[hist]:
                if key not in rowsLbls:
                    if prevkey is None:
                        rowsLbls.insert(lpos,key)
                        lpos += 1
                    else:
                        rowsLbls.insert(rowsLbls.index(prevkey)+1,key)
                prevkey = key
        # label columns with histograms
        common,hLbl = histLabels(G2frame)
        HAPSizer.Add((-1,-1))
        for hist in hLbl:
            HAPSizer.Add(wx.StaticText(panel,label=f"\u25A1 = {hist}"),0,
                             wx.ALIGN_CENTER)
        for row in rowsLbls:
            HAPSizer.Add(wx.StaticText(panel,label=row),0,WACV)
            for hist in HAPtable:
                if row not in HAPtable[hist]:
                    HAPSizer.Add((-1,-1))
                elif 'val' in HAPtable[hist][row] and 'ref' in HAPtable[hist][row]:
                    valrefsiz = wx.BoxSizer(wx.HORIZONTAL)
                    arr,indx = HAPtable[hist][row]['ref']
                    valrefsiz.Add(G2G.G2CheckBox(panel,'',arr,indx),0,WACV)
                    arr,indx = HAPtable[hist][row]['val']
                    valrefsiz.Add(G2G.ValidatedTxtCtrl(panel,
                                        arr,indx,size=(75,-1)),0,WACV)
                    HAPSizer.Add(valrefsiz,0,
                                     wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
                elif 'val' in HAPtable[hist][row]:
                    arr,indx = HAPtable[hist][row]['val']
                    HAPSizer.Add(G2G.ValidatedTxtCtrl(panel,
                                        arr,indx,size=(75,-1)),0,
                                     wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)

                elif 'ref' in HAPtable[hist][row]:
                    arr,indx = HAPtable[hist][row]['ref']
                    HAPSizer.Add(G2G.G2CheckBox(panel,'',arr,indx),0,
                                     wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
                elif 'str' in HAPtable[hist][row]:
                    HAPSizer.Add(wx.StaticText(panel,
                                    label=HAPtable[hist][row]['str']),0,
                                    wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER)
#                                    wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
                else:
                    print('Should not happen',HAPtable[hist][row],hist,row)
        panel.SetSizer(HAPSizer)
        HAPSizer.Fit(panel)
        panel.SetScrollbars(1, 1, HAPSizer.GetMinSize().width,
                                  HAPSizer.GetMinSize().height)

        #breakpoint()
        
        #wx.CallAfter(FillDDataWindow,page)

    G2frame.dataWindow.ClearData()
    
    topSizer = G2frame.dataWindow.topBox
    topParent = G2frame.dataWindow.topPanel
    botSizer = G2frame.dataWindow.bottomBox
    botParent = G2frame.dataWindow.bottomPanel
    common,_ = histLabels(G2frame)
    topSizer.Add(wx.StaticText(topParent,label=f'HAP parameters for group "{common}"'),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(topParent,helpIndex=G2frame.dataWindow.helpKey))

    # based on GSASIIddataGUI.MakeHistPhaseWin
    midPanel = G2frame.dataWindow
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    G2G.HorizontalLine(mainSizer,midPanel)
    midPanel.SetSizer(mainSizer)
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree() # reindex
    if not Phases:
        mainSizer.Add(wx.StaticText(midPanel,label='There are no phases in use'))
        G2frame.dataWindow.SetDataSize()
        return
    HAPBook = G2G.GSNoteBook(parent=midPanel)
    mainSizer.Add(HAPBook,1,wx.ALL|wx.EXPAND,1)
    HAPtabs = []
    phaseList = []
    for phaseName in Phases:
        phaseList.append(phaseName)
        HAPtabs.append(wx.ScrolledWindow(HAPBook))
        HAPBook.AddPage(HAPtabs[-1],phaseName)
    HAPBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)
    #G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.DataMenu)
    # fill the 'Select tab' menu
    # mid = G2frame.dataWindow.DataMenu.FindMenu('Select tab')
    # menu = G2frame.dataWindow.DataMenu.GetMenu(mid)
    # items = menu.GetMenuItems()
    # for item in items:
    #      menu.Remove(item)
    # if len(phaseList) == 0: return
    # for i,page in enumerate(phaseList):
    #     Id = wx.NewId()
    #     if menu.FindItem(page) >= 0: continue # is tab already in menu?
    #     menu.Append(Id,page,'')
    #     TabSelectionIdDict[Id] = page
    #     G2frame.Bind(wx.EVT_MENU, OnSelectPage, id=Id)
    page = 0
    HAPBook.SetSelection(page)
    OnPageChanged(None)
    #wx.CallAfter(FillDDataWindow,page)
    # done 
    G2frame.dataWindow.SetDataSize()

def getHAPvals(G2frame,phase=None):
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
        return

    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    groupName = G2frame.groupName
    groupDict = Controls.get('Groups',{}).get('groupDict',{})
    if groupName not in groupDict:
        print(f'Unexpected: {groupName} not in groupDict')
        return
    # loop over histograms in group
    parmDict = {}
    for hist in groupDict[groupName]:
        parmDict[hist] = makeHAPtbl(G2frame,phase,PhaseData,hist,Controls)
        #showTable(phase,hist,parmDict[hist])
    return parmDict

def makeHAPtbl(G2frame,phase,PhaseData,hist,Controls):
    '''Construct a dict pointing to the contents of the HAP
    variables. The contents of the dict will be:

       'label' : innerdict

    where innerdict can contain the following elements:

       'val' : (array, key)
       'ref' : (array, key)
       'str' : string

    One of these will be present. 

    The 'str' value is something that cannot be edited; If 'str' is
    present, it will be the only entry in innerdict. 

    The 'val' tuple provides a reference to the float value for the 
    defined quantity, array[key]

    The 'ref' tuple provides a reference to the bool value, array[key]
    for the refine flag defined quantity

    Both 'ref' and 'val' are usually defined together, but either may 
    occur alone.

    :return: the dict, as described above.
    '''
    SGData = PhaseData['General']['SGData']
    cell = PhaseData['General']['Cell'][1:]
    Amat,Bmat = G2lat.cell2AB(cell[:6])
    G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, hist)
    data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
    HAPdict = PhaseData['Histograms'][hist]

    parmDict = {}
    # phase fraction
    parmDict['Phase fraction'] = {
        'val' : (HAPdict['Scale'],0),
        'ref' : (HAPdict['Scale'],1),}

    parmDict['LeBail extraction'] = {
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
        parmDict['muStrain'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
    elif arr[0] == 'uniaxial':
        parmDict['muStrain/Eq'] = {
            'val' : (arr[1],0),
            'ref' : (arr[2],0),}
        parmDict['muStrain/Ax'] = {
            'val' : (arr[1],1),
            'ref' : (arr[2],1),}
        parmDict['muStrain/dir'] = {
            'str' : ','.join([str(i) for i in arr[3]])}
    else:
        Snames = G2spc.MustrainNames(SGData)
        for i,lbl in enumerate(Snames):
            if i >= len(arr[4]): break
            parmDict[f'muStrain/{lbl}'] = {
                'val' : (arr[4],i),
                'ref' : (arr[5],i),}
        muMean = G2spc.MuShklMean(SGData,Amat,arr[4][:len(Snames)])
        parmDict['muStrain/mean'] = {
            'str' : f'{muMean:.2f}'}
    parmDict['muStrain LGmix'] = {
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
        parmDict['Layer displacement'] = {
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

def showTable(phase,hist,parmDict):
    # show the variables and values
    print(phase,hist)
    for sel in parmDict:
        arr = parmDict[sel]
        val = 'N/A'
        if 'val' in arr:
            val = arr['val'][0] [arr['val'][1]]
        if 'str' in arr:
            val = f'"{arr['str']}"'
        ref = 'N/A'
        if 'ref' in arr:
            ref = arr['ref'][0] [arr['ref'][1]]
        print(f'{sel!r:20s}  {val}  {ref}')
    print('\n')    
    #break
