"""
Created on Sun Jul  6 14:27:24 2025

@author: vondreele
"""

# -*- coding: utf-8 -*-
# GSASII - phase data display routines
"""

Routines for Phase/RMC follow. Only Update routines are here
all others are in GSASIIphsGUI.py
"""
import copy
import os
import shutil
import time

import numpy as np
import wx
import wx.grid as wg

from .. import GSASIIElem as G2elem
from .. import GSASIIlattice as G2lat
from .. import GSASIImath as G2mth
from .. import GSASIIplot as G2plt
from .. import GSASIIpwd as G2pwd
from .. import GSASIIspc as G2spc
from ..config import atmdata
from . import GSASIIctrlGUI as G2G
from . import GSASIIdataGUI as G2gd
from . import GSASIIphsGUI as G2phsG

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

try:
    VERY_LIGHT_GREY = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
    WHITE = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
    BLACK = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT)
    RED = wx.Colour(255, 0, 0)
    WACV = wx.ALIGN_CENTER_VERTICAL
except:
    pass
mapDefault = G2elem.mapDefault
TabSelectionIdDict = {}
# trig functions in degrees
sind = lambda x: np.sin(x * np.pi / 180.0)
tand = lambda x: np.tan(x * np.pi / 180.0)
cosd = lambda x: np.cos(x * np.pi / 180.0)
asind = lambda x: 180.0 * np.arcsin(x) / np.pi
acosd = lambda x: 180.0 * np.arccos(x) / np.pi
atan2d = lambda x, y: 180.0 * np.arctan2(y, x) / np.pi
is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
sqt2 = np.sqrt(2.0)
sqt3 = np.sqrt(3.0)

# previous rigid body selections
prevResId = None
prevVecId = None
prevSpnId = None

GkDelta = chr(0x0394)
Angstr = chr(0x00C5)

RMCmisc = {}
ranDrwDict = {}
ranDrwDict["atomList"] = []
DrawStyleChoice = [
    " ",
    "lines",
    "vdW balls",
    "sticks",
    "balls & sticks",
    "ellipsoids",
    "polyhedra",
]
#### UpdateRMC GUI ################################################################################
# # fullrmc stuff TODO:
# #  1) need to implement swapping in scripts
# #  2) fullrmc tutorials


def UpdateRMC(G2frame, data):
    """Present the controls for running fullrmc, RMCProfile or PDFfit"""
    global runFile

    def OnRMCselect(event):
        G2frame.RMCchoice = RMCsel.GetStringSelection()
        wx.CallLater(200, UpdateRMC, G2frame, data)

    def GetAtmChoice(pnl, RMCPdict):
        Indx = {}

        def OnAtSel(event):
            Obj = event.GetEventObject()
            itype = Indx[Obj.GetId()]
            tid = RMCPdict["atSeq"].index(Obj.GetStringSelection())
            if itype < nTypes:
                if itype == tid:
                    tid += 1
                RMCPdict["atSeq"] = G2lat.SwapItems(RMCPdict["atSeq"], itype, tid)
            Pairs = []
            atSeq = RMCPdict["atSeq"]
            lenA = len(atSeq)
            for pair in [
                [" %s-%s" % (atSeq[i], atSeq[j]) for j in range(i, lenA)]
                for i in range(lenA)
            ]:
                #                for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
                Pairs += pair
            RMCPdict["Pairs"] = {pairs: [0.0, 0.0, 0.0] for pairs in Pairs}
            if RMCPdict["useBVS"]:
                BVSpairs = []
                for pair in [
                    [" %s-%s" % (atSeq[i], atSeq[j]) for j in range(i + 1, lenA)]
                    for i in range(lenA)
                ]:
                    BVSpairs += pair
                RMCPdict["BVS"] = {pairs: [0.0, 0.0, 0.0, 0.0] for pairs in BVSpairs}
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnValSel(event):
            Obj = event.GetEventObject()
            itype = Indx[Obj.GetId()]
            RMCPdict["Oxid"][itype][0] = Obj.GetStringSelection()
            wx.CallAfter(UpdateRMC, G2frame, data)

        nTypes = len(RMCPdict["aTypes"])
        atmChoice = wx.FlexGridSizer(nTypes + 1, 5, 5)
        atmChoice.Add(wx.StaticText(pnl, label="atom ordering: "), 0, WACV)
        for iType in range(nTypes):
            atChoice = RMCPdict["atSeq"][iType:]
            atmSel = wx.ComboBox(
                pnl, choices=atChoice, style=wx.CB_DROPDOWN | wx.TE_READONLY
            )
            atmSel.SetStringSelection(RMCPdict["atSeq"][iType])
            atmSel.Bind(wx.EVT_COMBOBOX, OnAtSel)
            Indx[atmSel.GetId()] = iType
            atmChoice.Add(atmSel, 0, WACV)
        if RMCPdict["useBVS"]:
            atmChoice.Add(wx.StaticText(pnl, label="Valence: "), 0, WACV)
            for itype in range(nTypes):
                valChoice = atmdata.BVSoxid[RMCPdict["atSeq"][itype]]
                valSel = wx.ComboBox(
                    pnl, choices=valChoice, style=wx.CB_DROPDOWN | wx.TE_READONLY
                )
                try:
                    valSel.SetStringSelection(RMCPdict["Oxid"][itype][0])
                except IndexError:
                    RMCPdict["Oxid"].append([RMCPdict["atSeq"][itype], 0.0])
                valSel.Bind(wx.EVT_COMBOBOX, OnValSel)
                Indx[valSel.GetId()] = itype
                atmChoice.Add(valSel, 0, WACV)
            atmChoice.Add(wx.StaticText(pnl, label="BVS weight: "), 0, WACV)
            for itype in range(nTypes):
                atmChoice.Add(
                    G2G.ValidatedTxtCtrl(pnl, RMCPdict["Oxid"][itype], 1, xmin=0.0),
                    0,
                    WACV,
                )
        if G2frame.RMCchoice == "RMCProfile":
            atmChoice.Add(wx.StaticText(pnl, label="max shift: "), 0, WACV)
            for iType in range(nTypes):
                atId = RMCPdict["atSeq"][iType]
                atmChoice.Add(
                    G2G.ValidatedTxtCtrl(
                        pnl, RMCPdict["aTypes"], atId, xmin=0.0, xmax=1.0
                    ),
                    0,
                    WACV,
                )
            atmChoice.Add(wx.StaticText(pnl, label="Isotope: "), 0, WACV)
            for iType in range(nTypes):
                atId = RMCPdict["atSeq"][iType]
                try:
                    lbl = RMCPdict["Isotope"][atId]
                except:
                    lbl = "?"
                atmChoice.Add(wx.StaticText(pnl, label=lbl), 0, WACV)
        return atmChoice

    def GetSwapSizer(RMCPdict):
        def OnDelSwap(event):
            Obj = event.GetEventObject()
            swap = Indx[Obj.GetId()]
            del RMCPdict["Swaps"][swap]
            wx.CallAfter(UpdateRMC, G2frame, data)

        Indx = {}
        atChoice = RMCPdict["atSeq"]
        # if G2frame.RMCchoice == 'fullrmc':
        #     atChoice = atNames
        swapSizer = wx.FlexGridSizer(6, 5, 5)
        swapLabels = [" ", "Atom-A", "Atom-B", " Swap prob.", " ", "delete"]
        for lab in swapLabels:
            swapSizer.Add(wx.StaticText(G2frame.FRMC, label=lab), 0, WACV)
        for ifx, swap in enumerate(RMCPdict["Swaps"]):
            swapSizer.Add((20, -1))
            for i in [0, 1]:
                if swap[i] not in atChoice:
                    swap[i] = atChoice[0]
                atmSel = G2G.EnumSelector(G2frame.FRMC, swap, i, atChoice)
                swapSizer.Add(atmSel, 0, WACV)
            swapSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC, swap, 2, xmin=0.01, xmax=0.5, size=(50, 25)
                ),
                0,
                WACV,
            )
            swapSizer.Add((20, -1))
            delBtn = wx.Button(G2frame.FRMC, label="Del", style=wx.BU_EXACTFIT)
            delBtn.Bind(wx.EVT_BUTTON, OnDelSwap)
            Indx[delBtn.GetId()] = ifx
            swapSizer.Add(delBtn, 0, WACV)
        return swapSizer

    def GetPairSizer(pnl, RMCPdict):
        pairSizer = wx.FlexGridSizer(len(RMCPdict["Pairs"]) + 1, 5, 5)
        pairSizer.Add((5, 5), 0)
        for pair in RMCPdict["Pairs"]:
            pairSizer.Add(wx.StaticText(pnl, label=pair), 0, WACV)
        if G2frame.RMCchoice == "RMCProfile":
            pairSizer.Add(wx.StaticText(pnl, label="%14s" % " Hard min: "), 0, WACV)
            for pair in RMCPdict["Pairs"]:
                pairSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        pnl,
                        RMCPdict["Pairs"][pair],
                        0,
                        xmin=0.0,
                        xmax=10.0,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            pairSizer.Add(wx.StaticText(pnl, label="%14s" % " Search from: "), 0, WACV)
        elif G2frame.RMCchoice == "fullrmc":
            pairSizer.Add(wx.StaticText(pnl, label="%14s" % " Distance min: "), 0, WACV)
        for pair in RMCPdict["Pairs"]:
            pairSizer.Add(
                G2G.ValidatedTxtCtrl(
                    pnl, RMCPdict["Pairs"][pair], 1, xmin=0.0, xmax=10.0, size=(50, 25)
                ),
                0,
                WACV,
            )
        pairSizer.Add(wx.StaticText(pnl, label="%14s" % "to: "), 0, WACV)
        for pair in RMCPdict["Pairs"]:
            pairSizer.Add(
                G2G.ValidatedTxtCtrl(
                    pnl, RMCPdict["Pairs"][pair], 2, xmin=0.0, xmax=10.0, size=(50, 25)
                ),
                0,
                WACV,
            )
        return pairSizer

    def GetMetaSizer(RMCPdict, metalist):
        metaSizer = wx.FlexGridSizer(0, 2, 5, 5)
        for item in metalist:
            metaSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" Metadata item: " + item + " "),
                0,
                WACV,
            )
            metaSizer.Add(
                G2G.ValidatedTxtCtrl(G2frame.FRMC, RMCPdict["metadata"], item), 0, WACV
            )
        return metaSizer

    def SetRestart(invalid, value, tc):
        RMCPdict["ReStart"] = [True, True]

    def GetSuperSizer(RMCPdict, Xmax):
        superSizer = wx.BoxSizer(wx.HORIZONTAL)
        axes = ["X", "Y", "Z"]
        for i, ax in enumerate(axes):
            superSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" %s-axis: " % ax), 0, WACV
            )
            superSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC,
                    RMCPdict["SuperCell"],
                    i,
                    xmin=1,
                    xmax=Xmax,
                    size=(50, 25),
                    OnLeave=SetRestart,
                ),
                0,
                WACV,
            )
        return superSizer

    def FileSizer(RMCPdict):
        def OnFileSel(event):
            Obj = event.GetEventObject()
            fil = Indx[Obj.GetId()]
            G2frame.OnFileSave(event)
            dlg = wx.FileDialog(
                G2frame.FRMC,
                "Choose " + fil,
                G2G.GetImportPath(G2frame),
                style=wx.FD_OPEN,
                wildcard=fil + "(*.*)|*.*",
            )
            if dlg.ShowModal() == wx.ID_OK:
                fpath, fName = os.path.split(dlg.GetPath())
                if os.path.exists(
                    fName
                ):  # is there a file by this name in the current directory?
                    RMCPdict["files"][fil][0] = fName
                else:  # nope, copy it
                    # TODO: is G2frame.LastGPXdir the right choice here or
                    #       do I want the current working directory (same?)
                    shutil.copy(dlg.GetPath(), os.path.join(G2frame.LastGPXdir, fName))
                if not os.path.exists(fName):  # sanity check
                    print(
                        f"Error: file {fName} not found in .gpx directory ({G2frame.LastGPXdir})"
                    )
                    return
                G2frame.LastImportDir = fpath  # set so next file is found in same place
                dlg.Destroy()
                RMCPdict["ReStart"][0] = True
                if G2frame.RMCchoice == "PDFfit":
                    start = 0
                    XY = np.empty((1, 2))
                    while XY.shape[0] == 1:
                        try:
                            XY = np.loadtxt(fName, skiprows=start)
                        except ValueError:
                            start += 1
                    name = "Ndata"
                    if "X" in fil:
                        name = "Xdata"
                    RMCPdict[name]["Datarange"][0] = np.min(XY.T[0])
                    RMCPdict[name]["Datarange"][1] = np.max(XY.T[0])
                    RMCPdict[name]["Fitrange"][1] = np.max(XY.T[0])
            else:
                dlg.Destroy()

            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnFileFormat(event):
            Obj = event.GetEventObject()
            fil = Indx[Obj.GetId()]
            RMCPdict["files"][fil][3] = Obj.GetStringSelection()

        def OnPlotBtn(event):
            Obj = event.GetEventObject()
            fil = Indx[Obj.GetId()]
            fileItem = RMCPdict["files"][fil]
            start = 0
            XY = np.empty((1, 2))
            while XY.shape[0] == 1:
                try:
                    XY = np.loadtxt(fileItem[0], skiprows=start)
                except ValueError:
                    start += 1
                    if start > 500:  # absurd number of header lines!
                        wx.MessageBox(
                            "WARNING: %s has bad data at end;\n RMCProfile may fail to read it"
                            % fileItem[0],
                            style=wx.ICON_ERROR,
                        )
                        break
            Xlab = "Q"
            if "G(R)" in fileItem[2].upper():
                Xlab = "R"
            G2plt.PlotXY(
                G2frame,
                [
                    XY.T[:2],
                ],
                labelX=Xlab,
                labelY=fileItem[2],
                newPlot=True,
                Title=fileItem[0],
                lines=True,
            )

        def OnCorrChk(event):
            Obj = event.GetEventObject()
            fil = Indx[Obj.GetId()]
            RMCPdict["files"][fil][3] = not RMCPdict["files"][fil][3]

        def OnDelBtn(event):
            Obj = event.GetEventObject()
            fil = Indx[Obj.GetId()]
            RMCPdict["files"][fil][0] = "Select"
            RMCPdict["ReStart"][0] = True
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnRef(event):
            Obj = event.GetEventObject()
            name, item = Indx[Obj.GetId()]
            RMCPdict[name][item][1] = not RMCPdict[name][item][1]

        def OnRefSel(event):
            RMCPdict["refinement"] = reftype.GetStringSelection()
            wx.CallLater(100, UpdateRMC, G2frame, data)

        def OnDataSel(event):
            RMCPdict["SeqDataType"] = dataType.GetStringSelection()

        def OnSeqCopy(event):
            RMCPdict["SeqCopy"] = not RMCPdict["SeqCopy"]

        def OnSeqReverse(event):
            RMCPdict["SeqReverse"] = not RMCPdict["SeqReverse"]

        # --- FileSizer starts here
        Indx = {}
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if G2frame.RMCchoice == "PDFfit":
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            reftype = wx.RadioBox(
                G2frame.FRMC,
                label="PDFfit refinement type:",
                choices=["normal", "sequential"],
            )
            reftype.SetStringSelection(RMCPdict.get("refinement", "normal"))
            reftype.Bind(wx.EVT_RADIOBOX, OnRefSel)
            topSizer.Add(reftype)
            if "seq" in RMCPdict.get("refinement", "normal"):
                dataType = wx.RadioBox(
                    G2frame.FRMC, label="Seq data type:", choices=["X", "N"]
                )
                dataType.SetStringSelection(RMCPdict.get("SeqDataType", "X"))
                dataType.Bind(wx.EVT_RADIOBOX, OnDataSel)
                topSizer.Add(dataType)
                endSizer = wx.BoxSizer(wx.VERTICAL)
                seqcopy = wx.CheckBox(G2frame.FRMC, label=" Copy to next")
                seqcopy.SetValue(RMCPdict["SeqCopy"])
                seqcopy.Bind(wx.EVT_CHECKBOX, OnSeqCopy)
                endSizer.Add(seqcopy)
                seqreverse = wx.CheckBox(G2frame.FRMC, label=" Reverse processing")
                seqreverse.SetValue(RMCPdict["SeqReverse"])
                seqreverse.Bind(wx.EVT_CHECKBOX, OnSeqReverse)
                endSizer.Add(seqreverse)
                topSizer.Add(endSizer, 0, WACV)
            mainSizer.Add(topSizer)
        elif G2frame.RMCchoice == "fullrmc":
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add(
                wx.StaticText(
                    G2frame.FRMC,
                    label='  Select data for processing (files must be 2 columns w/headers preceeded by "#"; edit if needed)',
                )
            )
            mainSizer.Add(topSizer)
            Heads = ["Name", "File", "type", "Plot", "Delete"]
            fileSizer = wx.FlexGridSizer(5, 5, 5)
            Formats = ["RMC", "GUDRUN", "STOG"]
            for head in Heads:
                fileSizer.Add(wx.StaticText(G2frame.FRMC, label=head), 0, WACV)
            for fil in RMCPdict["files"]:
                fileSizer.Add(wx.StaticText(G2frame.FRMC, label=fil), 0, WACV)
                Rfile = RMCPdict["files"][fil][0]
                filSel = wx.Button(G2frame.FRMC, label=Rfile)
                filSel.Bind(wx.EVT_BUTTON, OnFileSel)
                Indx[filSel.GetId()] = fil
                fileSizer.Add(filSel, 0, WACV)
                if Rfile and os.path.exists(
                    Rfile
                ):  # in case .gpx file is moved away from G(R), F(Q), etc. files
                    # fileSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['files'][fil],1,size=(50,25)),0,WACV)
                    # patch
                    if len(RMCPdict["files"][fil]) < 4:
                        RMCPdict["files"][fil].append(0)
                    if len(RMCPdict["files"][fil]) < 5:
                        RMCPdict["files"][fil].append(True)
                    # end patch
                    if "G(r)" in fil:
                        choices = "G(r)-RMCProfile", "G(r)-PDFFIT", "g(r)"
                        if type(RMCPdict["files"][fil][3]) is bool:
                            RMCPdict["files"][fil][3] = 0
                        fmtTyp = G2G.G2ChoiceButton(
                            G2frame.FRMC, choices, RMCPdict["files"][fil], 3
                        )
                    elif "(Q)" in fil:
                        choices = "F(Q)-RMCProfile", "S(Q)-PDFFIT"
                        if type(RMCPdict["files"][fil][3]) is bool:
                            RMCPdict["files"][fil][3] = 0
                        fmtTyp = G2G.G2ChoiceButton(
                            G2frame.FRMC, choices, RMCPdict["files"][fil], 3
                        )
                    else:
                        fmtTyp = (-1, -1)
                    fileSizer.Add(fmtTyp, 0, WACV)
                    plotBtn = wx.Button(
                        G2frame.FRMC, label="Plot", style=wx.BU_EXACTFIT
                    )
                    plotBtn.Bind(wx.EVT_BUTTON, OnPlotBtn)
                    Indx[plotBtn.GetId()] = fil
                    fileSizer.Add(plotBtn, 0, WACV)
                    delBtn = wx.Button(G2frame.FRMC, label="Del", style=wx.BU_EXACTFIT)
                    delBtn.Bind(wx.EVT_BUTTON, OnDelBtn)
                    Indx[delBtn.GetId()] = fil
                    fileSizer.Add(delBtn, 0, WACV)
                    if "(Q)" in fil:
                        fileSizer.Add((-1, -1), 0)
                        corrChk = wx.CheckBox(
                            G2frame.FRMC, label="Apply sinc convolution? "
                        )
                        corrChk.SetValue(RMCPdict["files"][fil][4])
                        Indx[corrChk.GetId()] = fil
                        corrChk.Bind(wx.EVT_CHECKBOX, OnCorrChk)
                        fileSizer.Add(corrChk, 0, WACV)
                        # fileSizer.Add((-1,-1),0)
                        fileSizer.Add((-1, -1), 0)
                        fileSizer.Add((-1, -1), 0)
                        fileSizer.Add((-1, -1), 0)
                elif "Select" not in Rfile:  # file specified, but must not exist
                    RMCPdict["files"][fil][0] = "Select"  # set filSel?
                    fileSizer.Add(
                        wx.StaticText(
                            G2frame.FRMC,
                            label="Warning: file not found.\nWill be removed",
                        ),
                        0,
                    )
                    fileSizer.Add((-1, -1), 0)
                    fileSizer.Add((-1, -1), 0)
                else:
                    RMCPdict["files"][fil][0] = "Select"  # set filSel?
                    # fileSizer.Add((-1,-1),0)
                    fileSizer.Add((-1, -1), 0)
                    fileSizer.Add((-1, -1), 0)
                    fileSizer.Add((-1, -1), 0)
            mainSizer.Add(fileSizer, 0)
            return mainSizer

        if G2frame.RMCchoice == "PDFfit" and RMCPdict["refinement"] == "sequential":

            def OnAddPDF(event):
                """Add PDF G(r)s while maintanining original sequence"""
                usedList = RMCPdict["seqfiles"]
                PDFlist = [item[1:][0] for item in G2frame.GetFileList("PDF")]
                PDFdict = dict([item[1:] for item in G2frame.GetFileList("PDF")])
                PDFnames = [
                    item for item in PDFdict if item not in [itm[0] for itm in usedList]
                ]
                dlg = G2G.G2MultiChoiceDialog(
                    G2frame.FRMC,
                    "Add PDF dataset",
                    "Select G(r) data to use in seq. PDFfit",
                    PDFnames,
                )
                if dlg.ShowModal() == wx.ID_OK:
                    PDFuse = dlg.GetSelections()
                    for item in PDFuse:
                        pId = G2gd.GetGPXtreeItemId(
                            G2frame, G2frame.root, PDFnames[item]
                        )
                        data = G2frame.GPXtree.GetItemPyData(
                            G2gd.GetGPXtreeItemId(G2frame, pId, "PDF Controls")
                        )
                        try:
                            insrt = PDFlist.index(PDFnames[item]) - 1
                            RMCPdict["seqfiles"].insert(
                                insrt + 1, [PDFnames[item], data]
                            )
                        except ValueError:
                            RMCPdict["seqfiles"].append([PDFnames[item], data])
                dlg.Destroy()
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnDelPDF(event):
                usedList = [item[0] for item in RMCPdict["seqfiles"]]
                dlg = G2G.G2MultiChoiceDialog(
                    G2frame.FRMC,
                    "Delete PDF dataset",
                    "Select G(r) data to delete frpm seq. PDFfit",
                    usedList,
                )
                if dlg.ShowModal() == wx.ID_OK:
                    PDFdel = dlg.GetSelections()
                    PDFdel.reverse()
                    for item in PDFdel:
                        del RMCPdict["seqfiles"][item]
                dlg.Destroy()
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnSetColVal(event):
                parms = {
                    "Rmin": [0.01, 5.0],
                    "Rmax": [5.0, 30.0],
                    "dscale": [0.5, 2.0],
                    "qdamp": [0.0, 0.5],
                    "qbroad": [0.0, 0.1],
                    "Temp": 300,
                }
                c = event.GetCol()
                if c >= 0:
                    if c in [3, 5, 7]:
                        seqGrid.ClearSelection()
                        seqGrid.SelectCol(c, True)
                        if seqGrid.GetColLabelValue(c) != "refine":
                            return
                        choice = [
                            "Y - vary all",
                            "N - vary none",
                        ]
                        dlg = wx.SingleChoiceDialog(
                            G2frame,
                            "Select refinement option for "
                            + seqGrid.GetColLabelValue(c - 1),
                            "Refinement controls",
                            choice,
                        )
                        dlg.CenterOnParent()
                        if dlg.ShowModal() == wx.ID_OK:
                            sel = dlg.GetSelection()
                            varib = colLabels[c - 1]
                            if sel == 0:
                                for row in range(seqGrid.GetNumberRows()):
                                    RMCPdict["seqfiles"][row][1][varib][1] = True
                            else:
                                for row in range(seqGrid.GetNumberRows()):
                                    RMCPdict["seqfiles"][row][1][varib][1] = False
                    elif c in [0, 1, 2, 4, 6, 8]:
                        seqGrid.ClearSelection()
                        seqGrid.SelectCol(c, True)
                        parm = colLabels[c]
                        dlg = G2G.SingleFloatDialog(
                            G2frame,
                            "New value",
                            "Enter value for " + parm,
                            0.0,
                            parms[parm],
                        )
                        if dlg.ShowModal() == wx.ID_OK:
                            value = dlg.GetValue()
                            if c in [2, 4, 6]:
                                for row in range(seqGrid.GetNumberRows()):
                                    RMCPdict["seqfiles"][row][1][parm][0] = value
                            elif c == 8:
                                for row in range(seqGrid.GetNumberRows()):
                                    RMCPdict["seqfiles"][row][1][parm] = value
                            else:
                                for row in range(seqGrid.GetNumberRows()):
                                    RMCPdict["seqfiles"][row][1]["Fitrange"][c] = value
                    wx.CallAfter(UpdateRMC, G2frame, data)

            def OnSetVal(event):
                r, c = event.GetRow(), event.GetCol()
                if c >= 0:
                    if c in [3, 5, 7]:
                        varib = colLabels[c - 1]
                        RMCPdict["seqfiles"][r][1][varib][1] = bool(
                            seqGrid.GetCellValue(r, c)
                        )
                    elif c in [0, 1, 2, 4, 6, 8]:
                        parm = colLabels[c]
                        if c in [2, 4, 6]:
                            RMCPdict["seqfiles"][r][1][parm][0] = float(
                                seqGrid.GetCellValue(r, c)
                            )
                        elif c == 8:
                            RMCPdict["seqfiles"][r][1][parm] = float(
                                seqGrid.GetCellValue(r, c)
                            )
                        else:
                            RMCPdict["seqfiles"][r][1]["Fitrange"][c] = float(
                                seqGrid.GetCellValue(r, c)
                            )

            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add(
                wx.StaticText(G2frame.FRMC, label="  Select data for processing: ")
            )
            mainSizer.Add(topSizer)
            G2frame.GetStatusBar().SetStatusText(
                'NB: All PDFs used in sequential PDFfit must be the same type ("X" or "N") - there is no check',
                1,
            )
            if "seqfiles" not in RMCPdict:
                RMCPdict["seqfiles"] = []
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add(
                wx.StaticText(
                    G2frame.FRMC, label=" Sequential data list for PDFfit:  "
                ),
                0,
                WACV,
            )
            addPDF = wx.Button(G2frame.FRMC, label="Add PDF G(r) data sets")
            addPDF.Bind(wx.EVT_BUTTON, OnAddPDF)
            topSizer.Add(addPDF, 0, WACV)
            delPDF = wx.Button(G2frame.FRMC, label="Delete PDF G(r) data sets")
            delPDF.Bind(wx.EVT_BUTTON, OnDelPDF)
            topSizer.Add(delPDF, 0, WACV)
            mainSizer.Add(topSizer)
            table = [
                [
                    item[1]["Fitrange"][0],
                    item[1]["Fitrange"][1],
                    item[1]["dscale"][0],
                    item[1]["dscale"][1],
                    item[1]["qdamp"][0],
                    item[1]["qdamp"][1],
                    item[1]["qbroad"][0],
                    item[1]["qbroad"][1],
                    item[1].get("Temp", 300.0),
                ]
                for item in RMCPdict["seqfiles"]
            ]
            colLabels = [
                "Rmin",
                "Rmax",
                "dscale",
                "refine",
                "qdamp",
                "refine",
                "qbroad",
                "refine",
                "Temp",
            ]
            rowLabels = [item[0] for item in RMCPdict["seqfiles"]]
            Types = [
                wg.GRID_VALUE_FLOAT + ":10,2",
                wg.GRID_VALUE_FLOAT + ":10,2",
                wg.GRID_VALUE_FLOAT + ":10,4",
                wg.GRID_VALUE_BOOL,
                wg.GRID_VALUE_FLOAT + ":10,4",
                wg.GRID_VALUE_BOOL,
                wg.GRID_VALUE_FLOAT + ":10,4",
                wg.GRID_VALUE_BOOL,
                wg.GRID_VALUE_FLOAT + ":10,2",
            ]
            seqTable = G2G.Table(
                table, rowLabels=rowLabels, colLabels=colLabels, types=Types
            )
            seqGrid = G2G.GSGrid(G2frame.FRMC)
            seqGrid.SetTable(seqTable, True)
            seqGrid.AutoSizeColumns(True)
            seqGrid.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, OnSetColVal)
            seqGrid.Bind(wg.EVT_GRID_CELL_CHANGED, OnSetVal)
            mainSizer.Add(seqGrid)
            return mainSizer

        # begin FileSizer
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(
            wx.StaticText(G2frame.FRMC, label="  Select data for processing: ")
        )
        mainSizer.Add(topSizer)
        # RMCProfile & PDFfit (Normal)
        Heads = ["Name", "File", "Format", "Weight", "Plot", "Delete"]
        fileSizer = wx.FlexGridSizer(6, 5, 5)
        Formats = ["RMC", "GUDRUN", "STOG"]
        for head in Heads:
            fileSizer.Add(wx.StaticText(G2frame.FRMC, label=head), 0, WACV)
        for fil in RMCPdict["files"]:
            for head in Heads:
                fileSizer.Add(wx.StaticText(G2frame.FRMC, label=20 * "-"), 0, WACV)
            fileSizer.Add(wx.StaticText(G2frame.FRMC, label=fil), 0, WACV)
            Rfile = RMCPdict["files"][fil][0]
            filSel = wx.Button(G2frame.FRMC, label=Rfile)
            filSel.Bind(wx.EVT_BUTTON, OnFileSel)
            Indx[filSel.GetId()] = fil
            fileSizer.Add(filSel, 0, WACV)
            nform = 3
            Name = "Ndata"
            if "Xray" in fil:
                nform = 1
                Name = "Xdata"
            if Rfile and os.path.exists(
                Rfile
            ):  # incase .gpx file is moved away from G(R), F(Q), etc. files
                fileFormat = wx.ComboBox(
                    G2frame.FRMC,
                    choices=Formats[:nform],
                    style=wx.CB_DROPDOWN | wx.TE_READONLY,
                )
                fileFormat.SetStringSelection(RMCPdict["files"][fil][3])
                Indx[fileFormat.GetId()] = fil
                fileFormat.Bind(wx.EVT_COMBOBOX, OnFileFormat)
                fileSizer.Add(fileFormat, 0, WACV)
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(G2frame.FRMC, RMCPdict["files"][fil], 1),
                    0,
                    WACV,
                )
                plotBtn = wx.Button(G2frame.FRMC, label="Plot?", style=wx.BU_EXACTFIT)
                plotBtn.Bind(wx.EVT_BUTTON, OnPlotBtn)
                Indx[plotBtn.GetId()] = fil
                fileSizer.Add(plotBtn, 0, WACV)
                delBtn = wx.Button(G2frame.FRMC, label="Del", style=wx.BU_EXACTFIT)
                delBtn.Bind(wx.EVT_BUTTON, OnDelBtn)
                Indx[delBtn.GetId()] = fil
                fileSizer.Add(delBtn, 0, WACV)
            else:
                RMCPdict["files"][fil][0] = "Select"
                fileSizer.Add((5, 5), 0)
                fileSizer.Add((5, 5), 0)
                fileSizer.Add((5, 5), 0)
                fileSizer.Add((5, 5), 0)
            if "Select" not in Rfile and "PDFfit" in G2frame.RMCchoice:
                fileSizer.Add(
                    wx.StaticText(G2frame.FRMC, label=" R-range (from/to)"), 0, WACV
                )
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict[Name]["Fitrange"],
                        0,
                        xmin=RMCPdict[Name]["Datarange"][0],
                        xmax=3.0,
                    ),
                    0,
                    WACV,
                )
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict[Name]["Fitrange"],
                        1,
                        xmin=10.0,
                        xmax=RMCPdict[Name]["Datarange"][1],
                    ),
                    0,
                    WACV,
                )
                fileSizer.Add(
                    wx.StaticText(G2frame.FRMC, label=" Scale factor: "), 0, WACV
                )
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, RMCPdict[Name]["dscale"], 0, xmin=0.001, xmax=20.0
                    ),
                    0,
                    WACV,
                )
                scaleref = wx.CheckBox(G2frame.FRMC, label="refine")
                scaleref.SetValue(RMCPdict[Name]["dscale"][1])
                Indx[scaleref.GetId()] = [Name, "dscale"]
                scaleref.Bind(wx.EVT_CHECKBOX, OnRef)
                fileSizer.Add(scaleref, 0, WACV)
                fileSizer.Add(wx.StaticText(G2frame.FRMC, label=" Qdamp "), 0, WACV)
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, RMCPdict[Name]["qdamp"], 0, xmin=0.0, xmax=1.0
                    ),
                    0,
                    WACV,
                )
                qdampref = wx.CheckBox(G2frame.FRMC, label="refine")
                qdampref.SetValue(RMCPdict[Name]["qdamp"][1])
                Indx[qdampref.GetId()] = [Name, "qdamp"]
                qdampref.Bind(wx.EVT_CHECKBOX, OnRef)
                fileSizer.Add(qdampref, 0, WACV)
                fileSizer.Add(wx.StaticText(G2frame.FRMC, label=" Qbroad "), 0, WACV)
                fileSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, RMCPdict[Name]["qbroad"], 0, xmin=0.0, xmax=1.0
                    ),
                    0,
                    WACV,
                )
                qbroadref = wx.CheckBox(G2frame.FRMC, label="refine")
                qbroadref.SetValue(RMCPdict[Name]["qbroad"][1])
                Indx[qbroadref.GetId()] = [Name, "qbroad"]
                qbroadref.Bind(wx.EVT_CHECKBOX, OnRef)
                fileSizer.Add(qbroadref, 0, WACV)

        mainSizer.Add(fileSizer, 0)

        return mainSizer

    def fullrmcSizer(RMCPdict):
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(
            wx.StaticText(
                G2frame.FRMC,
                label="""* "Atomic Stochastic Modeling & Optimization with fullrmc", B. Aoun, J. Appl. Cryst. 2022, 55(6) 1664-1676,
 DOI: 10.1107/S1600576722008536;
* "Fullrmc, a Rigid Body Reverse Monte Carlo Modeling Package Enabled with Machine Learning and Artificial
   Intelligence", B. Aoun, Jour. Comp. Chem. (2016), 37, 1102-1111. DOI: 10.1002/jcc.24304;
* www.fullrmc.com
 """,
            )
        )
        # if G2pwd.findfullrmc() is None:
        #     mainSizer.Add(wx.StaticText(G2frame.FRMC,
        #         label="\nsorry, fullrmc not installed or was not located"))
        #     return mainSizer
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC, True)
        # mainSizer.Add(wx.StaticText(G2frame.FRMC,label=' fullrmc big box starting pdb file preparation:'),0)

        # initialize fullrmc dictionary if needed
        RMCPdict = data["RMC"]["fullrmc"] = data["RMC"].get("fullrmc", {})
        # update, if atoms list has been updated
        Atypes = [
            atype.split("+")[0].split("-")[0] for atype in data["General"]["AtomTypes"]
        ]
        aTypes = dict(
            zip(
                Atypes,
                len(Atypes)
                * [
                    0.10,
                ],
                strict=False,
            )
        )
        if len(data["RMC"]["fullrmc"].get("aTypes", {})) != len(aTypes):
            # print('atypes has changed')
            atSeq = list(aTypes.keys())
            lenA = len(atSeq)
            Pairs = []
            for pair in [
                [
                    " %s-%s" % (atSeq[i], atSeq[j])
                    for j in range(i, lenA)
                    if "Va" not in atSeq[j]
                ]
                for i in range(lenA)
                if "Va" not in atSeq[i]
            ]:
                Pairs += pair
            Pairs = {pairs: [0.0, 0.0, 0.0] for pairs in Pairs}
            RMCPdict.update({"aTypes": aTypes, "atSeq": atSeq, "Pairs": Pairs})
        RMCPdict["files"] = RMCPdict.get(
            "files",
            {
                "Neutron real space data; G(r): ": ["Select", 1.0, "G(r)", 0, True],
                "Neutron reciprocal space data; S(Q)-1: ": [
                    "Select",
                    1.0,
                    "F(Q)",
                    0,
                    True,
                ],
                "Xray real space data; G(r): ": ["Select", 1.0, "G(r)", 0, True],
                "Xray reciprocal space data; S(Q)-1: ": [
                    "Select",
                    1.0,
                    "F(Q)",
                    0,
                    True,
                ],
            },
        )
        if "moleculePdb" not in RMCPdict:
            RMCPdict.update(
                {"moleculePdb": "Select", "targetDensity": 1.0, "maxRecursion": 10000}
            )
        if "Angles" not in RMCPdict:
            RMCPdict.update(
                {
                    "Angles": [],
                    "Angle Weight": 1.0e-5,
                    "Bond Weight": 1.0e-5,
                    "Torsions": [],
                    "Torsion Weight": 1.0e-5,
                }
            )
        for key, val in {
            "SuperCell": [1, 1, 1],
            "Box": [10.0, 10.0, 10.0],
            "ReStart": [False, False],
            "Cycles": 1,
            "Swaps": [],
            "useBVS": False,
            "FitScale": False,
            "AveCN": [],
            "FxCN": [],
            "min Contact": 1.5,
            "periodicBound": True,
        }.items():
            RMCPdict[key] = RMCPdict.get(key, val)

        def GetSuperSizer():
            def ShowRmax(*args, **kwargs):
                cell = data["General"]["Cell"][1:7]
                bigcell = np.array(cell) * np.array(RMCPdict["SuperCell"] + [1, 1, 1])
                bigG = G2lat.cell2Gmat(bigcell)[0]
                rmax = min(
                    [0.5 / np.sqrt(G2lat.calc_rDsq2(H, bigG)) for H in np.eye(3)]
                )
                rmaxlbl.SetLabel(f"  Rmax = {rmax:.1f}")

            superSizer = wx.BoxSizer(wx.HORIZONTAL)
            axes = ["X", "Y", "Z"]
            for i, ax in enumerate(axes):
                superSizer.Add(
                    wx.StaticText(G2frame.FRMC, label=" %s-axis: " % ax), 0, WACV
                )
                superSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict["SuperCell"],
                        i,
                        xmin=1,
                        xmax=20,
                        size=(50, 25),
                        OnLeave=ShowRmax,
                    ),
                    0,
                    WACV,
                )
            rmaxlbl = wx.StaticText(G2frame.FRMC, label=" Rmax=?")
            superSizer.Add(rmaxlbl, 0, WACV)
            ShowRmax()
            return superSizer

        def GetBoxSizer():
            boxSizer = wx.BoxSizer(wx.HORIZONTAL)
            axes = ["X", "Y", "Z"]
            for i, ax in enumerate(axes):
                boxSizer.Add(
                    wx.StaticText(G2frame.FRMC, label=" %s-axis: " % ax), 0, WACV
                )
                boxSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict["Box"],
                        i,
                        xmin=10.0,
                        xmax=50.0,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            return boxSizer

        def OnReStart(event):
            RMCPdict["ReStart"][0] = not RMCPdict["ReStart"][0]

        def OnAddSwap(event):
            RMCPdict["Swaps"].append(
                [
                    "",
                    "",
                    0.0,
                ]
            )
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnPdbButton(event):
            dlg = wx.FileDialog(
                G2frame.FRMC,
                "Choose molecule pdb file",
                G2frame.LastGPXdir,
                style=wx.FD_OPEN,
                wildcard="PDB file(*.pdb)|*.pdb",
            )
            if dlg.ShowModal() == wx.ID_OK:
                fpath, fName = os.path.split(dlg.GetPath())
                RMCPdict["moleculePdb"] = fName
                pdbButton.SetLabel(fName)

        def OnAddAngle(event):
            RMCPdict["Angles"].append(["", "", "", 0.0, 0.0, 0.0, 0.0])
            wx.CallAfter(UpdateRMC, G2frame, data)

        # def OnAddTorsion(event):
        #     RMCPdict['Torsions'].append(['','','','',0.,0.,0.,0.,0.,0.])
        #     wx.CallAfter(,G2frame,data,event)

        def GetAngleSizer():
            def OnDelAngle(event):
                Obj = event.GetEventObject()
                angle = Indx[Obj.GetId()]
                del RMCPdict["Angles"][angle]
                wx.CallAfter(UpdateRMC, G2frame, data)

            # def OnAngleAtSel(event):
            #     Obj = event.GetEventObject()
            #     angle,i = Indx[Obj.GetId()]
            #     RMCPdict['Angles'][angle][i] = Obj.GetStringSelection()

            def SetRestart1(invalid, value, tc):
                RMCPdict["ReStart"][1] = True

            Indx = {}
            atChoice = [atm for atm in RMCPdict["atSeq"] if "Va" not in atm]
            angleSizer = wx.GridBagSizer(0, 5)
            fxcnLabels1 = [
                " ",
                " ",
                "Central",
                "",
                None,
                "angle restraint values (deg)",
                None,
                "search distance (A)",
            ]
            fxcnLabels2 = [" ", "Atom-A", "Atom", "Atom-C", "min", "max", "from", "to"]
            for i in range(8):
                if fxcnLabels1[i]:
                    cspan = 1
                    coloff = 0
                    if fxcnLabels1[i - 1] is None:
                        cspan = 2
                        coloff = 1
                    angleSizer.Add(
                        wx.StaticText(G2frame.FRMC, label=fxcnLabels1[i]),
                        (0, i - coloff),
                        (1, cspan),
                    )
                if fxcnLabels2[i]:
                    angleSizer.Add(
                        wx.StaticText(
                            G2frame.FRMC,
                            wx.ID_ANY,
                            label=fxcnLabels2[i],
                            style=wx.CENTER,
                        ),
                        (1, i),
                    )
            row = 1
            for ifx, angle in enumerate(RMCPdict["Angles"]):
                row += 1
                angleSizer.Add((30, -1), (row, 0))
                for i in range(3):
                    if angle[i] not in atChoice:
                        angle[i] = atChoice[0]
                    atmSel = G2G.EnumSelector(G2frame.FRMC, angle, i, atChoice)
                    angleSizer.Add(atmSel, (row, 1 + i))
                for i in range(4):
                    if i == 0:
                        xmin, xmax = 0.0, 180.0
                    elif i == 2:
                        xmin, xmax = 0.1, 6.0
                    angleSizer.Add(
                        G2G.ValidatedTxtCtrl(
                            G2frame.FRMC,
                            angle,
                            3 + i,
                            xmin=xmin,
                            xmax=xmax,
                            OnLeave=SetRestart1,
                            size=(50, 25),
                        ),
                        (row, 4 + i),
                    )
                delBtn = wx.Button(G2frame.FRMC, label="Del", style=wx.BU_EXACTFIT)
                delBtn.Bind(wx.EVT_BUTTON, OnDelAngle)
                Indx[delBtn.GetId()] = ifx
                angleSizer.Add(delBtn, (row, 9))
            return angleSizer

        # def GetTorsionSizer():

        #     def OnDelTorsion(event):
        #         Obj = event.GetEventObject()
        #         angle = Indx[Obj.GetId()]
        #         del RMCPdict['Torsions'][angle]
        #         wx.CallAfter(UpdateRMC,G2frame,data)

        #     def OnTorsionAtSel(event):
        #         Obj = event.GetEventObject()
        #         torsion,i = Indx[Obj.GetId()]
        #         RMCPdict['Torsions'][torsion][i] = Obj.GetStringSelection()

        #     def SetRestart1(invalid,value,tc):
        #         RMCPdict['ReStart'][1] = True

        #     Indx = {}
        #     atChoice = [atm for atm in RMCPdict['atSeq'] if 'Va' not in atm]
        #     torsionSizer = wx.FlexGridSizer(11,5,5)
        #     fxcnLabels = [' ','Atom-A','Atom-B','Atom-C','Atom-D',' min angle1',' max angle1',' min angle2',' max angle2',' min angle3',' max angle3']
        #     for lab in fxcnLabels:
        #         torsionSizer.Add(wx.StaticText(G2frame.FRMC,label=lab),0,WACV)
        #     for ifx,torsion in enumerate(RMCPdict['Torsions']):
        #         delBtn = wx.Button(G2frame.FRMC,label='Delete')
        #         delBtn.Bind(wx.EVT_BUTTON,OnDelTorsion)
        #         Indx[delBtn.GetId()] = ifx
        #         torsionSizer.Add(delBtn,0,WACV)
        #         for i in [0,1,2,3]:
        #             atmSel = wx.ComboBox(G2frame.FRMC,choices=atChoice,style=wx.CB_DROPDOWN|wx.TE_READONLY)
        #             atmSel.SetStringSelection(torsion[i])
        #             atmSel.Bind(wx.EVT_COMBOBOX,OnTorsionAtSel)
        #             Indx[atmSel.GetId()] = [ifx,i]
        #             torsionSizer.Add(atmSel,0,WACV)
        #         for i in  [4,5,6,7,8,9]:
        #             torsionSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,torsion,i,xmin=0.,xmax=360.,OnLeave=SetRestart1,size=(50,25)),0,WACV)
        #     return torsionSizer

        generalData = data["General"]
        cx, ct, cs, cia = generalData["AtomPtrs"]
        # atomData = data['Atoms']
        # atNames = [atom[ct-1] for atom in atomData]
        # ifP1 = False
        # if generalData['SGData']['SpGrp'] == 'P 1':
        #     ifP1 = True
        ifBox = False
        if "macromolecular" in generalData["Type"]:
            ifBox = True
        lineSizer = wx.BoxSizer(wx.HORIZONTAL)
        if ifBox:
            lineSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" Big box dimensions, %s:" % Angstr),
                0,
                WACV,
            )
            lineSizer.Add(GetBoxSizer(), 0, WACV)
        #            elif ifP1:
        lineSizer.Add(
            wx.StaticText(G2frame.FRMC, label=" Lattice multipliers:"), 0, WACV
        )
        lineSizer.Add(GetSuperSizer(), 0, WACV)
        lineSizer.Add((5, -1))
        # Bachir suggests that w/o periodic boundaries, users are likely to use fullrmc wrong
        # lineSizer.Add(G2G.G2CheckBox(G2frame.FRMC,'Impose periodic boundaries',RMCPdict,'periodicBound'),
        #                  0,WACV)
        mainSizer.Add(lineSizer, 0)
        if ifBox:
            molecSizer = wx.BoxSizer(wx.HORIZONTAL)
            molecSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" Source molecule file "), 0, WACV
            )
            pdbButton = wx.Button(G2frame.FRMC, label=RMCPdict["moleculePdb"])
            pdbButton.Bind(wx.EVT_BUTTON, OnPdbButton)
            molecSizer.Add(pdbButton, 0, WACV)
            molecSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" target density, gm/cc "), 0, WACV
            )
            molecSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC, RMCPdict, "targetDensity", xmin=0.1, size=[60, 25]
                ),
                0,
                WACV,
            )
            molecSizer.Add(wx.StaticText(G2frame.FRMC, label=" max tries "), 0, WACV)
            molecSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC,
                    RMCPdict,
                    "maxRecursion",
                    xmin=1000,
                    xmax=1000000,
                    size=[60, 25],
                ),
                0,
                WACV,
            )
            mainSizer.Add(molecSizer, 0)
        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(
            wx.StaticText(G2frame.FRMC, label=" fullrmc run file preparation:")
        )
        resLine = wx.BoxSizer(wx.HORIZONTAL)
        resLine.Add(wx.StaticText(G2frame.FRMC, label=" Run "), 0, WACV)
        resLine.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC, RMCPdict, "Cycles", xmin=1, size=[60, 25]
            )
        )
        resLine.Add(
            wx.StaticText(G2frame.FRMC, label=" computation cycles of "), 0, WACV
        )
        RMCPdict["Steps/cycle"] = RMCPdict.get("Steps/cycle", 5000)
        resLine.Add(
            G2G.EnumSelector(
                G2frame.FRMC,
                RMCPdict,
                "Steps/cycle",
                ["1K", "5K", "10K", "50K"],
                [1000, 5000, 10000, 50000],
            ),
            0,
            WACV,
        )
        resLine.Add(wx.StaticText(G2frame.FRMC, label=" steps per cycle"), 0, WACV)
        mainSizer.Add(resLine, 0)
        resLine = wx.BoxSizer(wx.HORIZONTAL)
        resLine.Add(
            wx.StaticText(G2frame.FRMC, label=" Restart fullrmc Engine? "), 0, WACV
        )
        restart = wx.CheckBox(G2frame.FRMC, label="(will clear old result!) ")
        resLine.Add(restart, 0, WACV)

        restart.SetValue(RMCPdict["ReStart"][0])
        restart.Bind(wx.EVT_CHECKBOX, OnReStart)
        mainSizer.Add(resLine, 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(GetAtmChoice(G2frame.FRMC, RMCPdict), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        swapBox = wx.BoxSizer(wx.HORIZONTAL)
        swapBox.Add(
            wx.StaticText(G2frame.FRMC, label="Atom swap probabiities: "), 0, WACV
        )
        swapAdd = wx.Button(G2frame.FRMC, label="Add", style=wx.BU_EXACTFIT)
        swapAdd.Bind(wx.EVT_BUTTON, OnAddSwap)
        swapBox.Add(swapAdd, 0, WACV)
        mainSizer.Add(swapBox, 0)
        if len(RMCPdict["Swaps"]):
            mainSizer.Add(GetSwapSizer(RMCPdict), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(
            wx.StaticText(G2frame.FRMC, label="Geometry constraints && restraints"), 0
        )
        distBox = wx.BoxSizer(wx.HORIZONTAL)
        distBox.Add(wx.StaticText(G2frame.FRMC, label="Distance constraints"), 0, WACV)
        # weights removed for now
        # distBox.Add(wx.StaticText(G2frame.FRMC,label=', distance weight:'),0,WACV)
        # distBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Bond Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
        distBox.Add(wx.StaticText(G2frame.FRMC, label=" min contact dist: "), 0, WACV)
        distBox.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC, RMCPdict, "min Contact", xmin=0.0, xmax=4.0, size=(50, 25)
            ),
            0,
            WACV,
        )

        RMCPdict["useBondConstraints"] = RMCPdict.get("useBondConstraints", True)
        distBox.Add(
            wx.StaticText(G2frame.FRMC, label="  Use bond constraints? "), 0, WACV
        )
        distBox.Add(
            G2G.G2CheckBox(
                G2frame.FRMC,
                "",
                RMCPdict,
                "useBondConstraints",
                OnChange=lambda event: UpdateRMC(G2frame, data),
            ),
            0,
            WACV,
        )
        mainSizer.Add(distBox, 0)

        if RMCPdict["useBondConstraints"]:
            mainSizer.Add(GetPairSizer(G2frame.FRMC, RMCPdict), 0)
            mainSizer.Add((-1, 10))
        angBox = wx.BoxSizer(wx.HORIZONTAL)
        angBox.Add(wx.StaticText(G2frame.FRMC, label="A-B-C angle restraints"), 0, WACV)
        # weights removed for now
        # angBox.Add(wx.StaticText(G2frame.FRMC,label=', angle weight:'),0,WACV)
        # angBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Angle Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
        angBox.Add((20, -1))
        angAdd = wx.Button(G2frame.FRMC, label="Add", style=wx.BU_EXACTFIT)
        angAdd.Bind(wx.EVT_BUTTON, OnAddAngle)
        angBox.Add(angAdd, 0, WACV)
        mainSizer.Add(angBox, 0)
        if len(RMCPdict["Angles"]):
            mainSizer.Add(GetAngleSizer(), 0)
        RMCPdict["Groups"] = RMCPdict.get("Groups", [])

        def OnAddGroup(event):
            index = len(RMCPdict["Groups"])
            RMCPdict["Groups"].append([])
            GroupEditor(index)

        def OnDelGroup(event):
            index = event.EventObject.index
            del RMCPdict["Groups"][index]
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnEdtGroup(event):
            index = event.EventObject.index
            GroupEditor(index)

        def GroupEditor(index):
            cx, ct, cs, cia = data["General"]["AtomPtrs"]
            atomlbs = [a[ct - 1] for a in data["Atoms"]]
            dlg = G2G.G2MultiChoiceDialog(
                G2frame.FRMC,
                "Atom Selector",
                "Select atoms to include in group",
                atomlbs,
                selected=RMCPdict["Groups"][index],
            )
            if dlg.ShowModal() == wx.ID_OK:
                RMCPdict["Groups"][index] = dlg.GetSelections()
            dlg.Destroy()
            if len(RMCPdict["Groups"][index]) == 0:
                del RMCPdict["Groups"][index]
            wx.CallAfter(UpdateRMC, G2frame, data)

        if len(RMCPdict["Groups"]) == 0:
            grpAdd = wx.Button(
                G2frame.FRMC, label="Define atom group", style=wx.BU_EXACTFIT
            )
            grpAdd.Bind(wx.EVT_BUTTON, OnAddGroup)
            mainSizer.Add(grpAdd, 0)
        else:
            grpBox = wx.BoxSizer(wx.HORIZONTAL)
            grpBox.Add(wx.StaticText(G2frame.FRMC, label="Atom Groups:  "), 0, WACV)
            grpAdd = wx.Button(G2frame.FRMC, label="Add group", style=wx.BU_EXACTFIT)
            grpAdd.Bind(wx.EVT_BUTTON, OnAddGroup)
            RMCPdict["GroupMode"] = RMCPdict.get("GroupMode", 0)
            grpBox.Add(grpAdd, 0, WACV)
            grpBox.Add(
                wx.StaticText(G2frame.FRMC, label="  Group refinement mode: "), 0, WACV
            )
            grpBox.Add(
                G2G.EnumSelector(
                    G2frame.FRMC,
                    RMCPdict,
                    "GroupMode",
                    ("Rotate & Translate", "Rotate only", "Translate only"),
                    [0, 1, 2],
                ),
                0,
                WACV,
            )
            mainSizer.Add(grpBox, 0)
            for i, g in enumerate(RMCPdict["Groups"]):
                grpBox = wx.BoxSizer(wx.HORIZONTAL)
                grpBox.Add((20, -1))
                grpBox.Add(
                    wx.StaticText(G2frame.FRMC, label="Group #" + str(i + 1)), 0, WACV
                )
                grpBox.Add((4, -1))
                grpdel = wx.Button(G2frame.FRMC, label="Del", style=wx.BU_EXACTFIT)
                grpdel.Bind(wx.EVT_BUTTON, OnDelGroup)
                grpdel.index = i
                grpBox.Add(grpdel, 0, WACV)
                grpadd = wx.Button(G2frame.FRMC, label="Edit", style=wx.BU_EXACTFIT)
                grpadd.Bind(wx.EVT_BUTTON, OnEdtGroup)
                grpadd.index = i
                grpBox.Add(grpadd, 0, WACV)
                msg = " Contains atoms: "
                for i, n in enumerate(g):
                    if i + 1 == len(g):
                        msg += " && "
                    elif i > 0:
                        msg += ", "
                    msg += str(i)
                grpBox.Add(wx.StaticText(G2frame.FRMC, label=msg), 0, WACV)
                mainSizer.Add(grpBox, 0)

        RMCPdict["addThermalBroadening"] = RMCPdict.get("addThermalBroadening", False)
        mainSizer.Add((-1, 5))
        distBox = wx.BoxSizer(wx.HORIZONTAL)
        distBox.Add(
            wx.StaticText(G2frame.FRMC, label=" Add thermal broadening? "), 0, WACV
        )
        distBox.Add(
            G2G.G2CheckBox(
                G2frame.FRMC,
                "",
                RMCPdict,
                "addThermalBroadening",
                OnChange=lambda event: UpdateRMC(G2frame, data),
            ),
            0,
            WACV,
        )
        if RMCPdict["addThermalBroadening"]:
            distBox.Add((15, -1))
            distBox.Add(wx.StaticText(G2frame.FRMC, label="Uiso equiv."), 0, WACV)
            RMCPdict["ThermalU"] = RMCPdict.get("ThermalU", {})
            for atm in RMCPdict["aTypes"]:
                RMCPdict["ThermalU"][atm] = RMCPdict["ThermalU"].get(atm, 0.005)
                distBox.Add(
                    wx.StaticText(G2frame.FRMC, label="  " + atm + ":"), 0, WACV
                )
                distBox.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict["ThermalU"],
                        atm,
                        xmin=0.0001,
                        xmax=0.25,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
        mainSizer.Add(distBox, 0)
        if RMCPdict["addThermalBroadening"]:
            mainSizer.Add((-1, 5))

        # Torsions are difficult to implement. Need to be internal to a unit cell & named with fullrmc
        # atom labels. Leave this out, at least for now.
        # torBox = wx.BoxSizer(wx.HORIZONTAL)
        # torAdd = wx.Button(G2frame.FRMC,label='Add')
        # torAdd.Bind(wx.EVT_BUTTON,OnAddTorsion)
        # torBox.Add(torAdd,0,WACV)
        # torBox.Add(wx.StaticText(G2frame.FRMC,label=' A-B-C-D torsion angle restraints (intracell only), weight: '),0,WACV)
        # torBox.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict,'Torsion Weight',xmin=0.,xmax=100.,size=(50,25)),0,WACV)
        # mainSizer.Add(torBox,0)
        # if len(RMCPdict['Torsions']):
        #     mainSizer.Add(GetTorsionSizer(),0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(FileSizer(RMCPdict))
        return mainSizer

    def RMCProfileSizer(RMCPdict):
        def CheckAtms(Atypes):
            newAtm = False
            for atm in Atypes:
                if atm not in data["RMC"]["RMCProfile"].get("aTypes", {}):
                    newAtm = True
                    break
            for atm in data["RMC"]["RMCProfile"].get("aTypes", {}):
                if atm not in Atypes:
                    newAtm = True
                    break
            return newAtm

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        subSizer = wx.BoxSizer(wx.HORIZONTAL)
        subSizer.Add((-1, -1), 1, wx.EXPAND)
        subSizer.Add(wx.StaticText(G2frame.FRMC, label="RMCProfile setup"), 0, WACV)
        subSizer.Add((-1, -1), 1, wx.EXPAND)
        mainSizer.Add(subSizer)
        mainSizer.Add((5, 5))
        txt = wx.StaticText(
            G2frame.FRMC, label=f"Please cite: {G2G.GetCite('RMCProfile')}"
        )
        txt.Wrap(500)
        mainSizer.Add(txt)
        mainSizer.Add((5, 5))
        Atypes = [
            atype.split("+")[0].split("-")[0] for atype in data["General"]["AtomTypes"]
        ]
        aTypes = dict(
            zip(
                Atypes,
                len(Atypes)
                * [
                    0.10,
                ],
                strict=False,
            )
        )
        atSeq = list(aTypes.keys())
        lenA = len(atSeq)
        atOxid = [[atmdata.BVSoxid[atm][0], 0.001] for atm in atSeq]
        if CheckAtms(Atypes):
            oldPairs = data["RMC"]["RMCProfile"].get("Pairs", {})
            Pairs = {}
            #                for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
            for pairs in [
                [" %s-%s" % (atSeq[i], atSeq[j]) for j in range(i, lenA)]
                for i in range(lenA)
            ]:
                for pair in pairs:
                    if pair in oldPairs:
                        Pairs[pair] = oldPairs[pair]
                    else:
                        Pairs[pair] = [0.0, 0.0, 0.0]
            data["RMC"]["RMCProfile"].update(
                {
                    "aTypes": aTypes,
                    "atSeq": atSeq,
                    "Pairs": Pairs,
                    "Oxid": atOxid,
                }
            )

        if not data["RMC"]["RMCProfile"] or "metadata" not in RMCPdict:
            Pairs = {}
            #                for pairs in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
            for pairs in [
                [" %s-%s" % (atSeq[i], atSeq[j]) for j in range(i, lenA)]
                for i in range(lenA)
            ]:
                for pair in pairs:
                    Pairs[pair] = [0.0, 0.0, 0.0]
            BVSpairs = []
            if lenA > 1:
                #                    for pair in [[' %s-%s'%(atSeq[i],atSeq[j]) for j in range(i,lenA) if 'Va' not in atSeq[j]] for i in range(lenA) if 'Va' not in atSeq[i]]:
                for pair in [
                    [" %s-%s" % (atSeq[i], atSeq[j]) for j in range(i, lenA)]
                    for i in range(lenA)
                ]:
                    BVSpairs += pair
            BVS = {pairs: [0.0, 0.0, 0.0, 0.0] for pairs in BVSpairs}
            files = {
                "Neutron real space data; G(r): ": [
                    "Select",
                    0.05,
                    "G(r)",
                    "RMC",
                ],
                "Neutron reciprocal space data; F(Q): ": [
                    "Select",
                    0.05,
                    "F(Q)",
                    "RMC",
                ],
                "Neutron reciprocal space data; S(Q): ": [
                    "Select",
                    0.05,
                    "S(Q)",
                    "RMC",
                ],
                #                          'Xray real space data; G(r): ':['Select',0.01,'G(r)','RMC',],
                "Xray reciprocal space data; F(Q): ": [
                    "Select",
                    0.01,
                    "F(Q)",
                    "RMC",
                ],
            }
            runTimes = [10.0, 1.0]
            metadata = {
                "title": "none",
                "owner": "no one",
                "date": str(time.ctime()),
                "temperature": "300K",
                "material": "nothing",
                "phase": "vacuum",
                "comment": "none ",
                "source": "nowhere",
            }
            data["RMC"]["RMCProfile"].update(
                {
                    "SuperCell": [1, 1, 1],
                    "UseSampBrd": [True, True],
                    "aTypes": aTypes,
                    "histogram": ["", 1.0],
                    "files": files,
                    "metadata": metadata,
                    "FitScale": False,
                    "atSeq": atSeq,
                    "runTimes": runTimes,
                    "ReStart": [False, False],
                    "BVS": BVS,
                    "Oxid": atOxid,
                    "useBVS": False,
                    "Swaps": [],
                    "AveCN": [],
                    "FxCN": [],
                    "Potentials": {
                        "Angles": [],
                        "Angle search": 10.0,
                        "Stretch": [],
                        "Pairs": Pairs,
                        "Stretch search": 10.0,
                        "Pot. Temp.": 300.0,
                        "useGPU": False,
                    },
                }
            )

        #            data['RMC']['RMCProfile']['aTypes'] = {aTypes[atype] for atype in aTypes if atype in Atypes}
        data["RMC"]["RMCProfile"]["Isotope"] = copy.copy(data["General"]["Isotope"])
        data["RMC"]["RMCProfile"]["Isotopes"] = copy.deepcopy(
            data["General"]["Isotopes"]
        )
        data["RMC"]["RMCProfile"]["NoAtoms"] = copy.copy(data["General"]["NoAtoms"])
        RMCPdict = data["RMC"]["RMCProfile"]
        # patches
        if "FitScale" not in RMCPdict:
            RMCPdict["FitScale"] = False
        if "useGPU" not in RMCPdict:
            RMCPdict["useGPU"] = False

        # end patches

        def OnHisto(event):
            RMCPdict["histogram"][0] = histo.GetStringSelection()

        def OnSize(event):
            RMCPdict["UseSampBrd"][0] = samSize.GetValue()

        def OnStrain(event):
            RMCPdict["UseSampBrd"][1] = strain.GetValue()

        def OnFitScale(event):
            RMCPdict["FitScale"] = not RMCPdict["FitScale"]

        def SetRestart(invalid, value, tc):
            RMCPdict["ReStart"] = [True, True]

        def OnUseBVS(event):
            RMCPdict["useBVS"] = not RMCPdict["useBVS"]
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnAddSwap(event):
            RMCPdict["Swaps"].append(
                [
                    "",
                    "",
                    0.0,
                ]
            )
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnAddFxCN(event):
            RMCPdict["FxCN"].append(["", "", 0.5, 2.0, 6, 1.0, 0.00001])
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnAddAveCN(event):
            RMCPdict["AveCN"].append(["", "", 0.5, 2.0, 6.0, 0.00001])
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnAddAnglePot(event):
            RMCPdict["Potentials"]["Angles"].append(["", "", "", 0.0, 0.0, 0.0, 0.0])
            wx.CallAfter(UpdateRMC, G2frame, data)

        def OnAddBondPot(event):
            RMCPdict["Potentials"]["Stretch"].append(["", "", 0.0, 0.0])
            wx.CallAfter(UpdateRMC, G2frame, data)

        def GetTimeSizer():
            def OnUseGPU(event):
                RMCPdict["useGPU"] = not RMCPdict["useGPU"]

            timeSizer = wx.BoxSizer(wx.HORIZONTAL)
            timeSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" Total running time (min): "),
                0,
                WACV,
            )
            timeSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC, RMCPdict["runTimes"], 0, xmin=0.0, size=(70, 25)
                ),
                0,
                WACV,
            )
            timeSizer.Add(
                wx.StaticText(G2frame.FRMC, label=" Save interval time (min): "),
                0,
                WACV,
            )
            timeSizer.Add(
                G2G.ValidatedTxtCtrl(
                    G2frame.FRMC,
                    RMCPdict["runTimes"],
                    1,
                    xmin=0.1,
                    xmax=20.0,
                    size=(50, 25),
                ),
                0,
                WACV,
            )
            usegpu = wx.CheckBox(G2frame.FRMC, label=" use GPU?")
            usegpu.SetValue(RMCPdict["useGPU"])
            usegpu.Bind(wx.EVT_CHECKBOX, OnUseGPU)
            timeSizer.Add(usegpu, 0, WACV)
            return timeSizer

        # def GetSuperSizer(Xmax):
        #     superSizer = wx.BoxSizer(wx.HORIZONTAL)
        #     axes = ['X','Y','Z']
        #     for i,ax in enumerate(axes):
        #         superSizer.Add(wx.StaticText(G2frame.FRMC,label=' %s-axis: '%ax),0,WACV)
        #         superSizer.Add(G2G.ValidatedTxtCtrl(G2frame.FRMC,RMCPdict['SuperCell'],
        #             i,xmin=1,xmax=xamx,size=(50,25),OnLeave=SetRestart),0,WACV)
        #     return superSizer

        def GetBvsSizer(pnl):
            def OnResetBVS(event):
                Obj = event.GetEventObject()
                pair = Indx[Obj.GetId()]
                pId = [key for key in RMCPdict["BVS"]].index(pair) + 1
                nId = len(RMCPdict["BVS"]) + 1
                dist = G2elem.GetBVS(pair, RMCPdict["atSeq"], RMCPdict["Oxid"])
                if dist:
                    RMCPdict["BVS"][pair] = [dist, 0.37, 3.0]
                    bvsCh = bvsSizer.GetChildren()
                    addr = 2 * nId + pId
                    bvsCh[addr].Window.SetValue("%6.3f" % dist)
                    bvsCh[addr + nId].Window.SetValue("0.37")
                    bvsCh[addr + 2 * nId].Window.SetValue("3.00")

            bvsSizer = wx.FlexGridSizer(len(RMCPdict["BVS"]) + 1, 5, 5)
            bvsSizer.Add((5, 5), 0)
            for pair in RMCPdict["BVS"]:
                bvsSizer.Add(wx.StaticText(pnl, label=pair), 0, WACV)
            bvsSizer.Add(wx.StaticText(pnl, label=" Reset:"), 0, WACV)
            for pair in RMCPdict["BVS"]:
                reset = wx.Button(pnl, label="Yes")
                bvsSizer.Add(reset, 0, WACV)
                reset.Bind(wx.EVT_BUTTON, OnResetBVS)
                Indx[reset.GetId()] = pair
            bvsSizer.Add(wx.StaticText(pnl, label=" Bond length:"), 0, WACV)
            for pair in RMCPdict["BVS"]:
                bvsSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        pnl,
                        RMCPdict["BVS"][pair],
                        0,
                        xmin=0.0,
                        xmax=10.0,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            bvsSizer.Add(wx.StaticText(pnl, label=" B constant (0.37): "), 0, WACV)
            for pair in RMCPdict["BVS"]:
                bvsSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        pnl,
                        RMCPdict["BVS"][pair],
                        1,
                        xmin=0.0,
                        xmax=10.0,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            bvsSizer.Add(wx.StaticText(pnl, label=" Cut off: "), 0, WACV)
            for pair in RMCPdict["BVS"]:
                bvsSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        pnl,
                        RMCPdict["BVS"][pair],
                        2,
                        xmin=0.0,
                        xmax=10.0,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            return bvsSizer

        def GetFxcnSizer():
            def OnDelFxCN(event):
                Obj = event.GetEventObject()
                fxCN = Indx[Obj.GetId()]
                del RMCPdict["FxCN"][fxCN]
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnFxcnAtSel(event):
                Obj = event.GetEventObject()
                ifxCN, i = Indx[Obj.GetId()]
                RMCPdict["FxCN"][ifxCN][i] = Obj.GetStringSelection()

            fxcnSizer = wx.FlexGridSizer(8, 5, 5)
            atChoice = [atm for atm in RMCPdict["atSeq"] if "Va" not in atm]
            fxcnLabels = [
                " ",
                "Atom-1",
                "Atom-2",
                "min dist",
                "max dist",
                "CN",
                "fraction",
                "weight",
            ]
            for lab in fxcnLabels:
                fxcnSizer.Add(wx.StaticText(G2frame.FRMC, label=lab), 0, WACV)
            for ifx, fxCN in enumerate(RMCPdict["FxCN"]):
                delBtn = wx.Button(G2frame.FRMC, label="Delete")
                delBtn.Bind(wx.EVT_BUTTON, OnDelFxCN)
                Indx[delBtn.GetId()] = ifx
                fxcnSizer.Add(delBtn, 0, WACV)
                for i in [0, 1]:
                    atmSel = wx.ComboBox(
                        G2frame.FRMC,
                        choices=atChoice,
                        style=wx.CB_DROPDOWN | wx.TE_READONLY,
                    )
                    atmSel.SetStringSelection(fxCN[i])
                    atmSel.Bind(wx.EVT_COMBOBOX, OnFxcnAtSel)
                    Indx[atmSel.GetId()] = [ifx, i]
                    fxcnSizer.Add(atmSel, 0, WACV)
                fxcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 2, xmin=0.0, xmax=5.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                fxcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 3, xmin=0.0, xmax=5.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                fxcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 4, xmin=1, xmax=12, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                fxcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 5, xmin=0.0, xmax=1.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                fxcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 6, xmin=0.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
            return fxcnSizer

        def GetAvcnSizer():
            def OnDelAvCN(event):
                Obj = event.GetEventObject()
                fxCN = Indx[Obj.GetId()]
                del RMCPdict["AveCN"][fxCN]
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnAvcnAtSel(event):
                Obj = event.GetEventObject()
                ifxCN, i = Indx[Obj.GetId()]
                RMCPdict["AveCN"][ifxCN][i] = Obj.GetStringSelection()

            avcnSizer = wx.FlexGridSizer(7, 5, 5)
            atChoice = [atm for atm in RMCPdict["atSeq"] if "Va" not in atm]
            fxcnLabels = [
                " ",
                "Atom-1",
                "Atom-2",
                "min dist",
                "max dist",
                "CN",
                "weight",
            ]
            for lab in fxcnLabels:
                avcnSizer.Add(wx.StaticText(G2frame.FRMC, label=lab), 0, WACV)
            for ifx, fxCN in enumerate(RMCPdict["AveCN"]):
                delBtn = wx.Button(G2frame.FRMC, label="Delete")
                delBtn.Bind(wx.EVT_BUTTON, OnDelAvCN)
                Indx[delBtn.GetId()] = ifx
                avcnSizer.Add(delBtn, 0, WACV)
                for i in [0, 1]:
                    atmSel = wx.ComboBox(
                        G2frame.FRMC,
                        choices=atChoice,
                        style=wx.CB_DROPDOWN | wx.TE_READONLY,
                    )
                    atmSel.SetStringSelection(fxCN[i])
                    atmSel.Bind(wx.EVT_COMBOBOX, OnAvcnAtSel)
                    Indx[atmSel.GetId()] = [ifx, i]
                    avcnSizer.Add(atmSel, 0, WACV)
                avcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 2, xmin=0.0, xmax=5.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                avcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 3, xmin=0.0, xmax=5.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                avcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 4, xmin=1.0, xmax=12.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
                avcnSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, fxCN, 5, xmin=0.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
            return avcnSizer

        def GetAngleSizer():
            def OnDelAngle(event):
                Obj = event.GetEventObject()
                angle = Indx[Obj.GetId()]
                del RMCPdict["Potentials"]["Angles"][angle]
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnAngleAtSel(event):
                Obj = event.GetEventObject()
                angle, i = Indx[Obj.GetId()]
                RMCPdict["Potentials"]["Angles"][angle][i] = Obj.GetStringSelection()

            def SetRestart1(invalid, value, tc):
                RMCPdict["ReStart"][1] = True

            atChoice = [atm for atm in RMCPdict["atSeq"] if "Va" not in atm]
            angleSizer = wx.FlexGridSizer(8, 5, 5)
            fxcnLabels = [
                " ",
                "Atom-A",
                "Atom-B",
                "Atom-C",
                " ABC angle",
                "AB dist",
                "BC dist",
                "potential",
            ]
            for lab in fxcnLabels:
                angleSizer.Add(wx.StaticText(G2frame.FRMC, label=lab), 0, WACV)
            for ifx, angle in enumerate(RMCPdict["Potentials"]["Angles"]):
                delBtn = wx.Button(G2frame.FRMC, label="Delete")
                delBtn.Bind(wx.EVT_BUTTON, OnDelAngle)
                Indx[delBtn.GetId()] = ifx
                angleSizer.Add(delBtn, 0, WACV)
                for i in [0, 1, 2]:
                    atmSel = wx.ComboBox(
                        G2frame.FRMC,
                        choices=atChoice,
                        style=wx.CB_DROPDOWN | wx.TE_READONLY,
                    )
                    atmSel.SetStringSelection(angle[i])
                    atmSel.Bind(wx.EVT_COMBOBOX, OnAngleAtSel)
                    Indx[atmSel.GetId()] = [ifx, i]
                    angleSizer.Add(atmSel, 0, WACV)
                angleSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        angle,
                        3,
                        xmin=0.0,
                        xmax=180.0,
                        OnLeave=SetRestart1,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
                angleSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        angle,
                        4,
                        xmin=0.5,
                        xmax=5.0,
                        OnLeave=SetRestart1,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
                angleSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        angle,
                        5,
                        xmin=0.5,
                        xmax=5.0,
                        OnLeave=SetRestart1,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
                angleSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        angle,
                        6,
                        xmin=0.0,
                        OnLeave=SetRestart1,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
            return angleSizer

        def GetBondSizer():
            def OnDelBond(event):
                Obj = event.GetEventObject()
                bond = Indx[Obj.GetId()]
                del RMCPdict["Potentials"]["Stretch"][bond]
                wx.CallAfter(UpdateRMC, G2frame, data)

            def OnBondAtSel(event):
                Obj = event.GetEventObject()
                bond, i = Indx[Obj.GetId()]
                RMCPdict["Potentials"]["Stretch"][bond][i] = Obj.GetStringSelection()

            def SetRestart1(invalid, value, tc):
                RMCPdict["ReStart"][1] = True

            atChoice = [atm for atm in RMCPdict["atSeq"] if "Va" not in atm]
            bondSizer = wx.FlexGridSizer(5, 5, 5)
            fxcnLabels = [" ", "Atom-A", "Atom-B", " AB dist", "potential"]
            for lab in fxcnLabels:
                bondSizer.Add(wx.StaticText(G2frame.FRMC, label=lab), 0, WACV)
            for ifx, bond in enumerate(RMCPdict["Potentials"]["Stretch"]):
                delBtn = wx.Button(G2frame.FRMC, label="Delete")
                delBtn.Bind(wx.EVT_BUTTON, OnDelBond)
                Indx[delBtn.GetId()] = ifx
                bondSizer.Add(delBtn, 0, WACV)
                for i in [0, 1]:
                    atmSel = wx.ComboBox(
                        G2frame.FRMC,
                        choices=atChoice,
                        style=wx.CB_DROPDOWN | wx.TE_READONLY,
                    )
                    atmSel.SetStringSelection(bond[i])
                    atmSel.Bind(wx.EVT_COMBOBOX, OnBondAtSel)
                    Indx[atmSel.GetId()] = [ifx, i]
                    bondSizer.Add(atmSel, 0, WACV)
                bondSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        bond,
                        2,
                        xmin=0.0,
                        xmax=5.0,
                        OnLeave=SetRestart1,
                        size=(50, 25),
                    ),
                    0,
                    WACV,
                )
                bondSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC, bond, 3, xmin=0.0, size=(50, 25)
                    ),
                    0,
                    WACV,
                )
            return bondSizer

        Indx = {}

        mainSizer.Add(wx.StaticText(G2frame.FRMC, label=" Enter metadata items:"), 0)
        mainSizer.Add(
            GetMetaSizer(
                RMCPdict,
                [
                    "title",
                    "owner",
                    "material",
                    "phase",
                    "comment",
                    "source",
                ],
            ),
            0,
        )

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(GetTimeSizer(), 0)

        mainSizer.Add(
            wx.StaticText(
                G2frame.FRMC,
                label=" Lattice multipliers; if changed will force reset of atom positions:",
            ),
            0,
        )
        mainSizer.Add(GetSuperSizer(RMCPdict, 20), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)

        mSizer = wx.BoxSizer(wx.VERTICAL)
        mSizer.Add(wx.StaticText(G2frame.FRMC, label="Enter atom settings"), 0)
        mSizer.Add(GetAtmChoice(G2frame.FRMC, RMCPdict), 0)
        mSizer.Add(
            wx.StaticText(
                G2frame.FRMC,
                label=" N.B.: be sure to set cations first && anions last in atom ordering",
            )
        )
        mainSizer.Add(mSizer)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        swapBox = wx.BoxSizer(wx.HORIZONTAL)
        swapAdd = wx.Button(G2frame.FRMC, label="Add")
        swapAdd.Bind(wx.EVT_BUTTON, OnAddSwap)
        swapBox.Add(swapAdd, 0, WACV)
        swapBox.Add(
            wx.StaticText(G2frame.FRMC, label=" Atom swap probabilities: "), 0, WACV
        )
        mainSizer.Add(swapBox, 0)
        if len(RMCPdict["Swaps"]):
            mainSizer.Add(GetSwapSizer(RMCPdict), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)

        mSizer = wx.BoxSizer(wx.VERTICAL)
        mSizer.Add(
            wx.StaticText(
                G2frame.FRMC,
                label="Enter constraints && restraints via minimum && maximum distances for atom pairs:",
            ),
            0,
        )
        mSizer.Add(GetPairSizer(G2frame.FRMC, RMCPdict), 0)
        mainSizer.Add(mSizer)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        useBVS = wx.CheckBox(
            G2frame.FRMC,
            label=" Use bond valence sum restraints for (set to 0 for non-bonded ones):",
        )
        useBVS.SetValue(RMCPdict.get("useBVS", False))
        useBVS.Bind(wx.EVT_CHECKBOX, OnUseBVS)
        mainSizer.Add(useBVS, 0)
        if RMCPdict.get("useBVS", False):
            mSizer = wx.BoxSizer(wx.VERTICAL)
            mSizer.Add(GetBvsSizer(G2frame.FRMC), 0)
            mainSizer.Add(mSizer)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        fxcnBox = wx.BoxSizer(wx.HORIZONTAL)
        fxcnAdd = wx.Button(G2frame.FRMC, label="Add")
        fxcnAdd.Bind(wx.EVT_BUTTON, OnAddFxCN)
        fxcnBox.Add(fxcnAdd, 0, WACV)
        fxcnBox.Add(
            wx.StaticText(G2frame.FRMC, label=" Fixed coordination number restraint: "),
            0,
            WACV,
        )
        mainSizer.Add(fxcnBox, 0)
        if len(RMCPdict["FxCN"]):
            mainSizer.Add(GetFxcnSizer(), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        avcnBox = wx.BoxSizer(wx.HORIZONTAL)
        avcnAdd = wx.Button(G2frame.FRMC, label="Add")
        avcnAdd.Bind(wx.EVT_BUTTON, OnAddAveCN)
        avcnBox.Add(avcnAdd, 0, WACV)
        avcnBox.Add(
            wx.StaticText(
                G2frame.FRMC, label=" Average coordination number restraint: "
            ),
            0,
            WACV,
        )
        mainSizer.Add(avcnBox, 0)
        if len(RMCPdict["AveCN"]):
            mainSizer.Add(GetAvcnSizer(), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        pottempBox = wx.BoxSizer(wx.HORIZONTAL)
        pottempBox.Add(
            wx.StaticText(G2frame.FRMC, label=" Potential temperature (K): "), 0, WACV
        )
        pottempBox.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC,
                RMCPdict["Potentials"],
                "Pot. Temp.",
                xmin=0.0,
                xmax=1000.0,
                size=(50, 25),
            ),
            0,
            WACV,
        )
        mainSizer.Add(pottempBox, 0)
        bondpotBox = wx.BoxSizer(wx.HORIZONTAL)
        bondpotAdd = wx.Button(G2frame.FRMC, label="Add")
        bondpotAdd.Bind(wx.EVT_BUTTON, OnAddBondPot)
        bondpotBox.Add(bondpotAdd, 0, WACV)
        bondpotBox.Add(
            wx.StaticText(
                G2frame.FRMC,
                label=" A-B stretch potential restraints, search range (%): ",
            ),
            0,
            WACV,
        )
        bondpotBox.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC,
                RMCPdict["Potentials"],
                "Stretch search",
                xmin=0.0,
                xmax=100.0,
                size=(50, 25),
            ),
            0,
            WACV,
        )
        mainSizer.Add(bondpotBox, 0)
        if len(RMCPdict["Potentials"]["Stretch"]):
            mainSizer.Add(GetBondSizer(), 0)

        angpotBox = wx.BoxSizer(wx.HORIZONTAL)
        angpotAdd = wx.Button(G2frame.FRMC, label="Add")
        angpotAdd.Bind(wx.EVT_BUTTON, OnAddAnglePot)
        angpotBox.Add(angpotAdd, 0, WACV)
        angpotBox.Add(
            wx.StaticText(
                G2frame.FRMC,
                label=" A-B-C angle potential restraints, search range (%): ",
            ),
            0,
            WACV,
        )
        angpotBox.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC,
                RMCPdict["Potentials"],
                "Angle search",
                xmin=0.0,
                xmax=100.0,
                size=(50, 25),
            ),
            0,
            WACV,
        )
        mainSizer.Add(angpotBox, 0)
        if len(RMCPdict["Potentials"]["Angles"]):
            mainSizer.Add(GetAngleSizer(), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(wx.StaticText(G2frame.FRMC, label=" Select data:"), 0)
        histograms = data["Histograms"]
        histNames = list(histograms.keys())
        mainSizer.Add(
            wx.StaticText(
                G2frame.FRMC, label=" Select one histogram for Bragg processing:"
            ),
            0,
        )
        histoSizer = wx.BoxSizer(wx.HORIZONTAL)
        histo = wx.ComboBox(
            G2frame.FRMC, choices=histNames, style=wx.CB_DROPDOWN | wx.TE_READONLY
        )
        if RMCPdict["histogram"][0] == "":
            RMCPdict["histogram"][0] = histo.GetStringSelection()
        histo.SetStringSelection(RMCPdict["histogram"][0])
        histo.Bind(wx.EVT_COMBOBOX, OnHisto)
        histoSizer.Add(histo, 0, WACV)
        histoSizer.Add(wx.StaticText(G2frame.FRMC, label=" Weight "), 0, WACV)
        histoSizer.Add(
            G2G.ValidatedTxtCtrl(
                G2frame.FRMC,
                RMCPdict["histogram"],
                1,
                xmin=0.0,
                xmax=10000.0,
                size=(50, 25),
            ),
            0,
            WACV,
        )
        mainSizer.Add(histoSizer, 0)

        samSizer = wx.BoxSizer(wx.HORIZONTAL)
        samSize = wx.CheckBox(G2frame.FRMC, label=" Use size broadening?")
        samSize.SetValue(RMCPdict["UseSampBrd"][0])
        samSize.Bind(wx.EVT_CHECKBOX, OnSize)
        strain = wx.CheckBox(G2frame.FRMC, label=" Use mustrain broadening?")
        strain.SetValue(RMCPdict["UseSampBrd"][1])
        strain.Bind(wx.EVT_CHECKBOX, OnStrain)
        samSizer.Add(samSize, 0, WACV)
        samSizer.Add(strain, 0, WACV)
        fitscale = wx.CheckBox(G2frame.FRMC, label=" Fit scale factors?")
        fitscale.SetValue(RMCPdict["FitScale"])
        fitscale.Bind(wx.EVT_CHECKBOX, OnFitScale)
        samSizer.Add(fitscale, 0, WACV)
        mainSizer.Add(samSizer, 0)

        mainSizer.Add(FileSizer(RMCPdict))
        return mainSizer

    def PDFfitSizer(data):
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        Indx = {}

        def PDFParmSizer():
            def OnShape(event):
                RMCPdict["shape"] = shape.GetValue()
                wx.CallAfter(UpdateRMC, G2frame, data)

            parmSizer = wx.FlexGridSizer(3, 6, 5, 5)
            Names = ["delta1", "delta2", "sratio", "rcut", "spdiameter"]
            Names2 = [
                "stepcut",
            ]
            for name in Names:

                def OnRefine(event):
                    Obj = event.GetEventObject()
                    name = Indx[Obj.GetId()]
                    RMCPdict[name][1] = not RMCPdict[name][1]

                if name == "spdiameter" and RMCPdict.get("shape", "sphere") != "sphere":
                    pass
                else:
                    parmSizer.Add(wx.StaticText(G2frame.FRMC, label=name), 0, WACV)
                    if name == "rcut":
                        parmSizer.Add(
                            G2G.ValidatedTxtCtrl(
                                G2frame.FRMC, RMCPdict, name, xmin=0.0, size=(70, 25)
                            ),
                            0,
                            WACV,
                        )
                        parmSizer.Add((5, 5))
                        continue
                    parmSizer.Add(
                        G2G.ValidatedTxtCtrl(
                            G2frame.FRMC, RMCPdict[name], 0, xmin=0.0, size=(70, 25)
                        ),
                        0,
                        WACV,
                    )
                    refine = wx.CheckBox(G2frame.FRMC, label="Refine")
                    refine.SetValue(RMCPdict[name][1])
                    refine.Bind(wx.EVT_CHECKBOX, OnRefine)
                    Indx[refine.GetId()] = name
                    parmSizer.Add(refine, 0, WACV)
            parmSizer.Add(wx.StaticText(G2frame.FRMC, label=" Shape"), 0, WACV)
            shape = wx.ComboBox(
                G2frame.FRMC,
                choices=["sphere", "stepcut"],
                style=wx.CB_DROPDOWN | wx.TE_READONLY,
            )
            shape.SetStringSelection(RMCPdict.get("shape", "sphere"))
            shape.Bind(wx.EVT_COMBOBOX, OnShape)
            parmSizer.Add(shape, 0, WACV)
            if RMCPdict.get("shape", "sphere") == "stepcut":
                for name in Names2:
                    parmSizer.Add(wx.StaticText(G2frame.FRMC, label=name), 0, WACV)
                    parmSizer.Add(
                        G2G.ValidatedTxtCtrl(
                            G2frame.FRMC, RMCPdict, name, xmin=0.0, size=(70, 25)
                        ),
                        0,
                        WACV,
                    )

            return parmSizer

        def OnSpaceGroup(event):
            # try a lookup on the user-supplied name
            SpcGp = G2phsG.GetSpGrpfromUser(G2frame.FRMC, SpGrp)
            SGErr, SGData = G2spc.SpcGroup(SpcGp)
            if SGErr:
                text = [G2spc.SGErrors(SGErr) + "\nSpace Group set to previous"]
                SGTxt.SetLabel(RMCPdict.get("Target SpGrp", "P 1"))
                msg = "Space Group Error"
                Text = "\n".join(text)
                wx.MessageBox(Text, caption=msg, style=wx.ICON_EXCLAMATION)
            else:
                text, table = G2spc.SGPrint(SGData)
                RMCPdict["SGData"] = SGData
                SGTxt.SetLabel(RMCPdict["SGData"]["SpGrp"])
                msg = "Target Space Group Information"
                G2G.SGMessageBox(G2frame.FRMC, msg, text, table).Show()
            G2spc.UpdateSytSym(data)

        def OnCellRef(event):
            RMCPdict["cellref"] = not RMCPdict["cellref"]

        def AtomSizer():
            def OnSetVal(event):
                r, c = event.GetRow(), event.GetCol()
                if c > 0:
                    strval = atmGrid.GetCellValue(r, c).strip()
                    try:
                        #                            if strval == '' or ('@' in strval and int(strval.split('@')[-1]) >= 20):
                        if strval == "" or "@" in strval:
                            RMCPdict["AtomConstr"][r][c + 1] = strval
                        else:
                            raise ValueError
                    except ValueError:
                        atmGrid.SetCellValue(r, c, RMCPdict["AtomConstr"][r][c + 1])
                        wx.MessageBox(
                            'ERROR - atom constraints must be blank or have "@n" with n >= 20',
                            style=wx.ICON_ERROR,
                        )
                    wx.CallAfter(UpdateRMC, G2frame, data)

            def OnUisoRefine(event):
                RMCPdict["UisoRefine"] = uiso.GetValue()
                nextP = 80
                oType = ""
                oName = ""
                for atom in RMCPdict["AtomConstr"]:
                    if RMCPdict["UisoRefine"] == "No":
                        atom[6] = ""
                    elif "same" in RMCPdict["UisoRefine"]:
                        atom[6] = "@81"
                        RMCPdict["AtomVar"]["@81"] = 0.005
                    elif "type" in RMCPdict["UisoRefine"]:
                        if atom[1] != oType:
                            oType = atom[1]
                            nextP += 1
                        atom[6] = "@%d" % nextP
                        RMCPdict["AtomVar"]["@%d" % nextP] = 0.005
                    elif "parent" in RMCPdict["UisoRefine"]:
                        if atom[0] != oName:
                            oName = atom[0]
                            nextP += 1
                        atom[6] = "@%d" % nextP
                        RMCPdict["AtomVar"]["@%d" % nextP] = 0.005
                wx.CallAfter(UpdateRMC, G2frame, data)

            atmSizer = wx.BoxSizer(wx.VERTICAL)
            atmSizer.Add(
                wx.StaticText(
                    G2frame.FRMC,
                    label=' Atom Constraints; enter as e.g. "@n" or "0.5-@n"; n>=20 && "@n" should be at end',
                )
            )
            uisoSizer = wx.BoxSizer(wx.HORIZONTAL)
            uisoSizer.Add(wx.StaticText(G2frame.FRMC, label=" Refine Uiso? "), 0, WACV)
            uiso = wx.ComboBox(
                G2frame.FRMC,
                choices=["No", "by type", "by parent name", "all same"],
                style=wx.CB_DROPDOWN | wx.TE_READONLY,
            )
            uiso.SetValue(RMCPdict["UisoRefine"])
            uiso.Bind(wx.EVT_COMBOBOX, OnUisoRefine)
            uisoSizer.Add(uiso, 0, WACV)
            atmSizer.Add(uisoSizer)

            table = [item[1:] for item in RMCPdict["AtomConstr"]]
            addCol = False
            if len(RMCPdict["AtomConstr"][0]) > 6:
                addCol = True
            colLabels = [
                "Type",
                "x constraint",
                "y constraint",
                "z  constraint",
                "frac constr",
                "Uiso constr",
            ]
            rowLabels = [item[0] for item in RMCPdict["AtomConstr"]]
            Types = 6 * [
                wg.GRID_VALUE_STRING,
            ]
            if addCol:
                colLabels += [
                    "sym opr",
                ]
                Types = 7 * [
                    wg.GRID_VALUE_STRING,
                ]
            atmTable = G2G.Table(
                table, rowLabels=rowLabels, colLabels=colLabels, types=Types
            )
            atmGrid = G2G.GSGrid(G2frame.FRMC)
            atmGrid.SetTable(atmTable, True, useFracEdit=False)
            atmGrid.AutoSizeColumns(True)
            atmGrid.Bind(wg.EVT_GRID_CELL_CHANGED, OnSetVal)
            atmSizer.Add(atmGrid)
            return atmSizer

        def AtomVarSizer():
            atomVarSizer = wx.FlexGridSizer(0, 8, 5, 5)
            for item in RMCPdict["AtomVar"]:
                atomVarSizer.Add(wx.StaticText(G2frame.FRMC, label=item), 0, WACV)
                atomVarSizer.Add(
                    G2G.ValidatedTxtCtrl(
                        G2frame.FRMC,
                        RMCPdict["AtomVar"],
                        item,
                        xmin=-3.0,
                        xmax=3.0,
                        size=(70, 25),
                    ),
                    0,
                    WACV,
                )
            return atomVarSizer

        txt = wx.StaticText(
            G2frame.FRMC,
            label="For use of PDFfit, please cite: " + G2G.GetCite("PDFfit2"),
        )
        txt.Wrap(500)
        mainSizer.Add(txt)
        mainSizer.Add((5, 5))
        if (
            "PDFfit" not in data["RMC"]
            or not data["RMC"]["PDFfit"]
            or "delta1" not in data["RMC"]["PDFfit"]
        ):
            if "PDFfit" not in data["RMC"]:
                data["RMC"]["PDFfit"] = {}
            metadata = {
                "title": "none",
                "date": str(time.ctime()),
                "temperature": "300K",
                "doping": 0,
            }
            files = {
                "Neutron real space data; G(r): ": [
                    "Select",
                    0.05,
                    "G(r)",
                    "RMC",
                ],
                "Xray real space data; G(r): ": [
                    "Select",
                    0.01,
                    "G(r)",
                    "RMC",
                ],
            }
            data["RMC"]["PDFfit"].update(
                {
                    "files": files,
                    "ReStart": [False, False],
                    "metadata": metadata,
                    "delta1": [0.0, False],
                    "delta2": [0.0, False],
                    "spdiameter": [0.0, False],
                    "refinement": "normal",
                    "sratio": [1.0, False],
                    "rcut": 0.0,
                    "stepcut": 0.0,
                    "shape": "sphere",
                    "cellref": False,
                    "SeqDataType": "X",
                    "SeqCopy": True,
                    "SeqReverse": False,
                    "Xdata": {
                        "dscale": [1.0, False],
                        "Datarange": [0.0, 30.0],
                        "Fitrange": [0.0, 30.0],
                        "qdamp": [0.03, False],
                        "qbroad": [0.0, False],
                    },
                    "Ndata": {
                        "dscale": [1.0, False],
                        "Datarange": [0.0, 30.0],
                        "Fitrange": [0.0, 30.0],
                        "qdamp": [0.03, False],
                        "qbroad": [0.0, False],
                    },
                }
            )

        RMCPdict = data["RMC"]["PDFfit"]
        # patch
        if "AtomConstr" not in RMCPdict:  # keep this one
            RMCPdict["AtomConstr"] = []
        if "AtomVar" not in RMCPdict:
            RMCPdict["AtomVar"] = {}
        if "SGData" not in RMCPdict:
            RMCPdict["SGData"] = G2spc.SpcGroup("P 1")[1]
        if "refinement" not in RMCPdict:
            RMCPdict["refinement"] = "normal"
        if "cellref" not in RMCPdict:
            RMCPdict["cellref"] = False
        if "metadata" not in RMCPdict:
            RMCPdict["metadata"] = {
                "title": "none",
                "date": str(time.ctime()),
                "temperature": "300K",
                "doping": 0,
            }
        if "SeqDataType" not in RMCPdict:
            RMCPdict["SeqDataType"] = "X"
        if "SeqCopy" not in RMCPdict:
            RMCPdict["SeqCopy"] = False
            RMCPdict["SeqReverse"] = False
        if "UisoRefine" not in RMCPdict:
            RMCPdict["UisoRefine"] = "No"
        # end patch
        Atoms = data["Atoms"]
        cx, ct, cs, ci = G2mth.getAtomPtrs(data)
        if not RMCPdict["AtomConstr"]:
            for atom in Atoms:
                RMCPdict["AtomConstr"].append(
                    [atom[ct - 1], atom[ct], "", "", "", "", ""]
                )
        else:  # update name/type changes
            for iatm, atom in enumerate(Atoms):
                RMCPdict["AtomConstr"][iatm][:2] = atom[ct - 1 : ct + 1]

        mainSizer.Add(wx.StaticText(G2frame.FRMC, label=" Enter metadata items:"), 0)
        mainSizer.Add(
            GetMetaSizer(RMCPdict, ["title", "date", "temperature", "doping"]), 0
        )

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        SgSizer = wx.BoxSizer(wx.HORIZONTAL)
        SgSizer.Add(wx.StaticText(G2frame.FRMC, label=" Target space group: "), 0, WACV)

        mainSizer.Add(
            wx.StaticText(G2frame.FRMC, label="PDFfit phase structure parameters:")
        )

        SpGrp = RMCPdict["SGData"]["SpGrp"]
        SGTxt = wx.Button(G2frame.FRMC, wx.ID_ANY, SpGrp, size=(100, -1))
        SGTxt.Bind(wx.EVT_BUTTON, OnSpaceGroup)
        SgSizer.Add(SGTxt, 0, WACV)
        mainSizer.Add(SgSizer)

        cellref = wx.CheckBox(G2frame.FRMC, label=" Refine unit cell?")
        cellref.SetValue(RMCPdict["cellref"])
        cellref.Bind(wx.EVT_CHECKBOX, OnCellRef)
        mainSizer.Add(cellref)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(wx.StaticText(G2frame.FRMC, label="PDFfit atom parameters:"))
        mainSizer.Add(AtomSizer())

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(
            wx.StaticText(G2frame.FRMC, label="PDFfit starting atom variables:")
        )
        G2pwd.GetPDFfitAtomVar(data, RMCPdict)
        mainSizer.Add(AtomVarSizer())

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(
            wx.StaticText(G2frame.FRMC, label=" PDFfit phase profile coefficients:")
        )
        mainSizer.Add(PDFParmSizer(), 0)

        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        mainSizer.Add(FileSizer(RMCPdict))
        return mainSizer

    ####start of UpdateRMC
    G2frame.GetStatusBar().SetStatusText("", 1)
    G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC, False)
    G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC, False)
    if G2frame.RMCchoice == "RMCProfile":
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC, True)
        # G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_STOPRMC,False)
    elif G2frame.RMCchoice == "fullrmc":
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SETUPRMC, False)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_RUNRMC, False)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_VIEWRMC, False)
        # G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_STOPRMC,False)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_ATOMSRMC, True)
        G2frame.dataWindow.FRMCDataEdit.Enable(G2G.wxID_SUPERRMC, True)
    try:
        if G2frame.FRMC.GetSizer():
            G2frame.FRMC.GetSizer().Clear(True)
    except:  # wxAssertionError from C++
        pass
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    if not len(data["Atoms"]):
        mainSizer.Add(
            wx.StaticText(
                G2frame.FRMC, label="No atoms found - PDF fitting not possible"
            )
        )
    else:
        runFile = " "
        choice = ["RMCProfile", "fullrmc", "PDFfit"]
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        RMCsel = wx.RadioBox(G2frame.FRMC, -1, " Select RMC method:", choices=choice)
        RMCsel.SetStringSelection(G2frame.RMCchoice)
        RMCsel.Bind(wx.EVT_RADIOBOX, OnRMCselect)
        topSizer.Add(RMCsel, 0)
        topSizer.Add((20, 0))
        txt = wx.StaticText(
            G2frame.FRMC,
            label="NB: if you change any of the entries below, you must redo the Operations/Setup RMC step above to apply them before doing Operations/Execute",
        )
        txt.Wrap(250)
        topSizer.Add(txt, 0)
        mainSizer.Add(topSizer, 0)
        RMCmisc["RMCnote"] = wx.StaticText(G2frame.FRMC)
        mainSizer.Add(RMCmisc["RMCnote"])
        G2G.HorizontalLine(mainSizer, G2frame.FRMC)
        if G2frame.RMCchoice == "fullrmc":
            RMCPdict = data["RMC"]["fullrmc"]
            mainSizer.Add(fullrmcSizer(RMCPdict))

        elif G2frame.RMCchoice == "RMCProfile":
            RMCPdict = data["RMC"]["RMCProfile"]
            mainSizer.Add(RMCProfileSizer(RMCPdict))

        else:  # PDFfit
            mainSizer.Add(PDFfitSizer(data))

    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl = f"PDF fitting options {data['General']['Name']!r}"[:60]
    topSizer.Add(wx.StaticText(parent, label=lbl), 0, WACV)
    topSizer.Add((-1, -1), 1, wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent, helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    G2phsG.SetPhaseWindow(G2frame.FRMC, mainSizer)

    if G2frame.RMCchoice == "PDFfit" and not G2phsG.checkPDFfit(G2frame):
        RMCmisc["RMCnote"].SetLabel("PDFfit may not be installed or operational")
    elif G2frame.RMCchoice == "fullrmc" and G2pwd.findfullrmc() is None:
        msg = (
            "The fullrmc Python image is not found."
            + " Do you want it installed for you from "
            + "https://github.com/bachiraoun/fullrmc/tree/master/standalones?"
            + "\n\n40-50 Mb (download times vary)"
        )
        dlg = wx.MessageDialog(G2frame, msg, "Install fullrmc", wx.YES | wx.NO)
        try:
            dlg.CenterOnParent()
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_YES:
            wx.BeginBusyCursor()
            G2pwd.fullrmcDownload()
            wx.EndBusyCursor()
        else:
            RMCmisc["RMCnote"].SetLabel(
                "Note that fullrmc is not installed or was not located"
            )
