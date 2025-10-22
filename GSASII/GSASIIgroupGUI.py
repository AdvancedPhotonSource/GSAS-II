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
# from . import GSASIIspc as G2spc
# from . import GSASIIlattice as G2lat
# from . import GSASIIpwd as G2pwd
from . import GSASIIctrlGUI as G2G
from . import GSASIIpwdplot as G2pwpl
WACV = wx.ALIGN_CENTER_VERTICAL

def UpdateGroup(G2frame,item):
    G2gd.SetDataMenuBar(G2frame)
    G2frame.dataWindow.helpKey = "Groups/PWDR"
    topSizer = G2frame.dataWindow.topBox
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label=' Group edit goes here someday'),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    G2G.HorizontalLine(G2frame.dataWindow.GetSizer(),G2frame.dataWindow)
    #G2frame.dataWindow.GetSizer().Add(text,1,wx.ALL|wx.EXPAND)
    G2frame.groupName = G2frame.GPXtree.GetItemText(item)
    G2pwpl.PlotPatterns(G2frame,plotType='GROUP')
