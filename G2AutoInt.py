'''
*G2AutoInt: independent autointegration tool*
---------------------------------------------

Independent-running GSAS-II based auto-integration program with minimal
GUI, no visualization but intended to implement significant levels of 
parallelization.
'''
from __future__ import division, print_function
#import platform
import time
import math
#import random as ran
import copy
import sys
import os
import wx
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIautoInt as G2imG
import GSASIIfiles as G2fil
GSASIIpath.InvokeDebugOpts()

App = wx.App()
class dummyClass(object):
    def __init__(self):
        # find all the exporter files
        self.exporterlist = G2fil.LoadExportRoutines(self)
        self.Image = None
        self.GSASprojectfile = '/tmp/x.gpx'
        self.LastExportDir = ''
        self.LastGPXdir = ''
        self.PauseIntegration = False
    pass

if __name__ == "__main__":
    G2frame = dummyClass()
    frm = G2imG.AutoIntFrame(G2frame,5)
    App.GetTopWindow().Show(True)
    App.MainLoop()
