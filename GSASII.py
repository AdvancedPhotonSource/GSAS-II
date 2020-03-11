#!/usr/bin/env python
# -*- coding: utf-8 -*-
#GSASII
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSAS-II GUI*
=====================

This is the script to start the GSAS-II graphical user interface (GUI). 
This script imports GSASIIpath, which does some minor initialization
and then (before any wxPython calls can be made) creates a wx.App application. 
A this point :func:`GSASIIpath.SetBinaryPath` is called to establish
the directory where GSAS-II binaries are found. If the binaries 
are not installed or are incompatible with the OS/Python packages, 
the user is asked if they should be updated from the subversion site. 
The wxPython app is then passed to :func:`GSASIIdataGUI.GSASIImain`, 
which creates the GSAS-II GUI and finally the event loop is started.
'''

import sys
#import os
import platform
import wx
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

__version__ = '1.0.0'

class G2App(wx.App):
    '''Used to create a wx python application for the GUI for Mac.
    Customized to implement drop of GPX files onto app.
    '''
    startupMode = True
    def ClearStartup(self):
        '''Call this after app startup complete because a Drop event is posted 
        when GSAS-II is initially started.
        '''
        self.startupMode = False        
    def MacOpenFiles(self, filenames):
        if self.startupMode:
            return
        for project in filenames:
            #print("Start GSAS-II with project file "+str(project))
            GSASIIpath.MacStartGSASII(__file__,project)

if __name__ == '__main__':
    if sys.platform == "darwin": 
        application = G2App(0) # create the GUI framework
    else:
        application = wx.App(0) # create the GUI framework
    try:
        GSASIIpath.SetBinaryPath(True)
    except:
        print('Unable to run with current setup, do you want to update to the')
        try:
            if '2' in platform.python_version_tuple()[0]:            
                ans = raw_input("latest GSAS-II version? Update ([Yes]/no): ")
            else:
                ans = input("latest GSAS-II version? Update ([Yes]/no): ")                
        except:
            ans = 'no'
        if ans.strip().lower() == "no":
            import sys
            print('Exiting')
            sys.exit()
        print('Updating...')
        GSASIIpath.svnUpdateProcess()
    import GSASIIdataGUI as G2gd    
    G2gd.GSASIImain(application) # start the GUI
    if sys.platform == "darwin": 
        wx.CallLater(100,application.ClearStartup)
    GSASIIpath.InvokeDebugOpts()
    application.MainLoop()
