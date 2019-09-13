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

import platform
import wx
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

__version__ = '1.0.0'

if __name__ == '__main__':
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
    GSASIIpath.InvokeDebugOpts()
    import GSASIIdataGUI as G2gd    
    G2gd.GSASIImain(application) # start the GUI
    application.MainLoop()
