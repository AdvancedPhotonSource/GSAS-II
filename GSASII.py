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
*GSASII: GSAS-II GUI*
=====================

File GSASII.py is the script to start the GSAS-II graphical user 
interface (GUI). 
This script imports GSASIIpath, which does some minor initialization
and then (before any wxPython calls can be made) creates a wx.App application. 
A this point :func:`GSASIIpath.SetBinaryPath` is called to establish
the directory where GSAS-II binaries are found. If the binaries 
are not installed or are incompatible with the OS/Python packages, 
the user is asked if they should be updated from the subversion site. 
The wxPython app is then passed to :func:`GSASIIdataGUI.GSASIImain`, 
which creates the GSAS-II GUI and finally the event loop is started.

Keyboard Menu Shortcuts
----------------------------------------

Shortcuts for commonly-used menu commands are created by adding a 
menu command with a "\tctrl+" addition such as::

        item = parent.Append(wx.ID_ANY,'&Refine\tCtrl+R','Perform a refinement')

This will allow the above menu command to be executed with a "Control-R" 
keyboard command (on MacOS this will be "Command+R" rather than "Control-R") as well as using the menu to access that action. The following table lists the 
keyboard letters/numbers that have GSAS-II assigned actions.
are system assigned. Note that there are also plotting keyboard commands that are 
implemented in :mod:`GSASIIplot`. 
These can be discovered from the "K" button on the plot menu bar as they 
vary depending on the type of plot.

.. tabularcolumns:: |c|p{4in}|

==========  ====================================================
  key         explanation
==========  ====================================================
 O           Open project (File menu)
 E           Reopen recent (File menu)
 S           Save project (File menu)
 B           Project browser (File menu)
 Q           Quit (File menu). This is system assigned on MacOS  
 F4          Quit (File menu). This is system-assigned 
             on Windows  

 L           View LS parms (Calculate menu)
 R           Refine/Sequential Refine (Calculate menu)
 I           Parameter Impact (Calculate menu)

 U           Check for updates (Help menu)
 T           Tutorials (Help menu)
 F1          Help on current tree item (Help menu).
             This is system-assigned 

 P           Peakfit (Peak Fitting menu, requires selection of 
             Histogram Peak)

 M           Minimize GSAS-II windows (MacOS Windows menu).
             This is system-assigned 
==========  ====================================================

'''

import sys
#import os
import platform
import scipy.optimize # addresses problem with build for wx on Pi
try:
    import wx
# importing the following wx modules at the same time as wx seems to eliminate 
# the "Debug: Adding duplicate image handler for 'Windows bitmap file'"
# error message
    import wx.grid as wg
    import wx.aui
    import wx.lib.scrolledpanel as wxscroll
    import wx.html        # could postpone this for quicker startup
    import wx.lib.mixins.listctrl  as  listmix
    import wx.richtext as wxrt
    import wx.lib.filebrowsebutton as wxfilebrowse
except ImportError:
    pass
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
