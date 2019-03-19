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
This script imports GSASIIpath, does some minor initialization 
and then launches :func:`GSASIIdataGUI.GSASIImain`, 
which creates a wx.Application that in turns creates the GUI. 
If the GSAS-II binaries are not installed or are incompatible with
the OS/Python packages, the user is asked if they should be updated
from the subversion site. 
'''

import platform
import GSASIIpath

__version__ = '1.0.0'

if __name__ == '__main__':
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
    GSASIIpath.SetVersionNumber("$Revision$")
    GSASIIpath.InvokeDebugOpts()
    import GSASIIdataGUI as G2gd    
    G2gd.GSASIImain() # start the GUI
