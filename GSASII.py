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
*GSAS-II GUI Script*
=====================

Script to start the GSAS-II graphical user interface. This script imports GSASIIpath,
sets a few misc. values and then launches :func:`GSASIIdataGUI.GSASIImain`, which
creates a wx.Application which in turns creates the GUI. 
'''

import GSASIIpath

__version__ = '1.0.0'

if __name__ == '__main__':
    GSASIIpath.SetBinaryPath()
    GSASIIpath.SetVersionNumber("$Revision$")
    GSASIIpath.InvokeDebugOpts()
    import GSASIIdataGUI as G2gd    
    G2gd.GSASIImain() # start the GUI
