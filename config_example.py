# -*- coding: utf-8 -*-
#config.py - Variables used to set optional configuration options
########### SVN repository information ###################
# $Date: $
# $Author: toby $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*config.py: Configuration options*
----------------------------------

This file contains optional configuration options for GSAS-II. Note that
this file is not required to be present and code should be written to
provide the default behavior if the file is not present or if any configuration
variable is not set, but please do place a docstring here for every used
config variable explaining what it does. Access these variables using
:func:`GSASIIpath.GetConfigValue`.
'''

debug = False
'''Set to True to turn on debugging mode. This enables use of IPython on 
exceptions and on calls to :func:`GSASIIpath.IPyBreak`. Calls to
:func:`GSASIIpath.pdbBreak` will invoke pdb at that location.
If debug is false calls to :func:`GSASIIpath.IPyBreak` and
:func:`GSASIIpath.pdbBreak` are ignored.
'''

Enable_logging = None
'Set to True to enable use of command logging'

logging_debug = None
'Set to True to enable debug for logging'
