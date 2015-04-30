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

Enable_logging = False
'Set to True to enable use of command logging'

logging_debug = False
'Set to True to enable debug for logging'

Help_mode = None
'Set to "internal" to use a Python-based web viewer rather than a web browser'

Tutorial_location = None
'''Change this to place tutorials by in a different spot. If None, this defaults to
~/My Documents/G2tutorials (on windows) or ~/G2tutorials. If you want to use a different
location (such as to use the GSASII installation directory), this can be set here.
As an example, to always use ~/G2tutorials do this::

    import os.path
    Tutorial_location = os.path.join(os.path.expanduser('~'),'G2tutorials')

To install into the location where GSAS-II is installed, use this::

    import GSASIIpath
    Tutorial_location = GSASIIpath.path2GSAS2

'''

Starting_directory=None
'''Specifies a default location for starting GSAS-II
'''

Import_directory=None
'''Specifies a default location for starting GSAS-II
'''

wxInspector = None
'''If set to True, the wxInspector widget is displayed
'''
