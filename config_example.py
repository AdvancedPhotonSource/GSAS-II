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
the location where GSAS-II is loaded GSASIIpath.path2GSAS2. For installations where
G2 is installed by an administrator, it is a good idea to use something like this::

    import os.path
    Tutorial_location = os.path.join(os.path.expanduser('~'),'G2tutorials')

This will allow users to download tutorial files into their own file space.
'''

Starting_directory=None
'''Specifies a default location for starting GSAS-II
'''

Import_directory=None
'''Specifies a default location for starting GSAS-II
'''
