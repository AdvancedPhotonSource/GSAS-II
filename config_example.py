# -*- coding: utf-8 -*-
#config.py - Variables used to set optional configuration options
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*config_example.py: Configuration options*
-------------------------------------------

This file contains optional configuration options for GSAS-II. The variables
in this file can be copied to file config.py, which is imported if present.
Access these variables using :func:`GSASIIpath.GetConfigValue`, which returns
None if the variable is not set. Note that a config.py file need not
be present, but if in use it will typically be found with the GSAS-II source
directory (GSASIIpath.Path2GSAS2) or a directory for local GSAS-II
modifications (~/.G2local/ or /Documents and Settings/<User>/.G2local/).  

When defining new config variables for GSAS-II, define them here with a
default value: use None or a string for strings, or use integers or real
values. Include a doc string after each variable is defined to explain
what it does. Use names ending in _location or _directory for items
that will contain directory names.

For example::

    test_int = 0
    test_float = 0.0
    test_string = None (or)
    test_string = 'value'
'''

debug = False
'''Set to True to turn on debugging mode.This enables use of IPython on
exceptions and on calls to :func:`GSASIIpath.IPyBreak`. Calls to
:func:`GSASIIpath.pdbBreak` will invoke pdb at that location.

If debug is False, calls to :func:`GSASIIpath.IPyBreak` and
:func:`GSASIIpath.pdbBreak` are ignored.
'''

Clip_on = True
''' if True then line plots willl be clipped at plot border; 
if False line plots extend nto white space around plot frme
'''

Transpose = False
'Set to True to cause images to be Transposed when read (for code development)'

Enable_logging = False
'Set to True to enable use of command logging (under development.)'

logging_debug = False
'Set to True to enable debug for logging (under development.)'

Help_mode = "browser"
'''Set to "internal" to use a Python-based web viewer to display
help documentation and tutorials. If set to the default ("browser")
the default web browser is used.
'''

Tutorial_location = None
'''Change this to place tutorials by in a different spot. If None, this defaults to
<user>/My Documents/G2tutorials (on windows) or <user>/G2tutorials. If you want to
use a different location, this can be set here. To install into the location where
GSAS-II is installed, use this::

    Tutorial_location = GSASIIpath.path2GSAS2

As another example, to use ~/.G2tutorials do this::

    Tutorial_location = '~/.G2tutorials'

Note that os.path.expanduser is run on Tutorial_location before it is used.
Also note that GSASIIpath is imported inside config.py; other imports should be
avoided.
'''

Save_paths=False
'''When set to True, the last-used path for saving of .gpx and for
importing of input files is saved in the configuration file.
Note that since this causes the config.py file to be updated whenever files are
saved/imported, any temporary config settings can be saved to disk at that
point.
'''

Starting_directory=None
'''Specifies a default location for starting GSAS-II and where .gpx files
should be read from. Will be updated if Save_paths is True.
Note that os.path.expanduser is run on this before it is used, so the user's
home directory can be specified with a '~'.
'''

Import_directory=None
'''Specifies a default location for importing (reading) input files. Will be
updated if Save_paths is True.
Note that os.path.expanduser is run on this before it is used, so the user's
home directory can be specified with a '~'.
'''

wxInspector = False
'''If set to True, the wxInspector widget is displayed when
GSAS-II is started.
'''
