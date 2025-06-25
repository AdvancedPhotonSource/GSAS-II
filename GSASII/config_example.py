# -*- coding: utf-8 -*-
#config.py - Variables used to set optional configuration options
'''
This file contains optional configuration options for GSAS-II. The 
values for the variables named here will be set in file ~/.GSASII/config.ini
which is read on startup by :func:`GSASIIpath.LoadConfig`. To check
if a configuration variable has been set use 
:func:`GSASIIpath.GetConfigValue`, which returns
None if the variable is not set.
Values are typically changed using :func:`GSASIIctrlGUI.SelectConfigSetting`
which uses :func:`GSASIIctrlGUI.SaveConfigVars` to write the 
~/.GSASII/config.ini file. 

To define new config variables for GSAS-II, define them here with a
default value: use None or a string for strings. If an integer or real
values is used as a default, the routines will ensure that this type
is preserved for any user setting. 
Always include a doc string after defining each variable. This definition 
will be shown in the GUI to explain what the variable does. 

If a name ends with a particular keyword, then specialized edit 
routines are used.

* Names ending in _location or _directory are for path items
* Names ending in _exec for executable files (.exe on windows).
* Names ending in _color for colors, to be specified as RGBA values 
  (note that Contour_color is restricted to color maps). 
* Names ending in _pos or _Size are integer tuples for wx sizes or positions.

For example::

    int_config = 0
    float_config = 0.0
    string_config = None (or)
    string_config = 'value'
'''

debug = False
'''Set to True to turn on debugging mode. This enables use of IPython on
exceptions and on calls to :func:`GSASIIpath.IPyBreak` or breakpoint(). 
Calls to :func:`GSASIIpath.pdbBreak` will invoke pdb at that location.
%%
If debug is False, calls to :func:`GSASIIpath.IPyBreak`, breakpoint() and
:func:`GSASIIpath.pdbBreak` are ignored.
%%
From inside Spyder, calls to breakpoint() invoke the Spyder debugger, 
independent of the setting of debug. 
%%
Restart GSAS-II for the setting of debug to take effect.
'''

Clip_on = True
''' if True then line plots willl be clipped at plot border; 
if False line plots extend nto white space around plot frme
'''

Transpose = False
'Set to True to cause images to be Transposed when read (for code development)'

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

Spot_mask_diameter = 1.0
'''Specifies the default diameter for creation of spot masks. Default is 1.0 mm
'''

Ring_mask_thickness = 0.1
'''Specifies the default thickness for creation of ring and arc masks.
Default is 0.1 degrees 2-theta.
'''

Arc_mask_azimuth = 10.0
'''Specifies the default azimuthal range for creation of arc masks.
Default is 10.0 degrees 2-theta.
'''

Autoint_PollTime = 30.
'''Specifies the frequency, in seconds that AutoInt checks for new files.
Default is 30 seconds
'''

Autoscale_ParmNames = ['userComment2',r'extraInputs\1\extraInputs','Ion_Chamber_I0',]
'''Gives the possible selection of incident monitor names as found in an image metadata file.
Used in AutoIntegration
'''
DefaultAutoScale = "userComment2"
'''DefaultAutoScale selects one of the AutoScale_ParmNames.
Used in AutoIntegration
'''
Main_Size = (700,450)
'''Main window size (width, height) - initially uses wx.DefaultSize but will updated
and saved as the user changes the window.
This is used internally by GSAS-II and would not normally be changed by a user. 
'''
Main_Pos = (100,100)
'''Main window location - will be updated & saved when user moves
it. If position is outside screen then it will be repositioned to default.
This is used internally by GSAS-II and would not normally be changed by a user. '''
Plot_Size = (700,600)
'''Plot window size (width, height) - initially uses wx.DefaultSize but will 
updated and saved as the user changes the window.
This is used internally by GSAS-II and would not normally be changed by a user. 
'''
Plot_Pos = (200,200)
'''Plot window location - will be updated & saved when user moves it
these widows. If position is outside screen then it will be repositioned to default.
This is used internally by GSAS-II and would not normally be changed by a user. 
'''

Tick_length = 8.0
'''Specifies the length of phase tick marks in pixels. Default is 8.'''

Tick_width = 1.0
'''Specifies the width of phase tick marks in pixels.
Fractional values do seem to produce an effect. Default is 1.'''

Contour_color = 'GSPaired'
''' Specifies the color map to be used for contour plots (images, pole figures, etc.)
will be applied for new images and if Saved for a new start of GSAS-II
'''

Movie_fps = 10
''' Specifies movie frames-per-second; larger number will make smoother modulation movies but larger files.
'''

Movie_time = 5
''' Specifices time in sec for one modulation loop; larger number will give more frames for same fps'
'''
fullIntegrate = True
''' If True then full image integration is default; False otherwise
'''

Multiprocessing_cores = 0
''' Specifies the number of cores to use when performing multicore computing. A number less
than zero causes the recommended number of cores [using multiprocessing.cpu_count()/2]
to be used. Setting this number to 0 or 1 avoids use of the multiprocessing module: all
computations are performed in-line. 
'''

Show_timing = False
'''If True, shows various timing results.'''

Column_Metadata_directory = None
'''When specified and when images are read, GSAS-II will read metadata from a 1-ID
style .par and a .EXT_lbls (EXT = image extension) or .lbls file. See :func:`GSASIIfiles.readColMetadata` for
information on how this is done.
'''

Instprm_default = False
'''when True, GSAS-II instprm file are shown as default; when False, old GSAS stype prm, etc files are default
'''

Ref_Colors = 'b r c g m k'
'''The colors for reflection tick marks by phase. 
Use one of 'k'-black, 'r'-red, 'b'-blue, 'g'-green, 'm'-magenta, 'c'-cyan 
for the line colors, or any other valid matplotlib color name or hex code. 
'''

Obs_color = '0000ffff'
'''The color for plotting the observed powder diffraction pattern.
Colors are specified as hex RGBA values, as used in Matplotlib (without 
preceding #). 
The default is 0000ffff, which sets the color to blue. 
'''

Calc_color = '008000ff'
'''The color for plotting the computed powder diffraction pattern.
Colors are specified as hex RGBA values, as used in Matplotlib (without 
preceding #). 
The default is 00ff00ff, which sets the color to green. 
'''

Diff_color = '00bfbfff'
'''The color for plotting the obs-calc powder diffraction pattern.
Colors are specified as hex RGBA values, as used in Matplotlib (without 
preceding #). 
The default is 00ffffff, which sets the color to cyan. 
'''

Bkg_color = 'ff0000ff'
'''The color for plotting the background powder diffraction pattern.
Colors are specified as hex RGBA values, as used in Matplotlib (without 
preceding #). 
The default is ff0000ff, which sets the color to red.
'''

PDF_Rmax = 100.
'''Maximum radius for G(r) calculations: range is from 10-200A; default is 100A
'''

previous_GPX_files = []
'''A list of previously used .gpx files. This is used internally by GSAS-II 
and would not normally be changed by a user. 
'''

Image_calibrant = ''
''' Specifies a default calibrant material for images. Will be applied for 
newly-read images, but if changed the specified material will be saved.
'''

Image_2theta_min = 5.0
''' Specifies a default 2-theta minimum used for calibration and integration
as the Inner 2-theta value. Will be applied for 
newly-read images, but if changed the new value will be saved.
'''

Image_2theta_max = 50.0
''' Specifies a default 2-theta maximum used for calibration and integration
as the Outer 2-theta value. Will be applied for 
newly-read images, but if changed the new value will be saved.
'''

enum_DrawAtoms_default = ['','lines','vdW balls','sticks','balls & sticks','ellipsoids',]
'choices for DrawAtoms_default'
DrawAtoms_default = ''
'''Allows selection of the default plotting mode for structures
in Draw Atoms. The only valid values are:
'lines', 'vdW balls', 'sticks', 'balls & sticks', 'ellipsoids'.
%% If a non-valid choice is used (the default)
'vdW balls' is used.
'''

show_gpxSize = False
'''When True, the sizes of the sections of the GPX file are listed
when the GPX file is opened. Default is False.
'''

fullrmc_exec = None
'''Defines the full path to a Python executable that has been configured 
with the fullrmc package. If None (the default), GSAS-II will see if fullrmc
can be imported into the current Python (which is unlikely to ever work). 
If that does not work, GSAS-II will search for an executable named fullrmc* 
(or fullrmc*.exe on Windows) in the Python ``sys.path`` search path,
which includes the GSAS-II binary directory.
'''

pdffit2_exec = None
'''Defines the full path to a Python executable that has been configured 
with the PDFfit2 (diffpy) package. If None (the default), GSAS-II will see 
if PDFfit2 can be imported into the current Python.
'''

lastUpdateNotice = 0
'''Defines the version number for the last update notice that has been 
shown. This should not need to be changed manually.
'''

SeparateHistPhaseTreeItem = False
'''When this is set to True, the parameters specific to each histogram 
and phase together (such as peak shapes & phase fractions) 
are shown as a 1st-level tree item rather than inside each Phase's
Data tab. After changing this, GSAS-II needs to be restarted for the 
change to take effect. Default is False.
'''

G2RefinementWindow = False
'''When True a custom progress window is displayed to track the 
progress of refinements. When False a generic wxpython supplied progress
dialog is used. 
'''

HDF5selection = 20
'''When an HDF5 file contains more than this number of images, a selection
window is offered to determine which images will be read. If negative, 
the selection window will never be used. If zero, it will always be used.
'''

FontSize_incr = 0
'''Specifies a point size to increase (or decrease if negative) the default 
font size for the GSAS-II windows. Default is 0. An increment much larger than 
~4 will likely cause some places where text no longer fits, but might be useful on high resolution monitors.
%%
Restart GSAS-II for this setting to take effect.
'''

G2FileBrowser = False
'''When set to True, the GSAS-II provided file browser is used to find
files when files are imported. For Linux the default is True, but
for Windows and Mac, the default is False
'''
