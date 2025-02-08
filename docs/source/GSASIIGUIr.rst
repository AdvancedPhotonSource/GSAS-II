*GSAS-II GUI Support Modules*
==============================

The modules documented here provide GUI or graphics capabilities that
are used in multiple sections of the GSAS-II GUI or graphics. 

---------------------------------------------
*GSASIIctrlGUI: Custom GUI controls*
---------------------------------------------

 .. py:currentmodule:: GSASIIctrlGUI

A library of specialized widgets (GUI controls) for use in the GSAS-II
data window or dialogs:

.. tabularcolumns:: |l|p{4in}|

================================  =================================================================
Class or function name             Description
================================  =================================================================
:class:`EnumSelector`              A combo box with a built-in call back routine that
                                   automatically sets a dict or list entry.
:class:`G2ChoiceButton`            A customized wx.Choice that automatically initializes to
                                   the initial value and saves the choice directly into a dict
                                   or list value. Optionally calls function when a
                                   choice is selected
:class:`G2CheckBox`                A customized wx.CheckBox that automatically initializes to
                                   the initial value and saves the choice directly into a dict
                                   or list value. Optionally calls function when a
                                   choice is selected
:func:`G2CheckBoxFrontLbl`         A version of :class:`G2CheckBox` that places the label
                                   for the check box in front. Otherwise works the same. 
:func:`G2RadioButtons`             Creates a series of grouped radio buttons.
:class:`G2SliderWidget`            A customized combination of a wx.Slider and a validated 
                                   wx.TextCtrl (see :class:`ValidatedTxtCtrl`).
:class:`G2Slider`                  A wrapped version of wx.Slider that implements scaling
:class:`G2SpinWidget`              A customized combination of a wx.SpinButton and a validated 
                                   wx.TextCtrl (see :class:`ValidatedTxtCtrl`).
:class:`G2MultiChoiceWindow`       Similar to :class:`G2MultiChoiceDialog` but provides
                                   a sizer that can be placed in a frame or panel.
:class:`HelpButton`                Creates a button labeled with a "?" that when pressed
                                   displays help text in a modal message window
                                   or web browser. 
:class:`OrderBox`                  Creates a wx.Panel with scrollbars where items can be
                                   ordered into columns.
:class:`SortableLstCtrl`           Creates a wx.Panel for a table of data that  
                                   can be sorted by clicking on a column label.
:class:`ValidatedTxtCtrl`          A text control with a built-in call back routine to set dict
                                   or list elements. Optionally validates input as float, int or
                                   for strings non-blank. Value is set when focus changes
:func:`HorizontalLine`             Places a line in a Frame or Dialog to separate sections.
:class:`ScrolledStaticText`        A wx.StaticText widget that fits a large string into a 
                                   small space by scrolling it
:func:`ReadOnlyTextCtrl`           A wx.TextCtrl widget to be used wx.StaticText 
                                   (no edits allowed) text appears in a box.
:func:`setColorButton`             A button for color selection as a replacement 
                                   for wx.ColourSelect
================================  =================================================================

GSAS-II-provided Dialog (full window) routines:

.. tabularcolumns:: |l|p{4in}|

================================  =================================================================
Class or function name             Description
================================  =================================================================
:class:`gpxFileSelector`           File browser dialog for opening existing .gpx files
:func:`askSaveFile`                Get a file name from user
:func:`askSaveDirectory`           Get a directory name from user
:class:`MultiColumnSelection`      A dialog that builds a multicolumn table, word wrapping
                                   is used for the 2nd, 3rd,... columns. 
:func:`MultiColMultiSelDlg`        Calls a dialog that allows multiple selection from a
                                   multicolumn table
:class:`MultiDataDialog`           Dialog to obtain multiple data values from user, 
                                   with optional range validation; items can be float, str or bool
:class:`MultiIntegerDialog`        Dialog to obtain multiple integer values from user, 
                                   with a description for each value and optional
                                   defaults.
:class:`MultiStringDialog`         Dialog to obtain multiple string values from user, 
                                   with a description for each value and optional
:class:`DisAglDialog`              Distance/Angle Controls input dialog. 
:class:`FlagSetDialog`             Dialog that provides a table of items along with a
                                   checkbox for each. 
                                   defaults.
:class:`G2ColumnIDDialog`          A dialog for matching column data to desired items; some
                                   columns may be ignored.
:class:`G2HistoDataDialog`         A dialog for global edits to histogram data globally
:class:`G2MultiChoiceDialog`       Dialog similar to wx.MultiChoiceDialog, but provides
                                   a filter to select choices and buttons to make selection
                                   of multiple items more simple.
:class:`G2SingleChoiceDialog`      Dialog similar to wx.SingleChoiceDialog, but provides
                                   a filter to help search through choices.
:class:`SGMagSpinBox`              Special version of wx.MessageBox that displays
                                   magnetic spin text
:class:`SGMessageBox`              Special version of wx.MessageBox that displays space group 
                                   & super space group text in two blocks
:class:`SingleFloatDialog`         Dialog to obtain a single float value from user, with
                                   optional range validation.
:class:`SingleIntDialog`           Dialog to obtain a single integer value from user,
                                   with optional range validation.
:class:`SingleStringDialog`        Dialog to obtain a single string value from user, 
                                   with optional an optional default value.
:func:`G2MessageBox`               Displays text typically used for errors or warnings. 
:func:`ShowScrolledInfo`           Dialog to display longer text where scrolling is
                                   possibly needed
:func:`GetItemOrder`               Dialog for ordering items into columns
:func:`GetImportFile`              Gets one ore more file from the appropriate import
                                   directory, which can be overridden. Arguments follow those
                                   of :func:`wx.FileDialog`
:class:`ScrolledMultiEditor`       Dialog for editing many dict- or list-contained items.
                                   with validation. Results are placed in dict or list.
:func:`CallScrolledMultiEditor`    Routine for editing many dict- or list-contained items.
                                   using the :class:`ScrolledMultiEditor` dialog
:func:`G2ScrolledGrid`             Displays a multicolumn table of information with 
                                   possible scroll bars
:func:`ShowScrolledColText`        Displays tabular text with scrolling where needed
:func:`ItemSelector`               Select a single item or multiple items from list of choices.
                                   Creates and then destroys a wx.Dialog and returns the
                                   selections(s).
:func:`SelectEdit1Var`             Select a variable from a list, then edit it and select
                                   histograms to copy it to.
:func:`BlockSelector`              Select a single block for instrument parameters
:func:`MultipleBlockSelector`      Select one or more blocks of data, used for 
                                   CIF powder histogram imports only
:func:`MultipleChoicesSelector`    Dialog for displaying fairly complex choices, used for 
                                   CIF powder histogram imports only
:func:`PhaseSelector`              Select a phase from a list (used for phase importers)
================================  =================================================================

Miscellaneous GUI support routines:

.. tabularcolumns:: |l|p{4in}|

================================  =================================================================
Function name                      Description
================================  =================================================================
:func:`Define_wxId`                Create a unique wx.Id symbol in _initMenus in :mod:`GSASIIdataGUI`.
                                   Such symbols are needed when the menu item is defined in a 
                                   different location from the wx.Bind that links the menu item 
                                   to a function. This function allows all the menu Ids to be
                                   defined as the menus are created in one place and then can be 
                                   used in Bind elsewhere in the code.
:func:`openInNewTerm`              opens a Python routine (usually G2.py) in a
                                   new terminal window (works on all platforms)				   
================================  =================================================================

Other miscellaneous non-GUI routines that may be of use for GUI-related actions:

.. tabularcolumns:: |l|p{4in}|

================================  =================================================================
Function name                      Description
================================  =================================================================
:func:`StripIndents`               Regularizes the intentation from a string with multiple
                                   newline characters by removing spaces at the beginning
                                   of each line.
:func:`StripUnicode`               Removes unicode characters from strings 
:func:`GetImportPath`              Determines the default location to use for importing files.
                                   Tries sequentially :attr:`G2frame.TutorialImportDir`,
                                   config var ``Import_directory`` and
                                   :attr:`G2frame.LastImportDir`.
:func:`GetExportPath`              Determines the default location to use for writing files.
                                   Tries sequentially :attr:`G2frame.LastExportDir` and
                                   :attr:`G2frame.LastGPXdir`
================================  =================================================================


GSASIIctrlGUI Classes & Routines
---------------------------------------------

.. automodule:: GSASIIctrlGUI
    :members: 

---------------------------------------------
*GSASIImiscGUI: Misc I/O routines*
---------------------------------------------

Module with miscellaneous routines for input and output. Many
are GUI routines to interact with user.

Includes support for image reading.

GSASIImiscGUI Classes & Routines
---------------------------------------------
       
.. automodule:: GSASIImiscGUI
    :members: 

---------------------------------------------
*gltext: draw OpenGL text*
---------------------------------------------

.. automodule:: gltext
    :members: 

