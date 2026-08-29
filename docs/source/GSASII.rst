*GSASII: GSAS-II GUI*
===========================

The GSAS-II GUI is normally started using the `G2.py` file. This will
"hack" the PythonPath, if needed,  by adding the diectory containing the `GSASII`
directory to the path, allowing GSASII to be run without installing it
inside Python. The  `G2.py` invokes the `main()` routine inside
`GSASIIGUI.py` which actually starts the GUI.

When GSAS-II has been installed inside Python, either of these two commands is
sufficient to start the GSAS-II GUI:

``python -c "from GSASII.GSASIIGUI import main; main()"``

``python -m GSASII``

For documentation on the main GUI routines see the
:ref:`GSAS-II GUI Components <GSASIIGUI_chapter>` chapter.

Keyboard Menu Shortcuts
----------------------------------------

Shortcuts for commonly-used menu commands are created by adding a 
menu command with a "\\tCtrl+" addition such as::

        item = parent.Append(wx.ID_ANY,'&Refine\tCtrl+R','Perform a refinement')

This will allow the above menu command to be executed with a "Control-R" 
keyboard command (on MacOS this will be "Command+R" rather than "Control-R") as well as using the menu to access that action. The following table lists the 
keyboard letters/numbers that have GSAS-II assigned actions.
are system assigned.
Note that there are also plotting keyboard shortcut commands that are 
implemented in :mod:`GSASIIplot`. 
These can be discovered from the "K" button on the plot menu bar, as they 
vary depending on the type of plot that is shown.

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
