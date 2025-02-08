*GSASII: GSAS-II GUI*
 ===========================

Script G2.py 
----------------------------------------

File `G2.py` can be used to start the GSAS-II graphical user 
interface (GUI), particularly when GSAS-II has been installed
in a location outside of Python and thus  requires changing the Python
path. When GSAS-II is installed in a location that is on the
default Python path, this command is
sufficient to start the GSAS-II GUI:

`python -c "from GSASII.GSASIIGUI import main; main()"`

The `G2.py` script checks to see if GSAS-II is on the path. If not,
the directory where the `G2.py` file is located is placed into the
Python path. At this point the func:`GSASIIGUI.main` routine is
called to start the GSAS-II GUI.

Module GSASIIGUI.py
----------------------------------------

The `GSASIIGUI.py` script imports GSASIIpath, which does some minor initialization
and then (before any wxPython calls can be made) creates a wx.App application. 
At this point :func:`GSASIIpath.SetBinaryPath` is called to establish
the directory where GSAS-II binaries are found. If the binaries 
are not installed or are incompatible with the OS/Python packages, 
the user is asked if they should be updated from the subversion site. 
The wxPython app is then passed to :func:`GSASIIdataGUI.GSASIImain`, 
which creates the GSAS-II GUI and finally the event loop is started.

Keyboard Menu Shortcuts
----------------------------------------

Shortcuts for commonly-used menu commands are created by adding a 
menu command with a "\\tCtrl+" addition such as::

        item = parent.Append(wx.ID_ANY,'&Refine\tCtrl+R','Perform a refinement')

This will allow the above menu command to be executed with a "Control-R" 
keyboard command (on MacOS this will be "Command+R" rather than "Control-R") as well as using the menu to access that action. The following table lists the 
keyboard letters/numbers that have GSAS-II assigned actions.
are system assigned. Note that there are also plotting keyboard commands that are 
implemented in :mod:`GSASIIplot`. 
These can be discovered from the "K" button on the plot menu bar, as they 
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

GSAS-II contents
----------------------------------------


.. automodule:: GSASIIGUI
    :members: 
    :private-members:
    :special-members:
