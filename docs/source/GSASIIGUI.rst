.. _GSASIIGUI_chapter:

*GSAS-II GUI Components*
========================

These modules are used to create different parts of the GSAS-II
graphical user interface (GUI).

--------------------------------------
*G2.py: Routine to start GSAS-II*
--------------------------------------

`G2.py` is a short Python file that can be used to start the GSAS-II
GUI, particularly when GSAS-II has been installed
in a location outside of Python and thus requires changing the Python
path.


Note that when GSAS-II is installed in a location that is on the
default Python path (or the current working directory is the one that
contains the `GSASII` directory) either of these two commands is
sufficient to start the GSAS-II GUI:

``python -c "from GSASII.GSASIIGUI import main; main()"``

``python -m GSASII``

.. automodule:: GSASII.G2
        :members:

--------------------------------------
*GSASIIGUI: Main GUI for GSAS-II*
--------------------------------------

The `GSASIIGUI.py` module imports GSASIIpath, which does some minor initialization
and then (before any wxPython calls can be made) creates a wx.App application. 
At this point :func:`GSASIIpath.SetBinaryPath` is called to establish
the directory where GSAS-II binaries are found. If the binaries 
are not installed or are incompatible with the OS/Python packages, 
the user is asked if they should be updated from the subversion site. 
The wxPython app is then passed to :func:`GSASIIdataGUI.GSASIImain`, 
which creates the GSAS-II GUI and finally the event loop is started.

.. automodule:: GSASII.GSASIIGUI
    :members: 
    :private-members:
    :special-members:


----------------------------------------
*GSASIIdataGUI: overall GUI components*
----------------------------------------

This module that defines GUI routines and classes for the main GUI Frame (window)
including the routines that define the GSAS-II tree panel and much of the
data editing panel. Most menus are defined here but are not invoked
until later to reduce the GUI startup time. 

Routine :func:`GSASIIdataGUI.SelectDataTreeItem` is called
to respond to selections in the data tree and invoke the routines that
display the appropriate information in the Data Window. 
Two of the simpler data tree items (Notebook and Controls) are
implemented from here, but most are defined in other modules. 

GSASIIdataGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIdataGUI
        :members: 

----------------------------------------
*GSASIIseqGUI: Sequential Results GUI*
----------------------------------------

Module that defines GUI routines and classes for the various
sequential result GUI Frames (window).
Also defines GUI routines for Cluster Analysis results. 
   
Note that there are seven types of sequential results that GSAS-II can produce 
and all are displayed/analyzed with the code in this module. They vary by title so that 
a project can hold one result of each type without a naming collision:

 * Rietveld (Title: Sequential results)
 * PDF (Title: Sequential PDFfit2 results)
 * Peak fit (Title: Sequential peak fit results)
 * Small angle (Title: Sequential SASD fit results) 
 * Reflectometry (Title: Sequential REFD results)
 * Image (strain) (Title: Sequential strain fit results)
 * Image (calibration) (Title: Sequential image calibration results)


       
GSASIIseqGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIseqGUI
        :members: 

---------------------------
*Phase GUI files*
---------------------------

Routines create the GUI for display of phase information
in the data display window when a phase is selected.
Phase information is stored in one or more
:ref:`Phase Tree Item <Phase_table>` objects.
Note that there are functions
that respond to some tabs in the phase GUI in other modules
(such as GSASIIddata). This is large enough that it has been split
into three files. 

GSASIIphsGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIphsGUI
    :members: 

GSASIIphsGUI2 Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIphsGUI2
    :members: 

GSASIIrmcGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIrmcGUI
    :members:

--------------------------------------------
*GSASIIddataGUI: Phase Diffraction Data GUI*
--------------------------------------------

Module to create the GUI for display of HAP items (where there is
an entry for each histogram & phase). This is shown when the
Phase "Data" tab is selected or may appear as if in a separate
data tree item (see SeparateHistPhaseTreeItem in config.py).

GSASIIddataGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIddataGUI
    :members: 

-------------------------------------------------------
*GSASIIElemGUI: GUI to select and delete element lists*
-------------------------------------------------------

Module to select elements from a periodic table and
to delete an element from a list of selected elements.

GSASIIElemGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIElemGUI
    :members: 

------------------------------------------
*GSASIIconstrGUI: Constraint GUI routines*
------------------------------------------

GUI routines to define constraints and rigid bodies.


GSASIIconstrGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIconstrGUI
    :members: 

----------------------------------------
*GSASIIrestrGUI: Restraint GUI routines*
----------------------------------------

GUI Routines used to define restraints.


GSASIIrestrGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIrestrGUI
    :members: 

-------------------------
*GSASIIimgGUI: Image GUI*
-------------------------

GUI Routines used to control image display and processing


GSASIIimgGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIimgGUI
    :members: 

-------------------------------------------
*GSASIIpwdGUI: Powder Pattern GUI routines*
-------------------------------------------

Used to define GUI controls for the routines that interact
with the powder histogram (PWDR) data tree items.


GSASIIpwdGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIpwdGUI
    :members: 

-------------------------------------
*GSASIIexprGUI: Expression Handling*
-------------------------------------

This module defines a class for defining an expression in terms of values
in a parameter dictionary via a wx.Dialog. The dialog creates a
:class:`GSASIIexprGUI.GSASII.ExpressionObj`
which is used to evaluate the expression against a supplied parameter dictionary.

The expression is parsed to find variables used in the expression and then
the user is asked to assign parameters from the dictionary to each variable.

Default expressions are read from file DefaultExpressions.txt using
:func:`GSASIIpath.LoadConfigFile`.


GSASIIexprGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIexprGUI
    :members: 

-------------------------------------------------
*GSASIIfpaGUI: Fundamental Parameters Routines*
-------------------------------------------------

This module contains GUI routines to accept Fundamental Parameters 
Approach (FPA) input used to run the NIST XRD Fundamental 
Parameters Code, computes a set of peaks with that code and fits
profile terms to the peaks. 
Also allows for plotting the convolutors generated by that code. 

GSASIIfpaGUI Classes & Routines
---------------------------------------

.. automodule:: GSASII.GSASIIfpaGUI
    :members: 

-------------------------------------------------
*fprime: compute f' & f"*
-------------------------------------------------

This module contains GUI routines to calculate and plot real and resonant X-ray
scattering factors to 250keV, based on Fortran program of Cromer &
Liberman corrected for Kissel & Pratt energy term; Jensen term not included
(D. T. Cromer and D. A. Liberman, Acta Cryst. (1981). A37, 267-268.)


fprime Classes & Routines
---------------------------------------

.. automodule:: GSASII.fprime
    :members: 


-------------------------------------------------
*Absorb: Compute X-ray Absorption*
-------------------------------------------------

This module contains GUI routines to calculate X-ray absorption
factors to 250keV for cylindrical powder samples with a specified
chemical composition and density. Uses same approach
as in :mod:`fprime`
(D. T. Cromer and D. A. Liberman, Acta Cryst. (1981). A37, 267-268.)

Absorb Classes & Routines
---------------------------------------

.. automodule:: GSASII.Absorb
    :members: 


-------------------------------------------------
*PlotXNFF: Show Form Factors*
-------------------------------------------------

.. automodule:: GSASII.PlotXNFF
    :members: 
