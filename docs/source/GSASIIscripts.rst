======================
*GSAS-II Misc Scripts*
======================

*testDeriv: Check derivative computation*
=========================================

Use this to check derivatives used in structure least squares
refinement against numerical values computed in this script.

.. automodule:: GSASII.testDeriv
    :members: 

*GSASIItestplot: Plotting for testDeriv*
========================================

Plotting module used for script testDeriv.

.. automodule:: GSASII.GSASIItestplot
    :members: 

*scanCCD: reduce data from scanning CCD*
========================================

Quickly prototyped routine for reduction of data from detector described in
B.H. Toby, T.J. Madden, M.R. Suchomel, J.D. Baldwin, and R.B. Von Dreele,
"A Scanning CCD Detector for Powder Diffraction Measurements".
Journal of Applied Crystallography. 46(4): p. 1058-63 (2013). This is
no longer being updated. 

.. automodule:: GSASII.scanCCD
    :members: 

*makeMacApp: Create MacOS Applet*
=================================

This script creates an AppleScript app bundle to launch GSAS-II. It is
called by bootstrap.py during the GSAS-II installation process. It
creates a "copy" of Python that is able to run wx.Python programs and
names this version of Python as GSAS-II so that items in the menus are
named correctly. 

.. automodule:: GSASII.install.makeMacApp
    :members: 

*makeBat: Create GSAS-II Batch File*
====================================

This script performs Windows specific installation steps to allow for
easy launching of GSAS-II. It is
called by bootstrap.py during the GSAS-II installation process.

.. automodule:: GSASII.install.makeBat
    :members: 

*makeLinux: Create Linux Shortcuts*
===================================

This script performs Linux specific installation steps that
allowscreates files allowing 
GSAS-II to be launched from a desktop icon or desktop manager menu.
Not all desktop managers will recognize these files. 
It is called by bootstrap.py during the GSAS-II installation process.

.. automodule:: GSASII.install.makeLinux
    :members: 

*makeVarTbl: Make Table of Variable Names*
============================================

This creates a table of variable names from the definitions supplied
in :func:`GSASIIobj.CompileVarDesc`. This table is used in the
Sphinx documentation as the :ref:`GSAS-II Variable Names table <VarNames_table>`.
This is run as part of the Sphinx build from inside ``docs/source/conf.py``.
 
.. automodule:: GSASII.install.makeVarTbl
    :members: 
       
       
*testSytSym: Test Site Symmetry*
========================================

A GUI program for testing the site symmetry generation routines. 
       
.. automodule:: GSASII.testSytSym
    :members: 

*testSSymbols: Test Superspace Group Symbols*
===============================================

A GUI program for testing the 3+1 superspace group symmetry generation routines. 
       
.. automodule:: GSASII.testSSymbols
    :members: 

*Self-test Modules*
===================================

A set of scripts that can be run to test a series of self-tests in GSAS-II. 

.. automodule:: tests.test_diffax
    :members: 

.. automodule:: tests.test_elm
    :members: 

.. automodule:: tests.test_image
    :members: 

.. automodule:: tests.test_kvec
    :members: 

.. automodule:: tests.test_lattice
    :members: 

.. automodule:: tests.test_nistlat
    :members: 

.. automodule:: tests.test_scriptref
    :members: 

.. automodule:: tests.test_spg
    :members: 

.. automodule:: tests.test_tofref
    :members: 

*Other scripts*
========================================

A few scripts are also placed in the GSAS-II auxiliary repositories 

gitstrap.py
----------------------------------------------------------

    ``gitstrap.py`` in ``GSASII-buildtools/install/`` is 
    used to install the GSAS-II package, including the appropriate
    binary files. May be used directly to install GSAS-II from inside
    Python in an appropriately configured Python installation, or
    is also used to obtain or update the GSAS-II files in a conda
    installation. 

gitcompile.py
----------------------------------------------------------

    ``gitcompile.py`` in ``GSASII-buildtools/install/`` is 
    used to install the GSAS-II package, but also compiles 
    the binary files. May be used directly to install GSAS-II from inside
    Python in an appropriately configured Python installation. 


makeGitTutorial.py
----------------------------------------------------------

   ``makeGitTutorial.py`` in 
   ``GSASII-tutorials/scripts/``
   provides a script to creates the HTML page
   (``GSASII/help/Tutorials.html``) that lists all the tutorials defined in
   variable :data:`GSASIIctrlGUI.tutorialIndex`. Run this after adding
   new tutorials to that catalog.


tag-version.py
----------------------------------------------------------

   ``tag-version.py``   in  ``GSASII/install/``
   creates a new numerical tag number (advancing from ``5898`` to
   ``5899``) for the most recent git check in and records that in
   the ``git_version.py`` file. This also advances the minor version number
   for the GSAS-II version number (from ``5.x.y`` to ``5.(x+1).0``).
   Use this when there is a significant change to GSAS-II
   functionality. 


incr-mini-version.py
----------------------------------------------------------

   ``incr-mini-version.py``   in  ``GSASII/install/``
   creates a new numerical tag number (advancing from ``5898`` to
   ``5899``) for the most recent git check in and records that in
   the ``git_version.py`` file. This also advances the "mini" version number
   for the GSAS-II version number (from ``5.x.y`` to ``5.x.(y+1)``).
   Use this to note a minor by noteworthy change to GSAS-II
   functionality, such as a bug fix where users should be aware that
   something has changed. 

macStartScript.py 
----------------------------------------------------------
   ``macStartScript.py``   in  ``GSASII/install/``
   creates a MacOS applet to start GSAS-II similar to what is in
   ``makeMacApp.py`` (where the app is taken from a tar file).
   This is not in regular use, as it seems to have some permissions
   problems, but may be needed to update ``makeMacApp.py``.
   
   
