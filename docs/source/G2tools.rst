*GSAS-II Independent Tools*
================================

The modules here are used for independent programs to be used as
tools within the GSAS-II package and run
independently of the main GSAS-II program.

 * :ref:`GSASIIIntPDFtool<GSASIIautoint>`:  Parallelized auto-integration/PDF program 
 * :ref:`G2compare<G2compare>`: Project Comparison program

Both are under development.

.. _GSASIIautoint:

*GSASIIIntPDFtool: autointegration routines*
---------------------------------------------

An auto-integration program based on GSAS-II but with
a minimal GUI and no visualization that runs independently from
the main GSAS-II program . This is intended to implement
significant levels of parallelization and require less of a memory footprint. 

.. automodule:: GSASII.GSASIIIntPDFtool
    :members:

.. _G2compare:

*G2compare: Tool for project comparison*
---------------------------------------------

This is intended to read in multiple GSAS-II projects and provide
graphics, tables of information and so on. Not much of this has been
written at present. 

.. automodule:: GSASII.G2compare
    :members: 
