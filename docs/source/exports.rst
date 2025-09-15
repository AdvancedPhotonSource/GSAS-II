*GSAS-II Export Modules*
====================================

Many of the data files written by GSAS-II for use by other software is
written using a layer of routines called exporters. Exporters usually require quite
simple code, so they can be written quickly to allow GSAS-II to
provide output in new data formats. 
The interface to the exporters 
is self-configuring, so all supplied exporters are available once the
exporter is added to the code base. This allows GSAS-II to be quite
flexible in adapting to use many data formats without need for
extensive coding. 

Exports are implemented by deriving a class from 
:class:`~GSASII.GSASIIfiles.ExportBaseclass` in module
:mod:`~GSASII.GSASIIfiles`.      
Export routines commonly access the GUI to determine information on
the contents of the file(s) to be written, but for powder diffraction
and phase data export, a routine named ``Writer()`` may be defined. If this is
present, the exporter can be used from: mod:`~GSASII.GSASIIscriptable`
without GUI access. Note that the arguments for the``Writer()`` method
include a histogram tree name as well as a file name to
be written.

A file containing one or more export routine can be placed
either in the ``GSASII/exports`` directory (which requires
modification of the ``__init__.py`` file or the file can be placed in
the ``~/.GSASII/exports`` directory.
(Note that ``~`` here is translated to the 
user's home directory; for Windows this is usually taken from the
USERPROFILE setting or a combination of HOMEPATH and HOMEDRIVE,
so this directory will usually have form
``C:\\Users\\YourUsername\\.GSASII\\exports``.
The next time GSAS-II is started,
the file will be loaded with all the other GSAS-II files and
the new data format will appear in the appropriate exporter menu. 

.. _export_routines: 

======================================
 Writing an Exporter Routine
======================================

When writing a exporter routine, one should create a new class derived
from class `~GSASII.GSASIIfiles.ExportBaseclass`.
The name of the class is arbitrary, but if more than one class is
placed in file, each class must have a different name. The same name
can be repeated if it is in different files. 
As described below, this class will implement
an ``__init__()`` and an ``Exporter()`` method, and many will supply a 
``Writer()`` method, too. The purpose of each of these
routines is described below. The easiest way to craft a new exporter
will be to use the other exporters of the same data type as a model
for where to find the data values that will be written, but the documentation 
for the parent class (`~GSASII.GSASIIfiles.ExportBaseclass`) provides
useful information on support routines that pull information from the
GSAS-II data structures into the exporter. 

__init__()
--------------
 
The ``__init__`` method will follow standard boilerplate largely independent
of the data type: 

.. code-block:: python

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            G2frame=G2frame,
            formatName = 'Format name for menu',
            extension='.extn',
            longFormatName = 'Longer more detailed format name for status line'
           )

The first line in the ``__init__`` method calls the parent class 
``__init__`` method with the following parameters: 

  * ``G2frame``: a reference to the main GSAS-II GUI window or None
    when run scripted.
  * ``formatName``: a string to be used in the menu. Should be short. 
  * ``extension``: a string to be used in the file name. All files
    produced by the exporter will have this extension.
  * ``longFormatName``: a longer string to be used to describe the
    format in help. 

In addition, one instance variables must be defined:

.. code-block:: python

        self.exporttype = ['phase']

The value for ``self.exporttype`` determines the type of export that will be performed
('project', 'phase', 'single', 'powder', 'image', 'map', 'sasd', 'refd' or (someday)
'pdf') and the menu where the exporter will be placed.  Note that
'project' exports are those that include all data from a
.gpx file (all phases, histograms, etc.)

Another item is optional. 

.. code-block:: python

        self.multiple = True

The value specified for ``self.multiple``
determines if only a single phase, data set, etc. can be exported at a
time (when False) or when True, a file can be produced with multiple
histograms, phases, etc.

Exporter()
--------------

The class must supply a ``Exporter`` method that will write a
file. Depending on the settings for ``self.exporttype`` and
``self.multiple`` and the contents of the Data Tree, will dictate 
what dialogs will be presented to the user to select what will be
written. 

Writer()
--------------

For powder and phase exports, if the class supplies a ``Writer()``
module, then the export format
will be available for scripted output (with
:mod:`GSASIIscriptable`). These modules are supplied a histogram name
or a phase name and should not attempt to access the GUI. It is not
required that that this method be supplied, but usually it is not hard to do, 
unless information from the user is required.

Note that for phase exports the ``Writer()`` should be declared as follows:

.. code-block:: python

       def Writer(self,hist,phasenam,mode='w'):

while for histogram exports the ``Writer()`` should be declared as follows:

.. code-block:: python

     def Writer(self,hist,filename=None,mode='w'):

-------------------------------------------
*Module G2export_examples: Examples*
-------------------------------------------

 .. py:currentmodule:: GSASII.exports.G2export_examples  
 
Code to demonstrate how GSAS-II data export routines are created. The
classes defined here, :class:`ExportPhaseText`, 
:class:`ExportSingleText`, :class:`ExportPowderReflText`, 
and :class:`ExportPowderText` each demonstrate a different type
of export. Also see
:class:`~GSASII.exports.G2export_map.ExportMapASCII` in :mod:`~GSASII.exports.G2export_map` for an example of a map export.

G2export_examples Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_examples
    :members: 
    :synopsis: Demonstrates sample code that exports a phase or dataset to
      a text file.

-------------------------------------------
*Module G2export_csv: Spreadsheet export*
-------------------------------------------

 .. py:currentmodule:: GSASII.exports.G2export_csv

Code to create .csv (comma-separated variable) files for
GSAS-II data export to a spreadsheet program, etc. Defines a number of
.csv exports:

 * :class:`ExportPhaseCSV`: phases
 * :class:`ExportPowderCSV`: powder data, includes instrument parameters as well as obs & calc patterns, etc.
 * :class:`ExportMultiPowderCSV`: multiple powder datasets in a single spreadsheet
 * :class:`ExportPowderReflCSV`: reflections from a powder fit
 * :class:`ExportSASDCSV`: small angle data set
 * :class:`ExportREFDCSV`: reflectometry data set
 * :class:`ExportSingleCSV`: single crystal reflection data
 * :class:`ExportStrainCSV`: reflectometry datasets

G2export_csv Classes and Routines
-------------------------------------------
      
.. automodule:: GSASII.exports.G2export_csv
    :members: 
    :synopsis: Exports a phase or dataset to a spreadsheet via a 
       comma-separated-variable (csv) format file.

--------------------------------------------
*Module G2export_PDB: Macromolecular export*
--------------------------------------------

.. py:currentmodule:: GSASII.exports.G2export_PDB  

Code to export a phase into the venerated/obsolete (pick one)
ASCII PDB format. Also defines exporter :class:`ExportPhaseCartXYZ`
which writes atom positions in orthogonal coordinates for a phase.

G2export_PDB Classes and Routines
-------------------------------------------
       
.. automodule:: GSASII.exports.G2export_PDB
    :members: 
    :synopsis: Cartesian coordinate export, including PDB format

------------------------------------------------------
*Module G2export_image: 2D Image data export*
------------------------------------------------------

Demonstrates how an image is retrieved and written. Uses
a SciPy routine to write a PNG format file. 
	       
G2export_image Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_image
    :members: 
    :synopsis: Exports images

-------------------------------------------
*Module G2export_map: Map export*
-------------------------------------------

.. py:currentmodule:: GSASII.exports.G2export_map

Code to write Fourier/Charge-Flip atomic density maps out in formats that
can be read by external programs. At present a GSAS format
that is supported by FOX and DrawXTL 
(:class:`ExportMapASCII`) and the CCP4 format that
is used by COOT (:class:`ExportMapCCP4`) are implemented.

G2export_map Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_map
    :members: 
    :synopsis: Export Fourier and charge-flip atomic density maps
	      
-------------------------------------------
*Module G2export_shelx: Examples*
-------------------------------------------

Code to export coordinates in the SHELX .ins format
(as best as we can make sense of it).

G2export_shelx Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_shelx 
    :members: 
    :synopsis: Export a phase in Shelx format

------------------------------------------------------
*Module G2export_CIF: CIF Exports*
------------------------------------------------------

.. py:currentmodule:: GSASII.exports.G2export_CIF
		      
This implements a complex set of CIF (Crystallographic Information
Framework) exporters. The base class, :class:`ExportCIF`, implement a
variety of export capabilities,
where extra parameters for :meth:`ExportCIF.MasterExporter` determine if a project,
single phase or data set are written. The subclasses of
:class:`ExportCIF`, as listed below, supply these different parameters
when calling that method. 

 *  :class:`ExportProjectCIF`: writes an entire project in a complete
    CIF intended for submission as a publication,
 *  :class:`ExportPhaseCIF`: writes a single phase in CIF
 *  :class:`ExportPwdrCIF`: writes a one powder diffraction dataset CIF
 *  :class:`ExportHKLFCIF`: writes a single crystal dataset


G2export_CIF Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_CIF
    :members: 
    :synopsis: Export a project in CIF format
	       
-------------------------------------------------
*Module G2export_pwdr: Export powder input files*
-------------------------------------------------

Creates files used by GSAS (FXYE) & TOPAS (XYE) as input

G2export_pwdr Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_pwdr
    :members: 
    :synopsis: Export powder data in GSAS and Topas formats

-------------------------------------------
*Module G2export_FIT2D: Fit2D "Chi" export*
-------------------------------------------

Code to create .chi (Fit2D like) files for GSAS-II powder data export 

G2export_FIT2d Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_FIT2D 
    :members: 
    :synopsis: Export powder data in Fit2D (.chi) format
       
------------------------------------------------------
*Module G2export_JSON: ASCII .gpx Export*
------------------------------------------------------

.. py:currentmodule:: GSASII.exports.G2export_JSON

This implements a fairly simple exporter, :class:`ExportJSON`, that can export the 
contents of an entire project as a sort-of human readable (JSON) ASCII file.
This provides a way to see the contents of a GSAS-II project
file. This does not provide a mechanism to change the contents of a .gpx file,
since there are no provisions to read this file back into GSAS-II, as
the  likelihood of breaking a data structure is too high. 
If you want to change the contents of a .gpx file, use :mod:`GSASIIscriptable` 
where you can access the native Python data structures and change things, 
with a good chance of getting things to work.

G2export_JSON Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_JSON
    :members:

------------------------------------------------------
*Module G2export_Bracket: ASCII .gpx Export*
------------------------------------------------------

.. py:currentmodule:: GSASII.exports.G2export_Bracket
		      
This provides to methods for tabulating GSAS-II parameters from a
project for use in manuscript preparation into an ASCII .csv
(spreadsheet) file.
The exporter, :class:`Exportbracket`, creates a text file with
standard uncertainties for values in crystallographic (e.g. "bracket")
notation:
*i.e.*: ``1.234(5)``, which
indicates a value of ``1.234`` with a standard uncertainty of ``0.005``. A
second method, :class:`Export3col`, provides the standard uncertainties
as a separate column. 

This module initially written by Conrad Gillard. For any enquiries please contact conrad.gillard@gmail.com.
       
G2export_Bracket Classes and Routines
-------------------------------------------

.. automodule:: GSASII.exports.G2export_Bracket
    :members: 
