--------------------------------------------------------
*GSAS-II Importer Modules*
--------------------------------------------------------

Almost all the types of data files that GSAS-II can accept are read
using a layer of routines called importers. (The one exception is the
reading of powder peak positions.)  Importers usually require quite
simple code, so they can be written quickly for new data formats. 
The interface to the importers 
is self-configuring, so all supplied importers are available once the
importer is added to the code base. This allows GSAS-II to be quite
flexible in adapting to use many data formats without need for
extensive coding. 

Imports are implemented by deriving a class from 
:class:`GSASIIobj.ImportPhase`, :class:`GSASIIobj.ImportStructFactor`,
:class:`GSASIIobj.ImportPowderData` ,
:class:`GSASIIobj.ImportSmallAngleData`,  
:class:`GSASIIobj.ImportReflectometryData`,  
:class:`GSASIIobj.ImportPDFData`,  
or :class:`GSASIIobj.ImportImage`. These classes are in turn 
derived from :class:`GSASIIobj.ImportBaseclass`.
Module file names (`G2phase_`, `G2pwd_` and `G2sfact_`, etc.) are used to
determine which type of data will be read and which menu an importer
routine should be placed into. (N.B. in retrospect this
naming was an unnecessary choice; importer types could have been determined
from the base class as is done for exporters.)

To implement the import of 
a phase, a single crystal or a powder dataset, etc., create a file
named with the appropriate file name prefix:

* 'img' (images), 
* 'pdf' (pair distribution function), 
* 'phase' (cell/coordinates), 
* 'pwd' (powder diffraction), 
* 'rfd' (reflectivity), 
* 'sad' (small-angle scattering) or 
* 'sfact' (single crystal).

The importer file can be placed
either in the ``GSASII/imports`` directory (which requires
modification of the ``__init__.py`` and the ``meson.build`` files
or the importer file can be placed in
the ``~/.GSASII/imports`` directory.
(Note that ``~`` here is translated to the 
user's home directory; for Windows this is usually taken from the
USERPROFILE setting or a combination of HOMEPATH and HOMEDRIVE,
so this directory will usually have form
``C:\\Users\\YourUsername\\.GSASII\\imports``.
The next time GSAS-II is started,
the file will be loaded with all the other GSAS-II files and
the new data format(s) will appear in the appropriate importer menu. 
The importer file may contain a
single importer class or several. 

Importers are documented below, separated by type. Importers tend to
be fairly simple files, where many are in the range of 50-100 lines,
and where more than half of those lines are directly copied from other
importers without any changes. Details on this are given in the
:ref:`Writing a Importer Routine<import_routines>` section,
immediately below.

.. _import_routines: 

======================================
 Writing an Importer Routine
======================================

When writing a importer routine, one should create a new class derived
from
:class:`GSASIIobj.ImportPhase`, :class:`GSASIIobj.ImportStructFactor`,
:class:`GSASIIobj.ImportPowderData` ,
:class:`GSASIIobj.ImportSmallAngleData`,  
:class:`GSASIIobj.ImportReflectometryData`,  
:class:`GSASIIobj.ImportPDFData`,  
or :class:`GSASIIobj.ImportImage`.
The name of the class is arbitrary, but if more than one class is
placed in file, each class must have a different name. The same name
can be repeated if it is in different files. 
As described below, 
to implement an importer class, you must implement
an ``__init__()`` and a ``Reader()`` method, and many will supply a 
``ContentsValidator()`` method, too. The purpose of each of these
routines is described below. The easiest way to craft a new importer
will be to use the other importers of the same data type as a model
for what values should be set inside each routine, with most of the
work needed to create a ``Reader()`` routine. The documentation 
for the parent class may also have some useful information. See the appropriate class 
for details on what values each type of ``Reader()`` should
set. 

__init__()
--------------
 
The ``__init__`` method will follow standard boilerplate largely independent
of the data type: 

.. code-block:: python

    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.ext1','ext2'),
            strictExtension=True,
            formatName = 'example image',
            longFormatName = 'A longer description that this is an example image format'
            )

The first line in the ``__init__`` method calls the parent class 
``__init__`` method with the following parameters: 

  * ``extensionlist``: a list of extensions that may be used for this type of file.
  * ``strictExtension``: Should be True if only files with extensions in
    ``extensionlist`` are allowed; False if all file types should be offered
    in the file browser. Also if False, the importer class will be
    used on all files when "guess from format" is tried, though 
    readers with matching extensions will be tried first. 
    It is a very good idea to supply  a :ref:`ContentsValidator <ContentsValidator>`
    method when ``strictExtension`` is False.
  * ``formatName``: a string to be used in the menu. Should be short. 
  * ``longFormatName``: a longer string to be used to describe the
    format in help. 

Note that if an importer detects a condition which prevents its use,
for example because a required Python package is not present, it can
set the value of ``self.UseReader`` to False. Another possible use for
this would be an importer that requires a network connection to a
remote site. Setting ``self.UseReader`` to False must be done in the 
``__init__`` method and will prevent the
importer from being accessed or included in the appropriate GUI menu. 

Reader()
--------------

The class must supply a ``Reader`` method that actually performs the
reading. All readers must have at a minimum these arguments::

    def Reader(self, filename, filepointer, ParentFrame, **unused):

where the arguments have the following uses: 

 * ``filename``: a string with the name of the file being read
 * ``filepointer``: a file object (created by :func:`open`) that accesses
   the file and is points to the beginning of the file when Reader is
   called. 
 * ``ParentFrame``: a reference to the main GSAS-II (tree) windows, for
   the unusual ``Reader`` routines that will create GUI windows to ask
   questions. For use with scripting, the Reader should do something
   reasonable such as assume a default if ``ParentFrame`` is None,
   which indicates that the GUI should not be accessed.
   
In addition, the following keyword parameters are defined that ``Reader``
routines may optionally use:

 * ``buffer``: a dict that can be used to retain information between repeated calls of the routine
 * ``blocknum``: counts the number of times that a reader is called, to
   be used with files that contain more than one set of data (e.g. GSAS
   .gsa/.fxye files with multiple banks or image files with multiple images.)
 * ``usedRanIdList``: a list of previously used random Id values that can be checked to determine that a value is unique. 

Note that a Reader is used to read only a single phase, image,
dataset, etc. and will be called repeatedly when used to read files
that contain multiple datasets, etc. The ``buffer`` dict can be used
to hold information that will speed repeated calls. 
As an example, the ``buffer`` dict is used in CIF reading to hold the parsed CIF file,
so that when reading multiple datasets or phases from a multi-block
CIF, the parsed information can be reused without having to reread and
reparse the file for subsequent calls. For multi-image files, the
indexing of where images are to be found is done once and saved. This
greatly speeds the time needed to process image file. 

Some additional information specific to on what a ``Reader()`` method
should do for images and single-crystal datasets can be found in the
documentation for :class:`~GSASIIobj.ImportImage` (images) and 
:class:`~GSASIIobj.ImportStructFactor`, respectively. 

Reader return values
______________________

The ``Reader`` routine should return the value of True if the file has been
read successfully. Optionally, use `self.warnings` to indicate any
problems. 

If the file cannot be read,  the ``Reader`` routine should
return False or raise an :meth:`GSASIIobj.ImportBaseclass.ImportException`
exception. (Why are both available? Sometimes an exception is the easiest way to
bail out of a called routine.) Place text in `self.errors` and/or use:: 

     ImportException('Error message')

to give the user information on what went wrong during the
reading. The following variables are used to indicate results from the reader:

self.warnings
^^^^^^^^^^^^^^^^^^^^^

Use `self.warnings` to indicate any information
that should be displayed to the user if the file is read successfully,
but perhaps not completely or additional settings will need to be
made. 

self.errors
^^^^^^^^^^^^^^^^^^^^^

Use `self.errors` to give the user information on where and why a read
error occurs in the file. Note that text supplied with the ``raise``
statement will be appended to ``self.errors``. 

self.repeat
^^^^^^^^^^^^^^^^^^^^^

Set `self.repeat` to True (the default is False) if a Reader should be
called again to after reading to indicate that more data may exist in
the file to be read. This is used for reading multiple powder
histograms or multiple images from a single file. Variable
`self.repeatcount` is used to keep track of the block numbers.

*Reader support routines*
____________________________________

Note that module :mod:`~GSASII.GSASIIctrlGUI` supplies three GUI routines, 
:meth:`~GSASII.GSASIIctrlGUI.BlockSelector`, 
:meth:`~GSASII.GSASIIctrlGUI.MultipleBlockSelector`, and 
:meth:`~GSASII.GSASIIctrlGUI.MultipleChoicesDialog` that are useful for 
selecting amongst one or more datasets (and perhaps phases) or data items for 
``Reader()`` routines that may encounter more than one set of information
in a file. 

.. _ContentsValidator:  

ContentsValidator()
--------------------

Defining a ``ContentsValidator`` method is optional, but is usually a
good idea, particularly if the file extension is not a reliable
identifier for the file type. The intent of this routine is to take a
superficial look at the file to see if it has the expected
characteristics of the expected file type. For example, are there
numbers in the expected places? 

This routine is passed a single argument:

* `filepointer`: a file object [created by Python intrinsic :func:`open`] that accesses
  the file and is points to the beginning of the file when ContentsValidator is
  called. 

Note that :meth:`~GSASIIobj.ImportBaseclass.CIFValidator`
found in :class:`GSASIIobj.ImportBaseclass` is a ContentsValidator
for validating CIF files.

ContentsValidator return values
________________________________

The ``ContentsValidator`` routine should return the value of True if
the file appears to match the type expected for the class. 

If the file cannot be read by this class,  the routine should
return False. Preferably one will also place text in `self.errors` 
to give the user information on what went wrong during the reading.

ReInitialize()
--------------------

Importer classes are substantiated only once and are used as needed.
This means that if something needs to be initialized before the
``Reader()`` will be called to read a new file, the initialization step must be coded. The
``ReInitialize()`` method is provided for this and it is always called
before the ``ContentsValidator`` method is called. Use care to call
the parent class ``ReInitialize()`` method, if this is overridden. 

======================================
 Phase Importer Routines
======================================

Phase importer routines are classes derived from
:class:`GSASIIobj.ImportPhase`.  
They must be found in files named `G2phase*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

*Module G2phase* (PDB, PDF, .EXP & JANA m40,m50)
---------------------------------------------------------

A set of short routines to read in phases in a number of formats, as
created by a number of external programs. The PDB (Protein Data Bank)
format uses coordinates in Angstroms. Other software includes
GSAS/EXPGUI .EXP, JANA .m40-.m50 and Powder Diffaction File entries from the ICDD.

.. automodule:: GSASII.imports.G2phase
    :members: 
    :synopsis: Coordinates/Cell import from common external software
	       

*Module G2phase_GPX: Import phase from GSAS-II project*
--------------------------------------------------------

Copies a phase from another GSAS-II project file into the
current project.	       

.. automodule:: GSASII.imports.G2phase_GPX
    :members: 
    :synopsis: Reads phase information from a GSAS-II project (.gpx) file
      a text file.

      
*Module G2phase_CIF: Coordinates from CIF*
------------------------------------------

Parses a CIF using  PyCifRW from James Hester (https://github.com/jamesrhester/pycifrw) and pulls out the
structural information.

If a CIF generated by ISODISTORT is encountered, extra information is
added to the phase entry and constraints are generated. 
      
.. automodule:: GSASII.imports.G2phase_CIF
    :members: 
    :synopsis: Reads phase information from a CIF

	       
*Module G2phase_INS: Import phase from SHELX INS file*
--------------------------------------------------------

Copies a phase from SHELX ins file into the current project.

.. automodule:: GSASII.imports.G2phase_INS
    :members: 

*Module G2phase_rmc6f: Import phase from RMCProfile*
--------------------------------------------------------

Copies a phase from a file written by RMCProfile into the current GSAS-II project.

.. automodule:: GSASII.imports.G2phase_rmc6f
    :members:

*Module G2phase_xyz: read coordinates from an xyz file*
-------------------------------------------------------

A short routine to read in a phase from an xyz Cartesian coordinate file

.. automodule:: GSASII.imports.G2phase_xyz
    :members: 

*Module G2phase_RRUFF: read from RRUFF database*
-------------------------------------------------------

A short routine to read in a phase from the RRUFF database
in the native format it uses.

.. automodule:: GSASII.imports.G2phase_RRUFF
    :members: 

======================================
 Powder Data Importer Routines
======================================

Powder data importer routines are classes derived from
:class:`GSASIIobj.ImportPowderData`. 
They must be found in files named `G2pwd*.py` that are in the Python path
and the class must override the ``__init__`` method and add a
``Reader`` method. 

The distributed powder data importers are:

*Module G2pwd_GPX: GSAS-II projects*
------------------------------------
Routine to importer powder data from GSAS-II .gpx files

.. automodule:: GSASII.imports.G2pwd_GPX
    :members: 
    :synopsis: Reads powder data from from a GSAS-II project (.gpx) file

*Module G2pwd_fxye: GSAS data files*
------------------------------------
Routine to read in powder data in a variety of formats
that were defined in the original GSAS/EXPGUI software suite.

.. automodule:: GSASII.imports.G2pwd_fxye
    :members: 
    :synopsis: Reads powder data in all of the GSAS formats

*Module G2pwd_xye: Topas & Fit2D data*
-------------------------------------------

Routine to read in powder data from a number of related formats
including ones used in Topas and Fit2D. Typical file extensions are 
.xye, .qye, .chi, and .qchi.

.. automodule:: GSASII.imports.G2pwd_xye
    :members: 
    :synopsis: Reads powder data from a Topas format file

*Module G2pwd_CIF: CIF powder data*
------------------------------------

Routine to read in powder data from a CIF. 
Parses a CIF using  PyCifRW from James Hester
(https://github.com/jamesrhester/pycifrw).

.. automodule:: GSASII.imports.G2pwd_CIF
    :members: 
    :synopsis: Reads powder data from a CIF

*Module G2pwd_BrukerRAW: Bruker .raw*
---------------------------------------------------

Routine to read in powder data from most Bruker versions 2, 3, or 4 .raw
files.
	       
.. automodule:: GSASII.imports.G2pwd_BrukerRAW
    :members: 
    :synopsis: Reads powder data from a Bruker .raw file

*Module G2pwd_BrukerBRML: Bruker .brml*
---------------------------------------------------

Routine to read in powder data from a Bruker .brml file.
	       
.. automodule:: GSASII.imports.G2pwd_BrukerBRML
    :members: 
    :synopsis: Reads powder data from a Brucker .brml file
               

*Module G2pwd_FP: FullProf .dat data*
-------------------------------------

Routine to read in powder data from a FullProf .dat file

.. automodule:: GSASII.imports.G2pwd_FP
    :members:

*Module G2pwd_Panalytical: Panalytical .xrdml data*
---------------------------------------------------

Routines to importer powder data from a Pananalytical (XML) .xrdm file. 

.. automodule:: GSASII.imports.G2pwd_Panalytical
    :members: 

*Module G2pwd_csv: Read Excel .csv data*
------------------------------------------

Routine to read in powder data from Excel type comma separated variable
column-oriented variable. The only allowed extensions for this are
.csv, .xy, or .XY.

.. automodule:: GSASII.imports.G2pwd_csv
    :members: 

*Module G2pwd_rigaku: powder data from a Rigaku .txt file*
----------------------------------------------------------------

.. automodule:: GSASII.imports.G2pwd_rigaku
    :members:

*Module G2pwd_MIDAS: Read integration results from MIDAS*
----------------------------------------------------------------

The MIDAS package can be used on to read and integrate images,
particularly at APS Sector 1. This importer reads the Zarr container
files created by MIDAS. 

.. automodule:: GSASII.imports.G2pwd_MIDAS
    :members:


======================================
 Single Crystal Data Importer Routines
======================================

Single crystal data importer routines are classes derived from
, :class:`GSASIIobj.ImportStructFactor`.
They must be found in files named `G2sfact*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

*Module G2sfact: simple HKL import*
-----------------------------------
Read structure factors from a number of hkl file types. Routines are
provided to read from files containing F or F\ :sup:`2` values from a
number of sources. 

.. automodule:: GSASII.imports.G2sfact
    :members: 
    :synopsis: Reads single crystal data from simple hkl files

*Module G2sfact_CIF: CIF import*
-----------------------------------
Read structure factors from a CIF reflection table CIF using
PyCifRW from James Hester (https://github.com/jamesrhester/pycifrw).
	       
.. automodule:: GSASII.imports.G2sfact_CIF
    :members: 
    :synopsis: Reads single crystal data from CIF files

=================================================
 Small Angle Scattering Data Importer Routines
=================================================

Small angle scattering data importer routines are classes derived from
, :class:`GSASIIobj.ImportSmallAngleData`.
They must be found in files named `G2sad*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are in:

*Module G2sad_xye: read small angle data*
------------------------------------------------

Routines to read in small angle data from an .xye type file, with
two-theta or Q steps. Expected extensions are .xsad, .xdat, .nsad, or .ndat.

.. automodule:: GSASII.imports.G2sad_xye
    :members: 
    :synopsis: Reads small angle scattering data from simple files

======================================
 Image Importer Routines
======================================

Image importer routines are classes derived from
:class:`GSASIIobj.ImportImage`. 
See :ref:`Writing a Importer Routine<import_routines>` for general
information on importers and the :class:`GSASIIobj.ImportImage` for
information on what class variables a reader should set. 
Image importers must be found in files named `G2img*.py` that are in the Python path
and the class must override the ``__init__`` method and add a
``Reader`` method.

The distributed routines are:

*Module G2img_ADSC: .img image file*
--------------------------------------

.. automodule:: GSASII.imports.G2img_ADSC
    :members: 

*Module G2img_EDF: .edf image file*
--------------------------------------

.. automodule:: GSASII.imports.G2img_EDF
    :members: 

*Module G2img_SumG2: Python pickled image*
------------------------------------------

Routine to read an image from GSAS-II that has been pickled in Python. Images
in this format are created by the "Sum image data" command. At least for
now, only one image is permitted per file.
    
.. automodule:: GSASII.imports.G2img_SumG2
    :members: 

*Module G2img_GE: summed GE image file*
---------------------------------------

Read data from General Electric angiography x-ray detectors,
primarily as used at APS 1-ID. 
This shows an example of an importer that will handle files with
more than a single image. 

.. automodule:: GSASII.imports.G2img_GE
    :members: 

*Module G2img_MAR: MAR image files*
--------------------------------------

.. automodule:: GSASII.imports.G2img_MAR
    :members: 

*Module G2img_Rigaku: .stl image file*
--------------------------------------

.. automodule:: GSASII.imports.G2img_Rigaku
    :members: 

*Module G2img_1TIF: Tagged-image File images*
--------------------------------------------------

Routine to read an image in Tagged-image file (TIF) format as well as a variety
of slightly incorrect pseudo-TIF formats used at instruments around the world.
This uses a custom reader that attempts to determine the instrument and detector
parameters from various aspects of the file, not always successfully alas. 

.. automodule:: GSASII.imports.G2img_1TIF
    :members: 

*Module G2img_PILTIF: Std Tagged-image File images*
-----------------------------------------------------

Routine to read an image in Tagged-image file (TIF) format using a standard
image library function in Pillow or the now obsolete PIL package.
This means that parameters such as the pixel size
(which is in the TIFF header but is almost never correct) 
and distance to sample, etc. are not correct unless specified in a 
separate metadata file. See below for more information on metadata
files. 

.. automodule:: GSASII.imports.G2img_PILTIF
    :members: 

*Module G2img_png: png image file*
---------------------------------------

Routine to read an image in .png (Portable Network Graphics) format.
For now, the only known use of this is with converted Mars Rover (CheMin)
tif files, so default parameters are for that.

.. automodule:: GSASII.imports.G2img_CheMin
    :members: 

*Module G2img_CBF: .cbf cif image file*
---------------------------------------

.. automodule:: GSASII.imports.G2img_CBF
    :members: 

*Module G2img_HDF5: summed HDF5 image file*
-------------------------------------------

Reads images found in a HDF5 file. If the file contains multiple
images, all are read. 

.. automodule:: GSASII.imports.G2img_HDF5
    :members: 

*Module G2img_SFRM: Brucker .sfrm image file*
-------------------------------------------------

.. automodule:: GSASII.imports.G2img_SFRM
    :members: 

*Module G2img_pixirad_1ID_16bit: Pixirad detector*
----------------------------------------------------------------

Reads images from the pixirad detector in use as APS Sector 1. 

.. automodule:: GSASII.imports.G2img_pixirad_1ID_16bit
    :members:


======================================================
 Pair Distribution Function (PDF) Importer Routines
======================================================

PDF importer routines are classes derived from
:class:`GSASIIobj.ImportPDFData`. 
See :ref:`Writing a Importer Routine<Import_Routines>` for general information on importers. 

The distributed routines are in:

*Module G2pdf_gr: read PDF G(R) data*
------------------------------------------------

Routines to read in G(R) data from a pdfGet/GSAS-II .gr or gudrun .dat
file (with :math:`\AA` steps) or S(Q) data from a .fq file.

.. automodule:: GSASII.imports.G2pdf_gr
    :members: 

======================================
 Reflectometry Importer Routines
======================================

Reflectometry importer routines are classes derived from
:class:`GSASIIobj.ImportReflectometryData`. 
See :ref:`Writing a Importer Routine<Import_Routines>` for general information on importers. 

The distributed routines are:

*Module G2rfd_xye: read reflectometry data*
------------------------------------------------

Routines to read in reflectometry data from an
.xrfd, .xdat, .xtrfd, .xtdat, .nrfd or .ndat type file, with
two-theta or Q steps. 

.. automodule:: GSASII.imports.G2rfd_xye
    :members: 

*Module G2rfd_Panalytical: read Panalytical reflectometry data*
-------------------------------------------------------------------

Routine to importer reflectivity data from a Panalytical .xrdm (xml)
file.

.. automodule:: GSASII.imports.G2rfd_Panalytical
    :members: 

*Module G2rdf_rigaku: reflectometry data from a Rigaku file*
----------------------------------------------------------------

.. automodule:: GSASII.imports.G2rfd_rigaku
    :members:
