*GSAS-II Import Modules*
====================================

Imports are implemented by deriving a class from 
:class:`GSASIIobj.ImportPhase`, :class:`GSASIIobj.ImportStructFactor`,
:class:`GSASIIobj.ImportPowderData` ,
:class:`GSASIIobj.ImportSmallAngleData`,  
:class:`GSASIIobj.ImportReflectometryData`,  
:class:`GSASIIobj.ImportPDFData`,  
or :class:`GSASIIobj.ImportImage` (which are in turn 
derived from :class:`GSASIIobj.ImportBaseclass`)
to implement import of 
a phase, a single crystal or a powder dataset, respectively. 
Module file names (`G2phase_`, `G2pwd_` and `G2sfact_`, etc.) are used to
determine which menu an import routine should be placed into. (N.B. this
naming was an unnecessary choice; importer types could be determined
from the base class.)

Most importers are listed below by type (provided this documentation is
up to date), but note that since modules
may be loaded from anywhere in the path, your installation could have
locally-defined importers as well.

.. _import_routines: 

Writing an Import Routine
--------------------------

When writing a import routine, one should create a new class derived
from
:class:`GSASIIobj.ImportPhase`, :class:`GSASIIobj.ImportStructFactor`,
:class:`GSASIIobj.ImportPowderData` ,
:class:`GSASIIobj.ImportSmallAngleData`,  
:class:`GSASIIobj.ImportReflectometryData`,  
:class:`GSASIIobj.ImportPDFData`,  
or :class:`GSASIIobj.ImportImage`. As described below, 
all these classes will implement
an ``__init__()`` and a ``Reader()`` method, and most will supply a 
``ContentsValidator()`` method, too.
See the appropriate class documentation 
for details on what values each type of ``Reader()`` should set. 

__init__()
~~~~~~~~~~~~~~
 
The ``__init__`` method will follow standard boilerplate: 

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
    ``extensionlist`` are allows; False if all file types should be offered
    in the file browser. Also if False, the import class will be
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
importer from being used or included in the expected menu. 

Reader()
~~~~~~~~~~~~~~

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
   questions. The Reader should do something reasonable such as take a
   reasonable default if ``ParentFrame`` is None, which indicates that
   GUI should not be accessed. 

In addition, the following keyword parameters are defined that ``Reader``
routines may optionally use:

 * ``buffer``: a dict that can be used to retain information between repeated calls of the routine
 * ``blocknum``: counts the number of times that a reader is called, to
   be used with files that contain more than one set of data (e.g. GSAS
   .gsa/.fxye files with multiple banks or image files with multiple images.)
 * ``usedRanIdList``: a list of previously used random Id values that can be checked to determine that a value is unique. 

As an example, the ``buffer`` dict is used in CIF reading to hold the parsed CIF file
so that it can be reused without having to reread the file from
scratch. 

Reader return values
______________________

The ``Reader`` routine should return the value of True if the file has been
read successfully. Optionally, use `self.warnings` to indicate any
problems. 

If the file cannot be read,  the ``Reader`` routine should
return False or raise an :class:`GSASIIobj.ImportBaseclass.ImportException`
exception. (Why either? Sometimes an exception is the easiest way to
bail out of a called routine.) Place text in `self.errors` and/or use:: 

     ImportException('Error message')

to give the user information on what went wrong during the reading.

self.warnings
_____________________

Use `self.warnings` to indicate any information
that should be displayed to the user if the file is read successfully,
but perhaps not completely or additional settings will need to be
made. 

self.errors
_____________________

Use `self.errors` to give the user information on where and why a read
error occurs in the file. Note that text supplied with the ``raise``
statement will be appended to ``self.errors``. 

self.repeat
_____________________

Set `self.repeat` to True (the default is False) if a Reader should be
called again to after reading to indicate that more data may exist in
the file to be read. This is used for reading multiple powder
histograms or multiple images from a single file. Variable
`self.repeatcount` is used to keep track of the block numbers.

*support routines*
_________________________

Note that GSASIIIO supplies three routines, 
:meth:`~GSASIIIO.BlockSelector` 
:meth:`~GSASIIIO.MultipleBlockSelector` and 
:meth:`~GSASIIIO.MultipleChoiceSelector` that are useful for 
selecting amongst one or more datasets (and perhaps phases) or data items for 
``Reader()`` routines that may encounter more than one set of information
in a file. 

.. _ContentsValidator:  

ContentsValidator()
~~~~~~~~~~~~~~~~~~~~

Defining a ``ContentsValidator`` method is optional, but is usually a
good idea, particularly if the file extension is not a reliable
identifier for the file type. The intent of this routine is to take a
superficial look at the file to see if it has the expected
characteristics of the expected file type. For example, are there
numbers in the expected places? 

This routine is passed a single argument:

* `filepointer`: a file object (created by :func:`open`) that accesses
  the file and is points to the beginning of the file when ContentsValidator is
  called. 

Note that :meth:`GSASIIobj.ImportBaseclass.CIFValidator` is a ContentsValidator
for validating CIF files. 


ReInitialize()
~~~~~~~~~~~~~~~~~~~~

Import classes are substantiated only once and are used as needed.
This means that if something needs to be initialized before the
``Reader()`` will be called to read a new file, it must be coded. The
``ReInitialize()`` method is provided for this and it is always called
before the ``ContentsValidator`` method is called. Use care to call
the parent class ``ReInitialize()`` method, if this is overridden. 


ContentsValidator return values
________________________________

The ``ContentsValidator`` routine should return the value of True if
the file appears to match the type expected for the class. 

If the file cannot be read by this class,  the routine should
return False. Preferably one will also place text in `self.errors` 
to give the user information on what went wrong during the reading.

Phase Import Routines
----------------------------------------
Phase import routines are classes derived from
:class:`GSASIIobj.ImportPhase`.  
They must be found in files named `G2phase*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

.. automodule:: G2phase
    :members: 
    :synopsis: Uses previously implemented code: PDB and GSAS .EXP

.. automodule:: G2phase_GPX
    :members: 
    :synopsis: Reads phase information from a GSAS-II project (.gpx) file
      a text file. 

.. automodule:: G2phase_CIF
    :members: 
    :synopsis: Reads phase information from a CIF

.. automodule:: G2phase_INS
    :members: 



Powder Data Import Routines
---------------------------------------------
Powder data import routines are classes derived from
:class:`GSASIIobj.ImportPowderData`. 
They must be found in files named `G2pwd*.py` that are in the Python path
and the class must override the ``__init__`` method and add a
``Reader`` method. 

The distributed routines are:

.. automodule:: G2pwd_GPX
    :members: 
    :synopsis: Reads powder data from from a GSAS-II project (.gpx) file

.. automodule:: G2pwd_fxye
    :members: 
    :synopsis: Reads powder data in all of the GSAS formats

.. automodule:: G2pwd_xye
    :members: 
    :synopsis: Reads powder data from a Topas format file

.. automodule:: G2pwd_CIF
    :members: 
    :synopsis: Reads powder data from a CIF

.. automodule:: G2pwd_BrukerRAW
    :members: 
    :synopsis: Reads powder data from a Brucker .raw file

.. automodule:: G2pwd_FP
    :members: 

.. automodule:: G2pwd_Panalytical
    :members: 

.. automodule:: G2pwd_csv
    :members: 

.. automodule:: G2pwd_rigaku
    :members: 



Single Crystal Data Import Routines
-----------------------------------------------------
Single crystal data import routines are classes derived from
, :class:`GSASIIobj.ImportStructFactor`.
They must be found in files named `G2sfact*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

.. automodule:: G2sfact
    :members: 
    :synopsis: Reads single crystal data from simple hkl files

.. automodule:: G2sfact_CIF
    :members: 
    :synopsis: Reads single crystal data from CIF files


Small Angle Scattering Data Import Routines
-----------------------------------------------------
Small angle scattering data import routines are classes derived from
, :class:`GSASIIobj.ImportSmallAngle`.
They must be found in files named `G2sad*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

.. automodule:: G2sad_xye
    :members: 
    :synopsis: Reads small angle scattering data from simple files

Image Import Routines
-----------------------------------------------------
Image import routines are classes derived from
:class:`GSASIIobj.ImportImage`. 
See :ref:`Writing a Import Routine<import_routines>` for general
information on importers and the :class:`GSASIIobj.ImportImage` for
information on what class variables a reader should set. 
Image importers must be found in files named `G2img*.py` that are in the Python path
and the class must override the ``__init__`` method and add a
``Reader`` method.

The distributed routines are:

.. automodule:: G2img_ADSC
    :members: 

.. automodule:: G2img_EDF
    :members: 

.. automodule:: G2img_SumG2
    :members: 

.. automodule:: G2img_GE
    :members: 

.. automodule:: G2img_MAR
    :members: 

.. automodule:: G2img_Rigaku
    :members: 

.. automodule:: G2img_1TIF
    :members: 

.. automodule:: G2img_CheMin
    :members: 

.. automodule:: G2img_CBF
    :members: 

.. automodule:: G2img_HDF5
    :members: 

.. automodule:: G2img_SFRM
    :members: 
       
PDF Import Routines
-----------------------------------------------------
PDF import routines are classes derived from
:class:`GSASIIobj.ImportPDFData`. 
See :ref:`Writing a Import Routine<Import_Routines>` for general information on importers. 

The distributed routines are:

.. automodule:: G2pdf_gr
    :members: 

Reflectometry Import Routines
-----------------------------------------------------
Reflectometry import routines are classes derived from
:class:`GSASIIobj.ImportReflectometryData`. 
See :ref:`Writing a Import Routine<Import_Routines>` for general information on importers. 

The distributed routines are:

.. automodule:: G2rfd_xye
    :members: 


