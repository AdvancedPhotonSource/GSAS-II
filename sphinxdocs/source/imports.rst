*GSAS-II Import Modules*
====================================

Imports are implemented by deriving a class from 
:class:`GSASIIIO.ImportPhase`, :class:`GSASIIIO.ImportStructFactor`
or :class:`GSASIIIO.ImportPowderData` (which are in turn 
derived from :class:`GSASIIIO.ImportBaseclass`)
to implement import of 
a phase, a single crystal or a powder dataset, respectively. 
Module file names (`G2phase_`, `G2pwd_` and `G2sfact_`) are used to
determine which menu an import routine should be placed into. (N.B. this
was an unnecessary choice; this could be done from the class used.)

Writing an Import Routine
--------------------------

.. _Import_routines:


When writing a import routine, one should create a new class derived
from :class:`GSASIIIO.ImportPhase`, :class:`GSASIIIO.ImportStructFactor`
or :class:`GSASIIIO.ImportPowderData`. As described below, 
all these classes will implement
an ``__init__()`` and a ``Reader()`` method, and most will supply a 
``ContentsValidator()`` method, too.
See the :class:`~GSASIIIO.ImportPhase`, :class:`~GSASIIIO.ImportStructFactor`
or :class:`~GSASIIIO.ImportPowderData` class documentation 
for details on what values each type of ``Reader()`` should set. 

__init__()
~~~~~~~~~~~~~~
 
The class should supply a
``__init__`` method which calls the parent ``__init__`` method and
specifies the following parameters: 

  * `extensionlist`: a list of extensions that may be used for this type of file.
  * `strictExtension`: Should be True if only files with extensions in
    `extensionlist` are allows; False if all file types should be offered
    in the file browser. Also if False, the import class will be
    used on all files when "guess from format" is tried, though 
    readers with matching extensions will be tried first.
  * `formatName`: a string to be used in the menu. Should be short. 
  * `longFormatName`: a longer string to be used to describe the format in help. 

Reader()
~~~~~~~~~~~~~~

The class must supply a ``Reader`` method that actually performs the
reading. All readers must have at a minimum these arguments::

    def Reader(self, filename, filepointer, ParentFrame, **unused):

where the arguments have the following uses: 

 * `filename`: a string with the name of the file being read
 * `filepointer`: a file object (created by :func:`open`) that accesses
   the file and is points to the beginning of the file when Reader is
   called. 
 * `ParentFrame`: a reference to the main GSAS-II (tree) windows, for
   the unusual ``Reader`` routines that will create GUI windows to ask
   questions. 

In addition, the following keyword parameters are defined that ``Reader``
routines may optionally use:

 * `buffer`: a dict that can be used to retain information between repeated calls of the routine
 * `blocknum`: counts the number of times that a reader is called 
 * `usedRanIdList`: a list of previously used random Id values that can be checked to determine that a value is unique. 

As an example, the `buffer` dict is used for CIF reading to hold the parsed CIF file
so that it can be reused without having to reread the file from
scratch. 

Reader return values
______________________

The ``Reader`` routine should return the value of True if the file has been
read successfully. Optionally, use `self.warnings` to indicate any
problems. 

If the file cannot be read,  the ``Reader`` routine should
return False or raise an :class:`GSASIIIO.ImportBaseclass.ImportException`
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
called again to read a second block from a file. Most commonly 
(only?) used for reading multiple powder histograms from a single
file. Variable `self.repeatcount` is used to keep track of the block
numbers. 

*support routines*
_________________________
Note that the base class (:class:`GSASIIIO.ImportBaseclass`) supplies two routines, 
:meth:`~GSASIIIO.ImportBaseclass.BlockSelector` and 
:meth:`~GSASIIIO.ImportBaseclass.MultipleBlockSelector` that are useful for 
selecting amongst one or more datasets (and perhaps phases) for 
``Reader()`` routines that may encounter more than one set of information
in a file. 
Likewise, when an operation will take some time to complete, 
use :meth:`~GSASIIIO.ImportBaseclass.ShowBusy` and 
:meth:`~GSASIIIO.ImportBaseclass.DoneBusy` to show the user
that something is happening. 


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

Note that :meth:`GSASIIIO.ImportBaseclass.CIFValidator` is a ContentsValidator
for validating CIF files. 


ContentsValidator return values
________________________________

The ``ContentsValidator`` routine should return the value of True if
the file appears to match the type expected for the class. 

If the file cannot be read by this class,  the routine should
return False. Preferably one will also place text in `self.errors` 
to give the user information on what went wrong during the reading.

Currently Defined Phase Import Routines
----------------------------------------
Phase import routines are classes derived from
:class:`GSASIIIO.ImportPhase`.  
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

Currently Defined Powder Data Import Routines
---------------------------------------------
Powder data import routines are classes derived from
:class:`GSASIIIO.ImportPowderData`. 
They must be found in files named `G2pwd*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
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

Currently Defined Single Crystal Data Import Routines
-----------------------------------------------------
Single crystal data import routines are classes derived from
, :class:`GSASIIIO.ImportStructFactor`.
They must be found in files named `G2sfact*.py` that are in the Python path
and the class must override the ``__init__`` method and add a ``Reader`` method.
The distributed routines are:

.. automodule:: G2sfact
    :members: 
    :synopsis: Reads single crystal data from simple hkl files

.. automodule:: G2sfact_CIF
    :members: 
    :synopsis: Reads single crystal data from CIF files
