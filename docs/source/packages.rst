Required Packages & Supported Platforms
==========================================

GSAS-II requires a standard Python interpreter to be installed, as
well as several separately-developed packages. GSAS-II is being
developed using Python 3.7 and 3.9. At this point we think that 
most sections of the code have been exercised in Python 3,
but  bugs are still expected (please report them). We are no longer
testing with Python 2.7 and strongly urge everyone to upgrade,
but if problems running GSAS-II in Python 2.7 are reported, we will 
consider making that the code compliant with both. 

Note that the packages listed below are not distributed as part of the Python standard
library. We have been depending on the free Anaconda
Python (https://www.anaconda.com/)
distribution (and provide installers based on that), but Anaconda does
not seem to be supplying up to date versions of the wxpython package
that the GUI requires and does not have versions for all supported
platforms. Use of the miniforge
(https://github.com/conda-forge/miniforge) distribution is recommended
where needed. 

There are many other Python distributions, such as Enthought Inc.'s Canopy and
Python(x,y), see here:
https://www.python.org/download/alternatives/. We are no longer using
any of them and are unsure of how they will function. Some very old
GSAS-II installations were based on the quite outdated Enthought Python Distribution
(EPD); this is known to have some problems with reading CIFs and we
strongly encourage updating from that.

More details on allowed and prefered package versions can be found in
the documentation for variable :attr:`GSASIIdataGUI.versionDict`.

GUI Requirements
----------------

When using the GSAS-II graphical user interface (GUI), the following
Python extension packages are required:

* wxPython (http://wxpython.org/docs/api/). Note that GSAS-II has been tested with wxPython 2.8.x, 3.0.x and 4.0.x. We encourage use of 3.0 with Python 2.7 and 4.x with Python 3.x. 
* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/),
* matplotlib (http://matplotlib.org/contents.html)  and
* PyOpenGL (http://pyopengl.sourceforge.net/documentation). Note: a copy of this is distributed with GSAS-II (at present) and will be installed if the Python setuptools package is present. 

Several packages are used in sections of the code, but are not
required. If these packages are not present, warning messages may be
generated if they would be needed, but the vast bulk of GSAS-II will function normally. 

* Pillow (https://pillow.readthedocs.org) or PIL (http://www.pythonware.com/products/pil/). This is used to save
  and read certain types of images.
* h5py is the HDF5 interface and hdf5 is the support package. These
  packages are (not surprisingly) required
  to import images from HDF5 files. If these libraries are not present,
  the HDF5 importer(s) will not appear in the import menu and a
  warning message appears on GSAS-II startup. 
* imageio is used to make movies. 
* svn: When using Anaconda we also encourage installation of the
  svn (subversion) conda package. This is not actually part of Python
  and can be installed directly into your system's configuration. It is used by
  GSAS-II to download updates to our code. This can be skipped if svn
  is installed directly (easy Linux, but a bit harder on MacOS and
  Windows). In conda-forge this package is called subversion, but at
  present is only available for Linux.
* pywin32 (windows only): this provides the win32com module that is
  used to install GSAS-II on windows machines. GSAS-II can be used on
  Windows without this, but the installation will offer less
  integration into Windows. 
* conda: the conda package allows access to conda features from
  inside Python. It will be used inceasingly by GSAS-II to
  self-install software. The conda package is installed by default in
  miniconda and anaconda but if you create an environment for GSAS-II
  (`conda create -n <env> package-list...`), it will not be added
  unless you request it specifically.  

*Conda command*:
  Here is a typical conda command used to install a GSAS-II compatible
  Python interpreter::

    conda install python=3.9 wxpython numpy scipy matplotlib pyopengl pillow h5py imageio subversion -c conda-forge
    
  or to put a Python configured for GSAS-II into a separate conda
  environment (here named ``g2python``, but any name can be used), use
  command::

    conda create -n g2python python=3.9 wxpython numpy scipy matplotlib pyopengl  pillow h5py imageio conda subversion -c conda-forge 

Note that at present the subversion is only available for Linux, so
that should be removed from the commands above. For windows add pywin32
Also, while there is no
reason why GSAS-II should not run with Python 3.10, we are not yet
providing binaries for this. 
   
Remember to activate using: ``<path>\Scripts\activate``  (windows); 
``source <path>/bin/activate`` (Mac/Linux). Note that one should add
``g2python`` (etc.) at the end if using a conda environment.

Note that svn seems to be unsupported these days by Anaconda. For
Linux and MacOS, use subversion in conda-forge rather than svn. No
good solution yet for Windows.

Scripting  Requirements
-----------------------

When using the GSAS-II scripting interface (:mod:`GSASIIscriptable`),
only the following Python extension packages are required:

* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/).

Note that a few sections of the code require matplotlib (http://matplotlib.org/contents.html), Pillow
(https://pillow.readthedocs.org) (or PIL,
http://www.pythonware.com/products/pil/), or h5py + hdf5, but none of
these are required to run scripts and the vast
majority of scripts will not need these packages.

Optional Packages
-----------------------

* Sphinx (https://www.sphinx-doc.org) is used to generate the
  documentation you are currently reading. Generation of documentation
  is not generally something needed by users or even most code developers.

 * SCons (https://scons.org/) is used to compile the small amount of
   Fortran code that is included with GSAS-II. Use of this is
   discussed in the last section of this chapter.


External and Supplied Programs
--------------------------------

GSAS-II provides interfaces to use a number of programs developed by
others. Some are included with GSAS-II and others must be installed
separately. When these programs are accessed, citation
information is provided. 

GSAS-II includes copies of these programs:

  **DIFFaX**
    Simulate layered structures with faulting
    
  **CifFile**
    A software library used to read data and structures from CIF
    
  **Shapes**
    Model small angle scattering with shaped particles
    
  **NIST FPA**
    Use Fundamental Parameters to determine GSAS-II profile function 

No additional steps beyond a standard installation
are needed to access their functionality.

**Bilboa Crystallographic Server**: GSAS-II directly access the
Bilboa Crystallographic Server (provided
the computer has internet access). This allows automated use of the
k-SUBGROUPSMAG, k-SUBGROUPS and PseudoLattice web utilities for
computation of space group subgroups, color (magnetic) subgroups &
lattice search.

At the request of the program authors, these programs are not included
with GSAS-II and must be installed separately:

  **RMCProfile**
    Large-box PDF & S(Q) fitting. We have heard from users that V6.7.7
    of RMCProfile is compatible with the input created by GSAS-II,
    but not V6.7.9.

  **fullrmc**
    A modern software toolkit for large-box PDF & S(Q) fitting. Use
    version 5.0 or later. The implementation for this is not
    completed. 

  **Dysnomia**
    Computes enhanced Fourier maps with Maximum Entropy estimated
    extension of reflection sphere

  **PDFfit2**
  Small-box fitting of PDFs. This code is no longer supported, but is
  still quite useful. It can be installed from conda into Python
  versions up to Python 3.7, but is supplied for Windows within
  GSAS-II for Python 3.7, 3.8 and 3.9 and for MacOS only with Python
  3.7.

  For other platforms/Python versions, it is probably best to use a
  separate Python interpreter. 
    
Supported Platforms
--------------------------------

It should be possible to run GSAS-II on any computer where Python 3.7+ and
the appropriate required packages are available. For many platforms,
binary versions of the Fortran code used in GSAS-II are supplied, but the
binaries must match the platform and the major versions of both Python and
numpy; even for supported platforms; not all combinations are
provided. Should one wish to run GSAS-II where binary files are not
supplied, compilation will be needed. This will require the GNU Fortran (gfortran)
compiler (https://gcc.gnu.org/fortran/) as well as the Python SCons
package. Instructions are supplied for a number of platforms (such as 
https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallLinux#CompilingFortranCode). Note
that there are prepackaged versions of GSAS-II for most common
platforms. These include Python, all required and most optional
packages and a version of all files needed to run GSAS-II -- albeit
not usually the current version. The
installation process will try to update to the current version, if the
computer where installation is occuring has internet access. 

At present the following platforms are directly supported:

* **Windows-10**: Installation kits are available for both 32-bit and
  64-bit windows. Running GSAS-II on older versions of Windows is
  likely possible, but to do so one must locate compatible versions of Python
  and packages. This is getting increasingly tough. We have not tried
  Windows-11, but expect the Windows-10 versions to run there.

* **MacOS**: We provide an installer for Macs with Intel
  processors. This can also be used on ARM-equipped Macs ("M1" or "Apple
  Silicon" processors) but native M1 code runs way
  faster. Installation on the native ARM code is more complex; our
  instructions (https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/MacM1Notes)
  require that the homebrew package installer be installed and then
  brew (homebrew) be used to install Python and a number of packages.
  Macs older than
  Catalina (10.15) may require older distributions of Python. 

* **Linux** (Intel-compatible): GSAS-II does not get a lot of testing in Linux by us, but is
  fairly widely used on this platform nonetheless.  One can use the
  installer that we provide, but compatibility with older and very new
  versions of OSes can be tough and may require compatibility
  libraries. At times it may be better to use the Linux distribution's
  versions of Python and packages. This is typically done with a
  software tool such as apt or yum. An example on how to do this is
  shown for the Raspberry Pi.

* **Raspberry Pi** (ARM) Linux: GSAS-II has been installed on both 32-bit
  and the experimental 64-bit version of the Raspberry Pi OS (formerly
  called Raspbian) and compiled binaries are provided. It should also
  run with Ubuntu Linux for this platform, but this has not been
  tried. It is necessary to use the Raspbian Linux distribution's
  versions of Python and its packages. Instructions are provided
  (https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallPiLinux). 
  The performance of GSAS-II on a Raspberry Pi is not blindingly fast,
  but one can indeed run GSAS-II on a computer that costs only $15!
 
