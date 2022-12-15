Required Packages & Supported Platforms
==========================================

GSAS-II requires a standard Python interpreter to be installed, as
well as several separately-developed packages, as described
below. While for some packages, we have not seen much dependence on
versions, for others we do find significant differences, for those we
are currently recommending the following versions:

 * Python 3.10 is recommended, but 3.7 or later is fine
 * wxPython 4.1.1, but 3.0 or later should be OK
 * NumPy 1.23 recommended, but anything from 1.17 on is likely fine
 * matplotlib 3.6 is recommended, 3.4 or later preferred. 

Note that GSAS-II is being developed using Python 3.7, 3.9 and 3.10. We are no longer
supporting Python 2.7 and <=3.6, and recomment upgrading. More details on problematic package versions can be found in
the documentation for variable :attr:`GSASIIdataGUI.versionDict`.

There are a number of ways to install Python along with the packages
needed by GSAS-II. We had been creating installers utilizing the Anaconda
Python (https://www.anaconda.com/)
distribution, but no longer provides modern versions that include
wsPython or subversion, so we are transitioning to the
community-supported conda-forge library of packages and the miniforge
distribution. On MacOS, homebrew can be used for Python and most
needed packages, while on Linux, the native package installers
(apt-get or yum, etc.) offer the same. Any packages not provided in
that fashion can be installed with Python's pip mechanism. 
Other alternatives Python packages include Enthought Inc.'s Canopy and
Python(x,y), see here:
https://www.python.org/download/alternatives/. We are no longer using
any of them and are unsure of how they will function. 


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
* requests: this package simplifies http access
  (https://requests.readthedocs.io/). It is used for access to
  webpages such as ISODISTORT and for some internal software downloads.
  
*Conda command*:
  Here is a typical conda command used to install a GSAS-II compatible
  Python interpreter::

    conda install python=3.9 wxpython numpy scipy matplotlib pyopengl pillow h5py imageio subversion requests -c conda-forge
    
  or to put a Python configured for GSAS-II into a separate conda
  environment (here named ``g2python``, but any name can be used), use
  command::

    conda create -n g2python python=3.9 wxpython numpy scipy matplotlib pyopengl  pillow h5py imageio conda subversion requests -c conda-forge 

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
    Simulate layered structures with faulting. https://www.public.asu.edu/~mtreacy/DIFFaX.html
    
  **PyCifRW**
    A software library used to read data and structures from
    CIF. https://bitbucket.org/jamesrhester/pycifrw
    
    
  **Shapes**
    Model small angle scattering with shaped particles. 
    
  **NIST FPA**
    Use Fundamental Parameters to determine GSAS-II profile function 

  **NIST*LATTICE**
   Searches for higher symmetry unit cells and possible relationships
   between unit cells. An API has been written and this will be
   integrated into the GSAS-II GUI. 

No additional steps beyond a standard installation
are needed to access their functionality.

  **Bilboa Crystallographic Server** (https://www.cryst.ehu.es):
  GSAS-II directly access the 
  Bilboa Crystallographic Server (provided
  the computer has internet access). This allows automated use of the
  k-SUBGROUPSMAG, k-SUBGROUPS and PseudoLattice web utilities for
  computation of space group subgroups, color (magnetic) subgroups &
  lattice search.

  **BYU ISOTROPY Software Suite**
  (https://stokes.byu.edu/iso/isotropy.php): GSAS-II directly access
  capabilities in the ISOTROPY Software Suite from Brigham Young
  University for representational analysis and magnetism analysis. 

At the request of the program authors, other programs that can be
access within GSAS-II are not included
as part of the GSAS-II distribution and must be installed separately:

  **RMCProfile**
    Large-box PDF & S(Q) fitting. The GSAS-II interface was originally
    written for use with release 6.7.7 of RMCProfile, but updates have
    been made for compatible with 6.7.9 as well.

    RMCProfile must be downloaded by the user from
    http://rmcprofile.org/Downloads or
    https://rmcprofile.pages.ornl.gov/nav_pages/download/

  **fullrmc**
    A modern software framework for large-box PDF & S(Q) fitting. Note
    that the GSAS-II implementation is not compatible with the last
    open-source version of fullrmc, but rather the version 5.0 must be
    used, which is distributed as a compiled versions for 64-bit
    Intel-compatible processors running Windows, Linux and MacOS from
    website
    https://github.com/bachiraoun/fullrmc/tree/master/standalones. GSAS-II
    will offer to install this software into the binary directory when the fullrmc
    option is selected on the Phase/RMC tab. 

  **Dysnomia**
    Computes enhanced Fourier maps with Maximum Entropy estimated
    extension of reflection sphere. See https://jp-minerals.org/dysnomia/en/.

  **PDFfit2**
  Small-box fitting of PDFs; see
  https://github.com/diffpy/diffpy.pdffit2#pdffit2. This code is no
  longer supported, but is 
  still quite useful. It can be installed from conda into Python
  versions up to Python 3.7, but is supplied for Windows within
  GSAS-II for Python 3.7, 3.8 and 3.9 and for MacOS only with Python
  3.7.

  For other platforms/Python versions, it is probably best to use a
  separate Python interpreter. If GSAS-II is installed with the conda
  package manager (the usual installation practice), the GUI will
  offer the option to install PDFfit2 with Python 3.7 in a separate
  environment when the option is selected on
  the Phase/RMC tab. 

Supported Platforms
--------------------------------

It should be possible to run GSAS-II on any computer where Python 3.7+ and
the appropriate required packages are available. GSAS-II requires that
some code must be compiled. For the following platforms, binary images
are provided:

  * Windows-10: 64-bit Intel-compatible processors 
  * MacOS:  Intel processors 
  * MacOS: Apple Silicon (M1, etc) processors 
  * Linux: 64-bit Intel-compatible processors
  * Linux: ARM processors (64-bit and 32-bit Raspberry Pi)

Note that these binaries must the major versions of both Python and
numpy; only a small number of combinations are provided.
Should one wish to run GSAS-II where binary files are not
supplied (such as 32-bin Windows or Linux) or with other versions of
Python/NumPy, compilation will be needed but the user.
This will require the GNU Fortran (gfortran)
compiler (https://gcc.gnu.org/fortran/) as well as the Python SCons
package. General instructions are supplied for a number of platforms (such as 
https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallLinux#CompilingFortranCode).

More details on platforms is discussed below:

* **Windows-10**: Installation kits are provided for 
  64-bit windows. An installation kit with older Python versions
  is provided for 32-bit Window; this cannot be updated but GSAS-II
  will be updated if installed on a computer with internet access.
  Running GSAS-II on older versions of Windows is
  likely possible, but to do so one must locate compatible versions of Python
  and packages. This is getting increasingly tough. We have not tried
  Windows-11, but expect the Windows-10 versions to run there.

* **MacOS**: GSAS-II can run natively on Intel or ARM ("M1" or "Apple
  Silicon") processors. With the native code, Mac ARM machines offer
  the highest performance on any platform. 
  
  For Intel processor Macs, we provide an installer. This can also be
  used on ARM-equipped Macs but native M1 code runs way
  faster. Native ARM code installation is more complex; 
  instructions are provided
  (https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/MacM1Notes)
  that require use of either the miniforge package or the homebrew
  package installer. 
  Macs older than Catalina (10.15) will likely require older
  distributions of Python.  

* **Linux** (Intel-compatible): GSAS-II does not get a lot of testing
  in Linux by us, but is fairly widely used on this platform
  nonetheless.  One can use the 
  installer that we provide, but compatibility with older and very new
  versions of Linux can be tough and may require compatibility
  libraries. At times it may be better to use the Linux distribution's
  versions of Python and packages. This is typically done with a
  software tool such as apt or yum. An example on how to do this is
  shown for the Raspberry Pi.

* **Raspberry Pi** (ARM) Linux: GSAS-II has been installed on both 32-bit
  and the 64-bit version of the Raspberry Pi OS (formerly
  called Raspbian) and compiled binaries are provided. Note that
  64-bit is preferred on the models where it can be run (currently
  including  models 3A+, 3B, 3B+, 4, 400, CM3, CM3+, CM4, and Zero
  2 W)  It should also 
  run with Ubuntu Linux for this platform, but this has not been
  tried. For 32-bit  Raspberry Pi OS, it is necessary to use the distribution's
  versions of Python and its packages. Instructions are provided
  (https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallPiLinux). 
  The performance of GSAS-II on a Raspberry Pi is not blindingly fast,
  but one can indeed run GSAS-II on a motherboard that costs only $15
  and uses <5 Watts!

