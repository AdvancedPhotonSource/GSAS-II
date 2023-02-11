GSAS-II Requirements and Options
==========================================

Supported Platforms
--------------------------------

It should be possible to run GSAS-II on any computer where Python 3.7+ and
the appropriate required packages are available, as discussed below,
but GSAS-II also requires that some code must be compiled.
For the following platforms, binary images for this compiled code are provided:

  * Windows-10: 64-bit Intel-compatible processors 
  * MacOS: Intel processors 
  * MacOS: ARM processors, aka Apple Silicon (M1, etc) 
  * Linux: 64-bit Intel-compatible processors
  * Linux: ARM processors (64-bit and 32-bit Raspberry Pi OS only)

Details for GSAS-II use on these specific platforms follows below:

* **Windows**: Installation kits are provided for 
  64-bit Windows-10. An installation kit with older Python versions
  is provided for 32-bit Windows-10; this installer cannot be updated
  to provide newer Python versions than the supplied versions but GSAS-II
  will be updated if installed on a computer with internet
  access.  Running GSAS-II on older versions of Windows is
  likely possible, but to do so one must locate compatible versions of Python
  and packages. This is getting increasingly tough. We have not tried
  Windows-11, but expect the Windows-10 distribution to run fine there.

* **MacOS**: GSAS-II can run natively on Intel or ARM ("M1" or "Apple
  Silicon") processors. With the native code, Mac ARM machines offer
  the highest performance seen on any platform. 
  
  For Intel processor Macs, we provide an installer. This can also be
  used on ARM-equipped Macs but native M1 code runs way
  faster. Native ARM code installation is a bit more complex; but 
  detailed instructions are provided
  (https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/MacM1Notes)
  that require use of either the miniforge package or the homebrew
  package installer. 
  Macs older than Catalina (10.15) will likely require older
  distributions of Python.  

* **Linux** (Intel-compatible): Note that GSAS-II does not get a lot of testing
  in Linux by us, but is used fairly widely on this platform
  nonetheless.  One can use the 
  installer that we provide, but compatibility with older and very new
  versions of Linux may require compatibility
  library installation, not always easy to do. It may be
  better to use your Linux distribution's versions of Python and
  packages (typically done with a software tool such as apt or yum.)
  You may possibly need to use pip as well. Adapt from the detailed
  example for how to do for the 32-bit Raspberry Pi OS:
  https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallPiLinux.  

* **Raspberry Pi** (ARM) Linux: GSAS-II has been installed on both 32-bit
  and the 64-bit version of the Raspberry Pi OS (formerly
  called Raspbian) and compiled binaries are provided.
  It may be possible to use these binaries with Ubuntu Linux for
  this platform, but this has not been tried.
  The performance of GSAS-II on a Raspberry Pi is not blindingly fast,
  but one can indeed run GSAS-II on a motherboard that costs only $15
  (perhaps even one that costs $5) and uses <5 Watts! 

  Note that the 64-bit OS is preferred on the models where it can be run
  (currently including models 3A+, 3B, 3B+, 4, 400, CM3, CM3+, CM4,
  and Zero 2 W) .  With the 32-bit Raspberry Pi OS, which does run on
  all Raspberry Pi models, it is necessary to use the OS distribution's
  versions of Python and its packages. Instructions are provided here:
  https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallPiLinux. 

Python Requirements
-----------------------

GSAS-II requires a standard Python interpreter to be installed, as
well as several separately-developed packages that are not supplied
with Python, as are described below.
While for some packages, we have not seen much dependence on
versions, for others we do find significant differences; this is also
discussed further below. The GSAS-II GUI will warn about Python and
packages versions that are believed to be problematic,
as defined in variable :attr:`GSASIIdataGUI.versionDict`,
but for new installations we are currently recommending the following
interpreter/package versions: 

 * Python 3.10 is recommended, but 3.7 or later is fine. 
 * wxPython 4.2 or later is recommended, but with Python <=3.9 any
   wx4.x version should be OK. However,
   expect problems with Py>=3.10 and anything older than wx4.2.0.
 * NumPy 1.23 recommended, but anything from 1.17 on is likely fine
 * matplotlib 3.6 is recommended, but 3.4 or later is preferred. 
 * pyOpenGL: no version-related problems have been seen.
 * SciPy: no version-related problems have been seen, but in at least one
   case multiple imports are tried to account for where function
   names have changed. 

For more details on problems noted with specific versions of Python
and Python packages, see comments below and details here:
:attr:`GSASIIdataGUI.versionDict`,
   
Note that GSAS-II is being developed using Python 3.7, 3.9 and
3.10. No testing has yet been done with Python 3.11.  We are no longer
supporting Python 2.7 and <=3.6, and strongly encourage that
systems running GSAS-II under these older Python versions reinstall
Python typically via a new GSAS-II installation. 

There are a number of ways to install Python plus the packages
needed by GSAS-II. We had been creating installers utilizing the Anaconda
Python (https://www.anaconda.com/)
distribution, but no longer provides modern versions that include
wxPython or subversion, so we are transitioning to the
community-supported conda-forge library of packages and the miniforge
distribution in order to continue using the conda package manager and
installation process. This approach is available for nearly all supported
platforms (see below.)

Alternately, on MacOS, homebrew can be used for Python and most
needed packages, while on Linux, the native package installers
(apt-get or yum, etc.) offer the same. Any packages not provided in
that fashion can be installed with Python's pip mechanism. 
Other alternative Python packaging methods include Enthought Inc.'s Canopy and
Python(x,y), see here:
https://www.python.org/download/alternatives/. We are no longer using
any of them and are unsure of how well they will function, but in
theory any mechanism that supplies an internally compatible Python
with the packages required by GSAS-II should work fine. 

Package requirements depend on how GSAS-II is run. More packages are
required for GUI use and and stil others may optionally be used, as described
in the following section. A server that will only run GSAS-II
via the scripting interface, will need far fewer packages, as is
discussed in the section immediately following that. 


GUI Requirements
----------------

When using the GSAS-II graphical user interface (GUI), the following
Python extension packages are required:

* wxPython (http://wxpython.org/docs/api/). Note that GSAS-II has been
  tested with various wxPython versions over the years.  We encourage
  use of 4.x with Python 3.x, but with Py>=3.10 you must use
  wxPython 4.2.0 or later.
* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/),
* matplotlib (http://matplotlib.org/contents.html)  and
* PyOpenGL (http://pyopengl.sourceforge.net/documentation). 

GSAS-II will not start if the above packages are not available. In
addition, several Python packages are referenced in sections of the
GUI code, but are not required. If these packages are not present, warning
messages may be generated if they would be needed, or menu items may
be omitted, but the vast bulk of GSAS-II will function normally. These
optional packages are:

* Pillow (https://pillow.readthedocs.org) or PIL (http://www.pythonware.com/products/pil/). This is used to read and save certain types of images.
* h5py is the HDF5 interface and hdf5 is the support package. These
  packages are (not surprisingly) required
  to import images from HDF5 files. If these libraries are not present,
  the HDF5 importer(s) will not appear in the import menu and a
  warning message appears on GSAS-II startup. 
* imageio is used to make movies. This is optional and is offered for plotting
  superspace (modulated) structures. 
* requests: this package simplifies http access
  (https://requests.readthedocs.io/). It is used for access to
  webpages such as ISODISTORT and for some internal software downloads.
* win32com (windows only): this module is
  used to install GSAS-II on windows machines. GSAS-II can be used on
  Windows without this, but the installation will offer less
  integration into Windows. Conda provides this under the name pywin32.
* conda: the conda package allows access to package installation,
  etc. features from  inside Python. It is not required but is helpful
  to have, as it allows GSAS-II to install some packages that are not
  supplied initially. The conda package is included by default in
  the base miniconda and anaconda installations, but if you create an
  environment for GSAS-II 
  (`conda create -n <env> package-list...`), it will not be added
  to that environment unless you request it specifically.  

The following conda package is used where possible in GSAS-II but it provides a
command-line tool rather than a Python package.
  
* svn: the GSAS-II code utilizes the subversion
  program for software installation and updates. GSAS-II can be manually
  installed without it, but updates will also need to be done
  manually. Thus, GSAS-II works much better when
  subversion is available. The Anaconda distribution had provided
  subversion in a package named svn, but this is so no longer being updated. With
  the conda-forge repository we now use, it is only available for
  Linux (where it really is not needed since it is easy to install
  there) and the package is named subversion. (For the Mac the
  supplied subversion package lacks the ability to reach the GSAS-II
  repository via the internet and is thus not used.) 
  For MacOS and Windows, the GSAS-II gsas2full self-installer now
  provides binaries for the svn program.
  
*Conda command*:
  Should you wish to install Python and the desired packages yourself,
  this is certainly possible. For Linux, ``apt`` or ``yum`` is an option, as is
  homebrew. Homebrew is a good option on MacOS. However, we recommend  use
  of the miniconda or mambaconda self installers from
  conda-forge. Here is a typical conda command used to install a GSAS-II compatible
  Python interpreter on Linux after
  miniconda/miniforge/mambaforge/anaconda has been installed::

    conda install python=3.10 wxpython numpy scipy matplotlib pyopengl pillow h5py imageio subversion requests -c conda-forge
    
  or to put a Python configured for GSAS-II into a separate conda
  environment (below named ``g2python``, but any name can be used), use
  command::

    conda create -n g2python python=3.10 wxpython numpy scipy matplotlib pyopengl  pillow h5py imageio conda subversion requests -c conda-forge 

 For Windows/Mac/Raspberry Pi, omit subversion from the previous
 commands are::

    conda install python=3.10 wxpython numpy scipy matplotlib pyopengl pillow h5py imageio requests -c conda-forge
   
 and::

    conda create -n g2python python=3.10 wxpython numpy scipy matplotlib pyopengl  pillow h5py imageio conda requests -c conda-forge 

Before starting GSAS-II under conda remember to activate using:
``<path>\Scripts\activate``  (windows);
``source <path>/bin/activate`` (Mac/Linux),
or when an environment is used, add that name, (such as ``g2python``),
such as 
``<path>\Scripts\activate g2python``  (windows);
``source <path>/bin/activate g2python`` (Mac/Linux),


Note that at present we are not suppling binaries for Python 3.11, but
we are not aware of any reason why GSAS-II will not run fine with
this. 
  
Scripting  Requirements
-----------------------

When using the GSAS-II scripting interface (:mod:`GSASIIscriptable`),
only the two commonly-used Python extension packages are required:

* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/).

Note that a few sections of the code optionally require matplotlib
(http://matplotlib.org/contents.html), Pillow 
(https://pillow.readthedocs.org) (or PIL,
http://www.pythonware.com/products/pil/), or h5py + hdf5, but none of
these are required to run scripts and the vast
majority of scripts will not need these packages.

Optional Python Packages
---------------------------

* Sphinx (https://www.sphinx-doc.org) is used to generate the
  documentation you are currently reading. Generation of this documentation
  is not generally something needed by users or even most code
  developers since the prepared documentation on
  https://gsas-ii.readthedocs.io is usually reasonably up to date.  

 * SCons (https://scons.org/) is used to compile the relatively small amount of
   Fortran code that is included with GSAS-II. Use of this is
   discussed in the next section of this chapter.

Required Binary Files
--------------------------------

As noted before, GSAS-II also requires that some code be compiled.
For the following platforms, binary images are provided:

  * Windows-10: 64-bit Intel-compatible processors. [Prefix `win_64_`\ ]
  * MacOS: Intel processors. [Prefix `mac_64_`\ ]
  * MacOS: ARM processors, aka Apple Silicon (M1, etc). [Prefix `mac_arm_`\ ]
  * Linux: 64-bit Intel-compatible processors. [Prefix `linux_64_`\ ]
  * Linux: ARM processors (64-bit and 32-bit Raspberry Pi OS only).
    [Prefixes `linux_arm32_` and `linux_arm64_`\ ]

Note that these binaries must match the major versions of both Python and
numpy; binaries for only a small number of combinations are provided.
A full list of what is available can be seen by looking at the
contents of the directory at web address
https://subversion.xray.aps.anl.gov/trac/pyGSAS/browser/Binaries,
noting that a subdirectory name will be `prefix`\ _p\ `X.X`\ _n\ `Y.Y` where
`prefix` is noted above and `X.X` is the Python version and `Y.Y` is the numpy
version.
Should one wish to run GSAS-II where binary files are not
supplied (such as 32-bit Windows or Linux) or with other combinations of
Python/NumPy, compilation will be need to be done by the user.
This will require the GNU Fortran (gfortran)
compiler (https://gcc.gnu.org/fortran/) as well as the Python SCons
package. General instructions are provided for Linux: 
https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallLinux#CompilingFortranCode;
Windows: https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/CompilingWindows
and MacOS:
https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallMacHardWay,
but these may be out of date or require adaptation. 

Supported Externally-Developed Software
----------------------------------------------------

GSAS-II provides interfaces to use a number of programs developed by
others. Some are included with GSAS-II and others must be installed
separately. When these programs are accessed, citation
information is provided as we hope that users will recognize the
contribution made by the authors of these programs and will honor those
efforts by citing that work in addition to GSAS-II. 

GSAS-II includes copies of the following programs. No additional steps
beyond a standard installation are needed to access their functionality.

  **DIFFaX**
    Simulate layered structures with faulting. https://www.public.asu.edu/~mtreacy/DIFFaX.html
    
  **PyCifRW**
    A software library that reads and writes files using the IUCr's 
    Crystallographic Information Framework (CIF).
    https://bitbucket.org/jamesrhester/pycifrw. GSAS-II uses this to
    read data and structures from CIF files, 
    
  **Shapes**
    Derives the shapes of particles from small angle scattering data.
    
  **NIST FPA**
    Use Fundamental Parameters to determine GSAS-II profile function 

  **NIST*LATTICE**
    Searches for higher symmetry unit cells and possible relationships
    between unit cells. An API has been written and this will be
    integrated into the GSAS-II GUI. 

The following web services can also be accessed from computers that
have internet access. All software needed for this access is included
with GSAS-II.

  **Bilboa Crystallographic Server** (https://www.cryst.ehu.es):
    GSAS-II can directly access the Bilboa Crystallographic Server to
    utilize the k-SUBGROUPSMAG, k-SUBGROUPS and PseudoLattice web utilities for
    computation of space group subgroups, color (magnetic) subgroups &
    lattice search.

  **BYU ISOTROPY Software Suite** (https://stokes.byu.edu/iso/isotropy.php):
    GSAS-II directly accesses capabilities in the ISOTROPY Software
    Suite from Brigham Young University for representational analysis
    and magnetism analysis.  

At the request of the program authors, other programs that can be
accessed within GSAS-II are not included
as part of the GSAS-II distribution and must be installed separately:

  **Dysnomia**
    Computes enhanced Fourier maps with Maximum Entropy estimated
    extension of the reflection sphere. See https://jp-minerals.org/dysnomia/en/.

  **RMCProfile**
    Provides large-box PDF & S(Q) fitting. The GSAS-II interface was originally
    written for use with release 6.7.7 of RMCProfile, but updates have
    been made for compatible with 6.7.9 as well.
    RMCProfile must be downloaded by the user from
    http://rmcprofile.org/Downloads or
    https://rmcprofile.pages.ornl.gov/nav_pages/download/

  **fullrmc**
    A modern software framework for large-box PDF & S(Q) fitting. Note
    that the GSAS-II implementation is not compatible with the last
    open-source version of fullrmc, but rather the version 5.0 must be
    used, which is distributed only as compiled versions and only for 64-bit
    Intel-compatible processors running Windows, Linux and
    MacOS. Download this as a single executable from website
    https://github.com/bachiraoun/fullrmc/tree/master/standalones. GSAS-II
    will offer to install this software into the binary directory when the fullrmc
    option is selected on the Phase/RMC tab. 

  **PDFfit2**
    For small-box fitting of PDFs; see
    https://github.com/diffpy/diffpy.pdffit2#pdffit2. This code is no 
    longer supported by the authors, but is still quite useful. It can
    only be run with older Python versions. It is supplied within
    GSAS-II for Windows with Python 3.7-3.9 and for MacOS only with Python 3.7.
    When running GSAS-II with later versions of Python, as is strongly
    encouraged, it is best to install a separate older Python
    interpreter specifically for PDFfit2. When GSAS-II is run from a
    Python installation that includes the conda package manager (the
    usual installation practice), the GUI will offer an option to
    install PDFfit2 via a separate Python 3.7 environment when the
    PDFfit2 option is selected on the Phase/RMC tab. 
