GSAS-II Requirements, Python Packages and External Software
==============================================================

Supported Platforms
--------------------------------

It should be possible to run GSAS-II on any computer where Python 3.7+ and
the appropriate required packages are available, as discussed below,
but GSAS-II also requires that some code must be compiled.
For the following platforms, binary images for this compiled code are
currently provided:

  * Windows-10: 64-bit Intel-compatible processors 
  * MacOS: Intel processors 
  * MacOS: ARM processors, aka Apple Silicon (M1, etc) 
  * Linux: 64-bit Intel-compatible processors
  * Linux: ARM processors (64-bit Raspberry Pi OS only)

Details for GSAS-II use on these specific platforms follows below:

* **Windows**: self-Installation kits are provided for 
  64-bit Windows-10 and -11
  `here
  <https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/latest>`_.
  Less testing has been done with
  Windows-11, but both appear to working interchangeably with respect
  to GSAS-II. 

  In theory it should be possible to run GSAS-II on older versions of
  Windows, including 32-bit OS versions, but no current installation kit
  can be provided. Installing GSAS-II will require locating a
  compatible version (or compiling) Python and the required
  packages. It may be necessary to recompile the GSAS-II binaries. 

* **MacOS**: GSAS-II can run natively on Intel (or ARM ("M1"-"M3" aka "Apple
  Silicon") processors with relatively current versions of MacOS, with
  self-installers that can be run from the command-line available for download `here
  <https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/latest>`_.
  The Intel version will run on both types of Mac processors, but the
  native ARM versions offer  the highest GSAS-II performance we see on
  any platform.
  
  It appears that this installer can be used with MacOS versions 11.0
  and later.  Macs older than Catalina (10.15) will likely require older
  distributions of Python.  

* **Intel Linux**: Note that GSAS-II does not get a lot of testing
  in Linux by us, but is used fairly widely on this platform
  nonetheless.  We provide an installer `here
  <https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/latest>`_
  that includes Python and
  needed packages for Intel-compatible Linuxes, but compatibility with
  older and very new versions of Linux can sometimes be tricky as
  compatibility libraries may be needed -- not always easy to do. It may be
  better to use your Linux distribution's versions of Python and
  packages (typically done with a software tool such as apt or yum or
  pip. See
  https://advancedphotonsource.github.io/GSAS-II-tutorials/install-pip.html
  for more information.  

* **Non-Intel Linux**:
  Will GSAS-II run on Linux with other types of CPUs? That will mostly
  depend on support for Python and wxPython on that CPU. If those can
  be used, you can likely build the GSAS-II binaries with gcc &
  gfortran. Expect to modify the SConstruct file. 
  
  **Raspberry Pi** (ARM) Linux: GSAS-II has been installed on both 32-bit
  and the 64-bit version of the Raspberry Pi OS (formerly
  called Raspbian) and compiled binaries are provided at present for
  both, but the 32-bit support may not continue. It is expected that
  these binaries will also function on Ubuntu Linux for Raspberry Pi,
  but this has not been tried. 
  The performance of GSAS-II on a Raspberry Pi is not blindingly fast,
  but one can indeed run GSAS-II on a motherboard that costs only $15
  (perhaps even one that costs $5) and uses <5 Watts! 

  Note that the 64-bit OS is preferred on the models where it can be run
  (currently including models 3A+, 3B, 3B+, 4, 400, CM3, CM3+, CM4,
  and Zero 2 W) .  With the 32-bit Raspberry Pi OS, which does run on
  all Raspberry Pi models, it is necessary to use the OS distribution's
  versions of Python and its packages, `see here   for more information
  <https://advancedphotonsource.github.io/GSAS-II-tutorials/install-pip.html>`_.
  With
  64-bit Pi OS it may be possible for us to provide a GSAS2FULL installer
  (which will need to include a custom-supplied wxPython wheel, since
  that is not available in conda-forge) or else pip must be used to
  download and build wx. Please let us know if you are intending to
  use GSAS-II on a Raspberry Pi for a classroom, etc and would need
  this. 

Version Control
-----------------------
The master version of the source code for GSAS-II resides on
GitHub at URL and the git 
version control system (VCS) is usually used to install the files needed by GSAS-II. When
GSAS-II is installed in this manner, the software can be easily
updated, as git commands can download only the changed sections of files
that need to be updated. It is likewise possible to use git to regress
to an older version of GSAS-II, though there are some limitations on
how far back older versions of GSAS-II will be with current versions
of Python. While git is not required for use of GSAS-II, special
procedures must be used to install GSAS-II without it and once
installed without git, updates of GSAS-II must be done manually. 

We are currently in a transition period to GitHub from
a previous subversion server. As we migrate to GitHub, updates will be
made in parallel to both servers.

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

 * Python 3.11 is recommended. GSAS-II should run with any Python
   version from 3.7 or later, but you will need to locate (from the
   old subversion server) or locate binaries to match that Python version. 
 * wxPython 4.2 or later is recommended, but with Python <=3.9 any
   wx4.x version should be OK. However, there may be problems with
   newer sections of the GUI with wx <4.0.
 * NumPy 1.26 recommended, but anything from 1.17 on is likely fine,
   but does need to approximately match the GSAS-II binaries. 
 * matplotlib 3.6 is recommended, but 3.4 or later is preferred. 
 * pyOpenGL: no version-related problems have been seen.
 * SciPy: no version-related problems have been seen, but in at least one
   case multiple imports are tried to account for where function
   names have changed. 

For more details on problems noted with specific versions of Python
and Python packages, see comments below and details here:
:attr:`GSASIIdataGUI.versionDict`,
   
Note that GSAS-II is currently being developed using Python 3.11. We
are seeing compilation problems with Python 3.12 that will be addressed
later via the build migration to meson.  We are no longer
supporting Python 2.7 and <=3.6, and strongly encourage that
systems running GSAS-II under these older Python versions reinstall
Python. Typically this is done by reinstalling GSAS-II from a current self-installer. 

There are a number of ways to install Python plus the packages
needed by GSAS-II. See
https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html
and links therein for a discussion of installation. 

Python package requirements depend on how GSAS-II will be run, as will be
discussed in the next section. In order to run
the GUI for GSAS-II, a much larger number of packages are
required. Several more packages are optional, but some functionally will
not be available without those optional packages.
Far fewer packages are required to run GSAS-II on a
compute server via the scripting interface
and without a GUI.

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

* gitpython: (https://gitpython.readthedocs.io and
  https://github.com/gitpython-developers/GitPython). This 
  this package provides a bridge between the git version control
  system and Python. It is required for the standard GSAS-II
  installation process and for GSAS-II to update itself from GitHub.
  If your computer does not already have git in the path, also include
  the git package to obtain that binary (if you are not sure, it does
  not hurt to do this anyway).
* requests: this package simplifies http access
  (https://requests.readthedocs.io/). It is used for access to
  webpages such as ISODISTORT and for some internal software
  downloads. It is required for support of git updating and installation.  
* Pillow (https://pillow.readthedocs.org) or PIL (http://www.pythonware.com/products/pil/). This is used to read and save certain types of images.
* h5py is the HDF5 interface and hdf5 is the support package. These
  packages are (not surprisingly) required
  to import images from HDF5 files. If these libraries are not present,
  the HDF5 importer(s) will not appear in the import menu and a
  warning message appears on GSAS-II startup. 
* imageio is used to make movies. This is optional and is offered for plotting
  superspace (modulated) structures. 
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
  
*Conda command*:
  Should you wish to install Python and the desired packages yourself,
  this is certainly possible. For Linux, ``apt`` or ``yum`` is an option, as is
  homebrew. Homebrew is a good option on MacOS. However, we recommend  use
  of the miniforge self-installers from
  conda-forge. Here is a typical conda command used to install a GSAS-II compatible
  Python interpreter after miniforge has been installed::

       conda install python=3.11  numpy=1.26 wxpython scipy matplotlib pyopengl pillow h5py imageio requests git gitpython -c conda-forge

  for development environments, it is useful to have build and
  debugging tools available, so here is a more extensive list of
  useful packages::
    
     conda create -n py311 python=3.11 numpy=1.26 matplotlib scipy wxpython  pyopengl imageio h5py hdf5 pillow requests ipython conda spyder-kernels scons sphinx sphinx-rtd-theme jupyter git gitpython -c conda-forge

To find out what packages have been directly installed in a conda
environment this command can be used::
  
  conda env export --from-history -n <env>

.. _ScriptingRequirements:

  
Scripting Requirements
-----------------------

The GSAS-II scripting interface (:mod:`GSASIIscriptable`) will not
run without two Python extension packages:

* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/).

These fortunately are common and are easy to install. There are
some relatively minor scripting capabilities that will only run when a few
additional packages are installed:
  
* matplotlib (http://matplotlib.org/contents.html),
* Pillow (https://pillow.readthedocs.org) and/or
* h5py and hdf5 

but none of these are required to run scripts and the vast
majority of scripts will not need these packages.

**Installing a minimal Python configuration**:

There are many ways to install a minimal Python configuration.
Below, I show some example commands used to install using the 
the free miniconda installer from Anaconda, Inc., but I now tend to
use the Conda-Forge miniforge distributions instead. 
However, there are also plenty of  other ways to install Python, Numpy
and Scipy, depending on if they will be used on Linux, Windows and MacOS.
For Linux, the standard Linux distributions provide these using
``yum`` or ``apt-get`` etc., but these often supply package versions
that are so new that they probably have not been tested with GSAS-II.

.. code-block::  bash

    bash ~/Downloads/Miniconda3-latest-<platform>-x86_64.sh -b -p /loc/pyg2script
    source /loc/pyg2script/bin/activate
    conda install numpy scipy matplotlib pillow h5py hdf5

Some discussion on these commands follows:

* the 1st command (bash) assumes that the appropriate version of Miniconda has been downloaded from https://docs.conda.io/en/latest/miniconda.html and ``/loc/pyg2script`` is where I have selected for python to be installed. You might want to use something like ``~/pyg2script``.
* the 2nd command (source) is needed to access Python with miniconda. 
* the 3rd command (conda) installs all possible packages that might be
  used by scripting, but matplotlib, pillow, and hdf5 are not commonly
  needed and could be omitted.

Once Python has been installed and is in the path, use these commands to install GSAS-II:

.. code-block::  bash

    git clone https://github.com/AdvancedPhotonSource/GSAS-II.git /loc/GSAS-II
    python /loc/GSAS-II/GSASII/GSASIIscriptable.py

Notes on these commands:

* the 1st command (git) is used to download the GSAS-II software. ``/loc/GSASII`` is the location where I decided to install the software. You can select something different. 
* the 2nd command (python) is used to invoke GSAS-II scriptable for the first time, which is needed to load the binary files from the server.


Optional Python Packages
---------------------------

* Sphinx (https://www.sphinx-doc.org) is used to generate the
  documentation you are currently reading. Generation of this documentation
  is not generally something needed by users or even most code
  developers, since the prepared documentation on
  https://gsas-ii.readthedocs.io is usually reasonably up to date.  

* SCons (https://scons.org/) is used to compile the relatively small amount of
   Fortran code that is included with GSAS-II. Use of this is
   discussed in the next section of this chapter.

Required Binary Files
--------------------------------

As noted before, GSAS-II also requires that some code be compiled.
For the following platforms, binary images are provided at
https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/latest
for Python 3.11 and NumPy 1.26:

  * Windows-10: 64-bit Intel-compatible processors.
  * MacOS: Intel processors.
  * MacOS: ARM processors, aka Apple Silicon (M1, etc). 
  * Linux: 64-bit Intel-compatible processors. 
  * Linux: ARM processors (64-bit and 32-bit Raspberry Pi OS and
    Ubuntu for Raspberry Pi).

Note that these binaries must match the major versions of both Python and
NumPy; 
Should one wish to run GSAS-II where binary files are not
supplied (such as 32-bit Windows or Linux) or with other combinations of
Python/NumPy, compilation will be need to be done by the user. See
the `compilation information <https://advancedphotonsource.github.io/GSAS-II-tutorials/compile.html>`_ for more information. 

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

  **pybaselines**
   Determines a background for a powder pattern in the "autobackground"
   option. See https://pybaselines.readthedocs.io for more
   information. 
    
The following web services can also be accessed from computers that
have internet access. All software needed for this access is included
with GSAS-II.

  **Bilbao Crystallographic Server** (https://www.cryst.ehu.es):
    GSAS-II can directly access the Bilbao Crystallographic Server to
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
    longer being updated by the authors, but is still quite useful.
    It is supplied within GSAS-II for Python 3.7.
    It is likely best to install a separate Python
    interpreter specifically for PDFfit2. When GSAS-II is run from a
    Python installation that includes the conda package manager (the
    usual installation practice), the GUI will offer an option to
    install PDFfit2 via a separate environment when the
    PDFfit2 option is selected on the Phase/RMC tab. 
