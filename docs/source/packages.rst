Required packages
====================

GSAS-II requires a standard Python interpreter to be installed, as
well as several separately-developed packages. GSAS-II is being
developed using Python 2.7, 3.6 and 3.7. At this point we think that 
most sections of the code have been exercised in Python 2 and 3,
but  bugs are still expected (please report them). Our
goal is to keep the code compliant with both Python 2.7 and 3.x for
the immediate future. 

Note that the packages listed below are not distributed as part of the Python standard
library. We use the free Anaconda Python (https://www.anaconda.com/)
distribution (and provide installers based on that), but there are
many other fine distributions, such as Enthought Inc.'s Canopy and
Python(x,y), see here: https://www.python.org/download/alternatives/. 
We do some testing using the older Enthought Python Distribution
(EPD); this is known to have some problems with reading CIFs and
encourage updating from that.

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
  GSAS-II to download updates to our code.

Scripting  Requirements
-----------------------

When using the GSAS-II scripting interface (:mod:`GSASIIscriptable`),
only the following Python extension packages are required:

* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/).

Note that some sections of the code may require matplotlib (http://matplotlib.org/contents.html), Pillow
(https://pillow.readthedocs.org) (or PIL,
http://www.pythonware.com/products/pil/), or h5py + hdf5 to function
but scripts will load and run without these. 


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
    version 5.0 or later. 

  **Dysnomia**
    Computes enhanced Fourier maps with Maximum Entropy estimated
    extension of reflection sphere
