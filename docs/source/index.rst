.. This lists all the sections of the documentation. Some of the below
.. reference a single file. Others reference multiple files. 

GSAS-II Developer's Documentation
=================================

.. toctree::

  GSASII.rst
  GSASIIobj.rst
  GSASIIutil.rst
  GSASIIGUIr.rst
  GSASIIGUI.rst
  GSASIIstruc.rst
  GSASIImapvars.rst
  GSASIIimage.rst
  GSASIImath.rst
  GSASIIindex.rst
  GSASIIplot.rst
  GSASIIpwd.rst
  SAS.rst
  GSASIIscriptable.rst
  GSASIIscripts.rst
  imports.rst 
  exports.rst

*Required packages*
--------------------

GSAS-II requires a standard Python interpreter to be installed, as
well as several separately-developed packages. GSAS-II is being
developed using both Python 2.7 and Python 3.6, but some sections of
the code have not been exercised in Python 3 so bugs are to be
expected (please report them). Our
goal is to keep the code compliant with both Python 2.7 and 3.x for
the immediate future. 

Note that GSAS-II requires the Python extension packages 

* wxPython (http://wxpython.org/docs/api/). Note that GSAS-II has been tested with wxPython >=2.8, 3.0.x and 4.0.x
* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/),
* matplotlib (http://matplotlib.org/contents.html)  and
* PyOpenGL (http://pyopengl.sourceforge.net/documentation). Note: a copy of this is distributed with GSAS-II (at present) and will be installed if the Python setuptools package is present. 

Several packages are used by some parts of the code, but are not
required. If these packages are not present warning messages may be
generated when needed, but the vast bulk of GSAS-II will function normally. 

* PIL (http://www.pythonware.com/products/pil/) or Pillow (https://pillow.readthedocs.org). This is used to save
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

Note that the packages listed above are not distributed as part of the Python standard
library. We use the free Anaconda Python (https://www.anaconda.com/)
distribution (and provide installers based on that), but there are
many other fine distributions, such as Enthought Inc.'s Canopy and
Python(x,y), see here: https://www.python.org/download/alternatives/. 
We do some testing using the older Enthought Python Distribution
(EPD); this is known to have some problems with reading CIFs and
encourage updating from that. 
