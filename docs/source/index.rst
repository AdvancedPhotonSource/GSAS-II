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
  GSASIIscripts.rst
  exports.rst
  imports.rst

*Required packages*
--------------------

Note that GSAS-II requires the Python extension packages 

* wxPython (http://wxpython.org/docs/api/), 
* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/),
* matplotlib (http://matplotlib.org/contents.html)  and
* PyOpenGL (http://pyopengl.sourceforge.net/documentation)

Two packages are used by some parts of the code, but are not required:

* PIL (http://www.pythonware.com/products/pil/) or Pillow (https://pillow.readthedocs.org). This is used to save
  and read certain types of images.
* h5py is the HDF5 support package. This is (not surprisingly) required
  to import images from HDF5 files. If this library is not present,
  the HDF5 importer(s) will not appear in the import menu and a
  warning message appears on GSAS-II startup. 

Note that the packages listed above are not distributed as part of the Python standard
library and must be obtained separately (or in a bundled Python
package such as the Enthought Inc.'s Canopy or
Continuum.io's anaconda; we also use the older Enthought Python
Distribution).  
One exception is the PyOpenGL package. This will be installed into
Python by GSAS-II if not found, so it does not need to be included in
the Python bundle, but the setuptools package
(https://pythonhosted.org/setuptools/) is needed by GSAS-II to install
PyOpenGL.
