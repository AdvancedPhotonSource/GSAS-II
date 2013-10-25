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
  GSASIIscripts.rst
  exports.rst
  imports.rst

*Required packages*
--------------------

Note that GSAS-II requires the Python extension packages 

* wxPython (http://wxpython.org/docs/api/), 
* NumPy (http://docs.scipy.org/doc/numpy/reference/), 
* SciPy (http://docs.scipy.org/doc/scipy/reference/),
* matplotlib (http://matplotlib.org/contents.html) and 
* PyOpenGL (http://pyopengl.sourceforge.net/documentation)

These are not distributed
as part of the Python standard library and must be obtained separately
(or in a bundled Python package such as the Enthought Python
Distribution/Canopy). 
The PyOpenGL package is installed into Python by GSAS-II if not found,
so it does not need to be included in the Python bundle.
