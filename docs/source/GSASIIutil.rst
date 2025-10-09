*GSAS-II Utility Modules*
=========================

---------------------------------
*GSASIIpath: locations & updates*
---------------------------------

Routines for dealing with file locations, etc.

Determines the location of the compiled (.pyd or .so) libraries.

Interfaces with git and gets version number.

Accesses configuration options, as defined in config.py

GSASIIpath Classes & Routines
------------------------------------

.. automodule:: GSASII.GSASIIpath
    :members: 

-------------------------------------------
*config_example.py: Configuration options*
-------------------------------------------

Configuration variables
------------------------------------

.. automodule:: GSASII.config_example
    :members: 

-----------------------------------------
*GSASIIElem: functions for element types*
-----------------------------------------

GSASIIElem Routines
------------------------------------

.. automodule:: GSASII.GSASIIElem
    :members: 

-----------------------------------------
*GSASIIlattice: Unit Cell Computations*
-----------------------------------------

Performs lattice-related computations

Note that as used here 
*G* is the reciprocal lattice tensor, and *g* is its inverse,
:math:`G = g^{-1}`, where 

  .. math::

   g = \left( \begin{matrix}
   a^2 & a b\cos\gamma & a c\cos\beta \\
   a b\cos\gamma & b^2 & b c \cos\alpha \\
   a c\cos\beta &  b c \cos\alpha & c^2
   \end{matrix}\right)

The "*A* tensor" terms are defined as
:math:`A = (\begin{matrix} G_{11} & G_{22} & G_{33} & 2G_{12} & 2G_{13} & 2G_{23}\end{matrix})` and *A* can be used in this fashion:
:math:`d^* = \sqrt {A_0 h^2 + A_1 k^2 + A_2 l^2 + A_3 hk + A_4 hl + A_5 kl}`, where
*d* is the d-spacing, and :math:`d^*` is the reciprocal lattice spacing, 
:math:`Q = 2 \pi d^* = 2 \pi / d`. 
Note that GSAS-II variables ``p::Ai`` (i = 0, 1,... 5) and ``p`` is a phase number are 
used for the *Ai* values. See :func:`~GSASIIlattice.A2cell`, :func:`~GSASIIlattice.cell2A` for interconversion between A and 
unit cell parameters; :func:`~GSASIIlattice.cell2Gmat` :func:`~GSASIIlattice.Gmat2cell` for G and cell parameters. 

When the hydrostatic/elastic strain coefficients (*Dij*, :math:`D_{ij}`) are used, they are added to the 
*A* tensor terms (Ai, :math:`A_{i}`) so that A is redefined 
:math:`A = (\begin{matrix} A_{0} + D_{11} & A_{1} + D_{22} & A_{2} + D_{33} & A_{3} + D_{12} & A_{4} + D_{13} & A_{5} + D_{23}\end{matrix})`. See :func:`~GSASIIlattice.cellDijFill`. 
Note that GSAS-II variables ``p:h:Dij`` (i,j = 1, 2, 3) and ``p`` is a phase number 
and ``h`` a histogram number are used for the *Dij* values.

GSASIIlattice Classes & Routines
------------------------------------

.. automodule:: GSASII.GSASIIlattice
    :members: 

       
-----------------------------------------
*GSASIIspc: Space Group Computations*
-----------------------------------------

Space group interpretation routines. Note that space group information is
stored in a :ref:`Space Group (SGData)<SGData_table>` object.

GSASIIspc Classes & Routines
------------------------------------

.. automodule:: GSASII.GSASIIspc
    :members: 

---------------------------------------------
*GSASIIfiles: data (non-GUI) I/O routines*
---------------------------------------------

Module with miscellaneous routines for input and output from files.

GSASIIfiles Classes & Routines
------------------------------------

.. automodule:: GSASII.GSASIIfiles
    :members: 

--------------------------------------------------
*GSASIImpsubs: routines used in multiprocessing*
--------------------------------------------------

GSASIImpsubs Classes & Routines
------------------------------------

.. automodule:: GSASII.GSASIImpsubs
    :members: 

---------------------------------------------------
*Module nistlat: NIST*LATTICE cell computations*
---------------------------------------------------

nistlat Classes & Routines
------------------------------------

.. automodule:: GSASII.nistlat
    :members: 

-----------------------------------------
*ReadMarCCDFrame: Read Mar Files*
-----------------------------------------

.. automodule:: GSASII.ReadMarCCDFrame
    :members: 

-----------------------------------------
*G2shapes: Compute SAS particle shapes*
-----------------------------------------

Program SHAPES from
"A New Algroithm for the Reconstruction of Protein Molecular Envelopes
from X-ray Solution Scattering Data", 
John Badger, J. Appl. Cryst. (2019) 52, 937-944.
(DOI: 10.1107/S1600576719009774) modified to run inside GSAS-II.

.. automodule:: GSASII.G2shapes
    :members:
       
--------------------------------------------
*tutorialIndex: index to GSAS-II tutorials*
--------------------------------------------

The file `tutorialIndex.py` provides information on what tutorials
are currently available in GSAS-II. The information on this is placed
in this file and all indices to tutorials are generated from this.

tutorialIndex Contents
------------------------------------

.. automodule:: GSASII.tutorialIndex
    :members: 

--------------------------------------------
*k_vector_search: k-vector utilities*
--------------------------------------------

This contains routines used to perform k-vector (magnetic scattering
lattice expansion) determination. 

k_vector_search Contents
------------------------------------

.. automodule:: GSASII.k_vector_search
    :members: 

--------------------------------------------
*pathHacking utilities*
--------------------------------------------

When GSAS-II is installed via Git (gsas2main for example, as opposed
to pixi), the location of GSAS-II is "outside" the Python file
structure. It then becomes necessary to set the Python `sys.path`
variable to find GSAS-II files. This practice is discouraged and
commands to do this are "exiled" to file `pathHacking.py`.

.. automodule:: GSASII.pathHacking
    :members: 

