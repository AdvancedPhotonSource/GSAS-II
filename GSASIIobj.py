# -*- coding: utf-8 -*-
#GSASIIobj - data objects for GSAS-II
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################

'''
*GSASIIobj: Data objects*
=========================

This module defines and/or documents the data structures used in GSAS-II, as well
as provides misc. support routines. 

Constraints Tree Item
----------------------

.. _Constraints_table:

.. index::
   single: Constraints object description
   single: Data object descriptions; Constraints

Constraints are stored in a dict, separated into groups.
Note that parameter are named in the following pattern,
p:h:<var>:n, where p is the phase number, h is the histogram number
<var> is a variable name and n is the parameter number.
If a parameter does not depend on a histogram or phase or is unnumbered, that
number is omitted.
Note that the contents of each dict item is a List where each element in the
list is a :ref:`constraint definition objects <Constraint_definitions_table>`.
The constraints in this form are converted in
:func:`GSASIIstrIO.ProcessConstraints` to the form used in :mod:`GSASIImapvars`

The keys in the Constraints dict are:

.. tabularcolumns:: |l|p{4.5in}|

==========  ====================================================
  key         explanation
==========  ====================================================
Hist        This specifies a list of constraints on 
            histogram-related parameters,
            which will be of form :h:<var>:n.
HAP         This specifies a list of constraints on parameters
            that are defined for every histogram in each phase
            and are of form p:h:<var>:n.            
Phase       This specifies a list of constraints on phase
            parameters,
            which will be of form p::<var>:n.
Global      This specifies a list of constraints on parameters
            that are not tied to a histogram or phase and
            are of form ::<var>:n
==========  ====================================================

.. _Constraint_definitions_table:

.. index::
   single: Constraint definition object description
   single: Data object descriptions; Constraint Definition

Each constraint is defined as an item in a list. Each constraint is of form::

[[<mult1>, <var1>], [<mult2>, <var2>],..., <fixedval>, <varyflag>, <constype>]

Where the variable pair list item containing two values [<mult>, <var>], where: 

  * <mult> is a multiplier for the constraint (float)
  * <var> a :class:`G2VarObj` object (previously a str variable name of form
      'p:h:name[:at]')

Note that the last three items in the list play a special role:

 * <fixedval> is the fixed value for a `constant equation` (``constype=c``)
   constraint or is None. For a `New variable` (``constype=f``) constraint,
   a variable name can be specified as a str (used for externally
   generated constraints)
 * <varyflag> is True or False for `New variable` (``constype=f``) constraints
   or is None. This will be implemented in the future to indicate if these variables
   should be refined.
 * <constype> is one of four letters, 'e', 'c', 'h', 'f' that determines the type of constraint:

    * 'e' defines a set of equivalent variables. Only the first variable is refined (if the
      appropriate refine flag is set) and and all other equivalent variables in the list
      are generated from that variable, using the appropriate multipliers. 
    * 'c' defines a constraint equation of form,
      :math:`m_1 \\times var_1 + m_2 \\times var_2 + ... = c`
    * 'h' defines a variable to hold (not vary). Any variable on this list is not varied,
      even if its refinement flag is set. Only one [mult,var] pair is allowed in a hold
      constraint and the mult value is ignored.
      This is of particular value when needing to hold one or more variables where a
      single flag controls a set of variables such as, coordinates, 
      the reciprocal metric tensor or anisotropic displacement parameter. 
    * 'f' defines a new variable (function) according to relationship 
      :math:`newvar = m_1 \\times var_1 + m_2 \\times var_2 + ...`

Covariance Tree Item
--------------------

.. _Covariance_table:

.. index::
   single: Covariance description
   single: Data object descriptions; Covariance

The Covariance tree item has results from the last least-squares run. They
are stored in a dict with these keys:

.. tabularcolumns:: |l|l|p{4in}|

=============  ===============  ====================================================
  key            sub-key        explanation
=============  ===============  ====================================================
newCellDict    \                dict with lattice parameters computed by
                                :func:`GSASIIstrMath.GetNewCellParms` (dict)
title          \                Name of gpx file(?) (str)
variables      \                Values for all N refined variables
                                (list of float values, length N,
                                ordered to match varyList)
sig            \                Uncertainty values for all N refined variables
                                (list of float values, length N,
                                ordered to match varyList)
varyList       \                List of directly refined variables 
                                (list of str values, length N)
newAtomDict    \                dict with atom position values computed in 
                                :func:`GSASIIstrMath.ApplyXYZshifts` (dict)
Rvals          \                R-factors, GOF, Marquardt value for last
                                refinement cycle (dict)
\              Nobs             Number of observed data points (int)
\              Rwp              overall weighted profile R-factor (%, float)
\              chisq            sum[w*(Iobs-Icalc)**2] for all data
                                note this is not the reduced chi squared (float)
\              lamMax           Marquardt value applied to Hessian diagonal
                                (float)
\              GOF              The goodness-of-fit, aka square root of
                                the reduced chi squared. (float)
covMatrix      \                The (NxN) covVariance matrix (np.array)
=============  ===============  ====================================================

Phase Tree Items
----------------

.. _Phase_table:

.. index::
   single: Phase object description
   single: Data object descriptions; Phase

Phase information is stored in the GSAS-II data tree as children of the
Phases item in a dict with keys:

.. tabularcolumns:: |l|l|p{4in}|

==========  ===============  ====================================================
  key         sub-key        explanation
==========  ===============  ====================================================
General         \            Overall information for the phase (dict)
  \         AtomPtrs         list of four locations to use to pull info
                             from the atom records (list)
  \         F000X            x-ray F(000) intensity (float)
  \         F000N            neutron F(000) intensity (float)
  \         Mydir            directory of current .gpx file (str)
  \         MCSA controls    Monte Carlo-Simulated Annealing controls (dict)
  \         Cell             List with 8 items: cell refinement flag (bool)
                             a, b, c, (Angstrom, float)
                             alpha, beta & gamma (degrees, float)
                             volume (A^3, float)
  \         Type             'nuclear' or 'macromolecular' for now (str)
  \         Map              dict of map parameters
  \         SH Texture       dict of spherical harmonic preferred orientation 
                             parameters
  \         Isotope          dict of isotopes for each atom type
  \         Isotopes         dict of scattering lengths for each isotope
                             combination for each element in phase  
  \         Name             phase name (str)
  \         SGData           Space group details as a :ref:`space group (SGData) object <SGData_table>`
                             as defined in :func:`GSASIIspc.SpcGroup`.
  \         Pawley neg wt    Restraint value for negative Pawley intensities
                             (float)
  \         Flip             dict of Charge flip controls 
  \         Data plot type   data plot type ('Mustrain', 'Size' or 
                             'Preferred orientation') for powder data (str)
  \         Mass             Mass of unit cell contents in g/mol
  \         POhkl            March-Dollase preferred orientation direction
  \         Z                dict of atomic numbers for each atom type 
  \         vdWRadii         dict of van der Waals radii for each atom type 
  \         Color            Colors for atoms (list of (r,b,g) triplets)
  \         AtomTypes        List of atom types
  \         AtomMass         List of masses for atoms
  \         doPawley         Flag for Pawley intensity extraction (bool)
  \         NoAtoms          Number of atoms per unit cell of each type (dict)
  \         Pawley dmin      maximum Q (as d-space) to use for Pawley 
                             extraction (float)
  \         BondRadii        Default radius for each atom used to compute 
                             interatomic distances (list of floats)
  \         AngleRadii       Default radius for each atom used to compute 
                             interatomic angles (list of floats)
  \         DisAglCtls       Dict with distance/angle search controls,
                             which has keys 'Name', 'AtomTypes',
                             'BondRadii', 'AngleRadii' which are as above
                             except are possibly edited. Also contains
                             'Factors', which is a 2 element list with
                             a multiplier for bond and angle search range
                             [typically (0.85,0.85)].
ranId           \            unique random number Id for phase (int)
pId             \            Phase Id number for current project (int).
Atoms           \            Atoms in phase as a list of lists. The outer list
                             is for each atom, the inner list contains varying
                             items depending on the type of phase, see
                             the :ref:`Atom Records <Atoms_table>` description.
                             (list of lists)
Drawing         \            Display parameters (dict)
\           ballScale        Size of spheres in ball-and-stick display (float)
\           bondList         dict with bonds
\           contourLevel     map contour level in e/A^3 (float)
\           showABC          Flag to show view point triplet (bool). True=show.
\           viewDir          cartesian viewing direction (np.array with three
                             elements)
\           Zclip            clipping distance in A (float)
\           backColor        background for plot as and R,G,B triplet
                             (default = [0, 0, 0], black).
                             (list with three atoms)
\           selectedAtoms    List of selected atoms (list of int values)
\           showRigidBodies  Flag to highlight rigid body placement
\           sizeH            Size ratio for H atoms (float) 
\           bondRadius       Size of binds in A (float)
\           atomPtrs         positions of x, type, site sym, ADP flag in Draw Atoms (list)
\           viewPoint        list of lists. First item in list is [x,y,z]
                             in fractional coordinates for the center of
                             the plot. Second item list of previous & current 
                             atom number viewed (may be [0,0])
\           showHydrogen     Flag to control plotting of H atoms.
\           unitCellBox      Flag to control display of the unit cell.
\           ellipseProb      Probability limit for display of thermal
                             ellipsoids in % (float).
\           vdwScale         Multiplier of van der Waals radius for
                             display of vdW spheres. 
\           Atoms            A list of lists with an entry for each atom
                             that is plotted.
\           Zstep            Step to de/increase Z-clip (float)
\           Quaternion       Viewing quaternion (4 element np.array)
\           radiusFactor     Distance ratio for searching for bonds. ? Bonds
                             are located that are within r(Ra+Rb) and (Ra+Rb)/r
                             where Ra and Rb are the atomic radii.
\           oldxy            previous view point (list with two floats)
\           cameraPos        Viewing position in A for plot (float)
\           depthFog         True if use depthFog on plot - set currently as False (bool)
RBModels        \            Rigid body assignments (note Rigid body definitions
                             are stored in their own main top-level tree entry.)
Pawley ref      \            Pawley reflections
Histograms      \            A dict of dicts. The key for the outer dict is
                             the histograms tied to this phase. The inner
                             dict contains the combined phase/histogram
                             parameters for items such as scale factors,
                             size and strain parameters. (dict)
MCSA            \            Monte-Carlo simulated annealing parameters (dict)
\           
==========  ===============  ====================================================

Rigid Body Objects
------------------

.. _RBData_table:

.. index::
   single: Rigid Body Data description
   single: Data object descriptions; Rigid Body Data
   
Rigid body descriptions are available for two types of rigid bodies: 'Vector' 
and 'Residue'. Vector rigid bodies are developed by a sequence of translations each
with a refinable magnitude and Residue rigid bodies are described as Cartesian coordinates
with defined refinable torsion angles.

.. tabularcolumns:: |l|l|p{4in}|

==========  ===============  ====================================================
  key         sub-key        explanation
==========  ===============  ====================================================
Vector      RBId             vector rigid bodies (dict of dict)
\           AtInfo           Drad, Color: atom drawing radius & color for each atom type (dict)
\           RBname           Name assigned by user to rigid body (str)
\           VectMag          vector magnitudes in A (list)
\           rbXYZ            Cartesian coordinates for Vector rigid body (list of 3 float)
\           rbRef            3 assigned reference atom nos. in rigid body for origin 
                             definition, use center of atoms flag (list of 3 int & 1 bool)
\           VectRef          refinement flags for VectMag values (list of bool)
\           rbTypes          Atom types for each atom in rigid body (list of str)
\           rbVect           Cartesian vectors for each translation used to build rigid body (list of lists)
\           useCount         Number of times rigid body is used in any structure (int)
Residue     RBId             residue rigid bodies (dict of dict)
\           AtInfo           Drad, Color: atom drawing radius & color for each atom type(dict)
\           RBname           Name assigned by user to rigid body (str)
\           rbXYZ            Cartesian coordinates for Residue rigid body (list of 3 float)
\           rbTypes          Atom types for each atom in rigid body (list of str)
\           atNames          Names of each atom in rigid body (e.g. C1,N2...) (list of str)
\           rbRef            3 assigned reference atom nos. in rigid body for origin 
                             definition, use center of atoms flag (list of 3 int & 1 bool)
\           rbSeq            Orig,Piv,angle,Riding (list): definition of internal rigid body
                             torsion; origin atom (int), pivot atom (int), torsion angle (float),
                             riding atoms (list of int)
\           SelSeq           [int,int] used by SeqSizer to identify objects
\           useCount         Number of times rigid body is used in any structure (int)
RBIds           \            unique Ids generated upon creation of each rigid body (dict)
\           Vector           Ids for each Vector rigid body (list)
\           Residue          Ids for each Residue rigid body (list)
==========  ===============  ====================================================

Space Group Objects 
-------------------

.. _SGData_table:

.. index::
   single: Space Group Data description
   single: Data object descriptions; Space Group Data

Space groups are interpreted by :func:`GSASIIspc.SpcGroup` 
and the information is placed in a SGdata object 
which is a dict with these keys:

.. tabularcolumns:: |l|p{4.5in}|

==========  ====================================================
  key         explanation
==========  ====================================================
SpGrp       space group symbol (str)
Laue        one of the following 14 Laue classes:
            -1, 2/m, mmm, 4/m, 4/mmm, 3R,
            3mR, 3, 3m1, 31m, 6/m, 6/mmm, m3, m3m (str)
SGInv       True if centrosymmetric, False if not (bool)
SGLatt      Lattice centering type. Will be one of
            P, A, B, C, I, F, R (str)
SGUniq      unique axis if monoclinic. Will be
            a, b, or c for monoclinic space groups.
            Will be blank for non-monoclinic. (str)
SGCen       Symmetry cell centering vectors. A (n,3) np.array
            of centers. Will always have at least one row:
            ``np.array([[0, 0, 0]])``
SGOps       symmetry operations as a list of form
            ``[[M1,T1], [M2,T2],...]``
            where :math:`M_n` is a 3x3 np.array
            and :math:`T_n` is a length 3 np.array.
            Atom coordinates are transformed where the
            Asymmetric unit coordinates [X is (x,y,z)]
            are transformed using
            :math:`X^\prime = M_n*X+T_n`
SGSys       symmetry unit cell: type one of
            'triclinic', 'monoclinic', 'orthorhombic',
            'tetragonal', 'rhombohedral', 'trigonal',
            'hexagonal', 'cubic' (str)
SGPolax     Axes for space group polarity. Will be one of
            '', 'x', 'y', 'x y', 'z', 'x z', 'y z',
            'xyz'. In the case where axes are arbitrary 
            '111' is used (P 1, and ?).
==========  ====================================================

Atom Records
------------

.. _Atoms_table:

.. index::
   single: Atoms record description
   single: Data object descriptions; Atoms record


If ``phasedict`` points to the phase information in the data tree, then
atoms are contained in a list of atom records (list) in
``phasedict['Atoms']``. Also needed to read atom information 
are four pointers, ``cx,ct,cs,cia = phasedict['General']['atomPtrs']``,
which define locations in the atom record, as shown below. Items shown are 
always present; additional ones for macromolecular phases are marked 'mm'

.. tabularcolumns:: |l|p{4.5in}|

==============   ====================================================
location         explanation
==============   ====================================================
ct-4              mm - residue number (str)
ct-3              mm - residue name (e.g. ALA) (str)
ct-2              mm - chain label (str)
ct-1              atom label (str)
ct                atom type (str)
ct+1              refinement flags; combination of 'F', 'X', 'U' (str)
cx,cx+1,cx+2      the x,y and z coordinates (3 floats)
cs                site symmetry (str)
cs+1              site multiplicity (int)
cia               ADP flag: Isotropic ('I') or Anisotropic ('A')
cia+1             Uiso (float)
cia+2...cia+7     U11, U22, U33, U12, U13, U23 (6 floats)
atom[cia+8]       unique atom identifier (int)

==============   ====================================================

Drawing Atom Records
--------------------

.. _Drawing_atoms_table:

.. index::
   single: Drawing atoms record description
   single: Data object descriptions; Drawing atoms record


If ``phasedict`` points to the phase information in the data tree, then
drawing atoms are contained in a list of drawing atom records (list) in
``phasedict['Drawing']['Atoms']``. Also needed to read atom information 
are four pointers, ``cx,ct,cs,ci = phasedict['Drawing']['AtomPtrs']``,
which define locations in the atom record, as shown below. Items shown are 
always present; additional ones for macromolecular phases are marked 'mm'

.. tabularcolumns:: |l|p{4.5in}|

==============   ====================================================
location         explanation
==============   ====================================================
ct-4              mm - residue number (str)
ct-3              mm - residue name (e.g. ALA) (str)
ct-2              mm - chain label (str)
ct-1              atom label (str)
ct                atom type (str)
cx,cx+1,cx+2      the x,y and z coordinates (3 floats)
cs-1              Sym Op symbol; sym. op number + unit cell id (e.g. '1,0,-1') (str)
cs                atom drawing style; e.g. 'balls & sticks' (str)
cs+1              atom label style (e.g. 'name') (str)
cs+2              atom color (RBG triplet) (int)
cs+3              ADP flag: Isotropic ('I') or Anisotropic ('A')
cs+4              Uiso (float)
cs+5...cs+11      U11, U22, U33, U12, U13, U23 (6 floats)
ci                unique atom identifier; matches source atom Id in Atom Records (int)
==============   ====================================================

Powder Diffraction Tree Items
-----------------------------

.. _Powder_table:

.. index::
   single: Powder data object description
   single: Data object descriptions; Powder Data

Every powder diffraction histogram is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "PWDR ". The
diffraction data for that information are directly associated with
that tree item and there are a series of children to that item. The
routines :func:`GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ===============  ====================================================
  key                      sub-key        explanation
======================  ===============  ====================================================
Comments                      \           Text strings extracted from the original powder 
                                          data header. These cannot be changed by the user; 
                                          it may be empty.
Limits                       \            A list of two two element lists, as [[Ld,Hd],[L,H]]
                                          where L and Ld are the current and default lowest
                                          two-theta value to be used and 
                                          where H and Hd are the current and default highest
                                          two-theta value to be used.
Reflection Lists              \           A dict with an entry for each phase in the
                                          histogram. The contents of each dict item
                                          is a dict containing reflections, as described in
                                          the :ref:`Powder Reflections <PowderRefl_table>`
                                          description.
Instrument Parameters         \           A list containing two dicts where the possible
                                          keys in each dict are listed below. The value
                                          for each item is a list containing three values:
                                          the initial value, the current value and a
                                          refinement flag which can have a value of
                                          True, False or 0 where 0 indicates a value that
                                          cannot be refined. The first and second
                                          values are floats unless otherwise noted.
                                          Items in the first dict are noted as [1]
\                         Lam             Specifies a wavelength in Angstroms [1]
\                         Lam1            Specifies the primary wavelength in
                                          Angstrom, when an alpha1, alpha2
                                          source is used [1]
\                         Lam2            Specifies the secondary wavelength in
                                          Angstrom, when an alpha1, alpha2
                                          source is used [1]
                          I(L2)/I(L1)     Ratio of Lam2 to Lam1 [1]           
\                         Type            Histogram type (str) [1]: 
                                           * 'PXC' for constant wavelength x-ray
                                           * 'PNC' for constant wavelength neutron
                                           * 'PNT' for time of flight neutron
\                         Zero            Two-theta zero correction in *degrees* [1]
\                         Azimuth         Azimuthal setting angle for data recorded
                                          with differing setting angles [1]
\                         U, V, W         Cagliotti profile coefficients
                                          for Gaussian instrumental broadening, where the
                                          FWHM goes as
                                          :math:`U \\tan^2\\theta + V \\tan\\theta + W` [1]
\                         X, Y            Cauchy (Lorentzian) instrumental broadening
                                          coefficients [1]
\                         SH/L            Variant of the Finger-Cox-Jephcoat asymmetric
                                          peak broadening ratio. Note that this is the
                                          average between S/L and H/L where S is
                                          sample height, H is the slit height and
                                          L is the goniometer diameter. [1]
\                         Polariz.        Polarization coefficient. [1]
wtFactor                      \           A weighting factor to increase or decrease
                                          the leverage of data in the histogram (float).
                                          A value of 1.0 weights the data with their
                                          standard uncertainties and a larger value
                                          increases the weighting of the data (equivalent
                                          to decreasing the uncertainties).
Sample Parameters             \           Specifies a dict with parameters that describe how
                                          the data were collected, as listed
                                          below. Refinable parameters are a list containing
                                          a float and a bool, where the second value
                                          specifies if the value is refined, otherwise
                                          the value is a float unless otherwise noted.
\                         Scale           The histogram scale factor (refinable) 
\                         Absorption      The sample absorption coefficient as
                                          :math:`\\mu r` where r is the radius
                                          (refinable). Only valid for Debye-Scherrer geometry.
\                         SurfaceRoughA   Surface roughness parameter A as defined by
                                          Surotti,J. Appl. Cryst, 5,325-331, 1972.(refinable - 
                                          only valid for Bragg-Brentano geometry)                                         
\                         SurfaceRoughB   Surface roughness parameter B (refinable - 
                                          only valid for Bragg-Brentano geometry)                                          
\                         DisplaceX,      Sample displacement from goniometer center
                          DisplaceY       where Y is along the beam direction and
                                          X is perpendicular. Units are :math:`\\mu m`
                                          (refinable).
\                         Phi, Chi,       Goniometer sample setting angles, in degrees.
                          Omega
\                         Gonio. radius   Radius of the diffractometer in mm
\                         InstrName       A name for the instrument, used in preparing
                                          a CIF (str).
\                         Force,          Variables that describe how the measurement
                          Temperature,    was performed. Not used directly in 
                          Humidity,       any computations. 
                          Pressure,
                          Voltage
\                         ranId           The random-number Id for the histogram
                                          (same value as where top-level key is ranId)
\                         Type            Type of diffraction data, may be 'Debye-Scherrer'
                                          or 'Bragg-Brentano' (str).
\                         Diffuse         not in use?
hId                           \           The number assigned to the histogram when
                                          the project is loaded or edited (can change)
ranId                         \           A random number id for the histogram
                                          that does not change
Background                    \           The background is stored as a list with where
                                          the first item in the list is list and the second
                                          item is a dict. The list contains the background
                                          function and its coefficients; the dict contains
                                          Debye diffuse terms and background peaks.
                                          (TODO: this needs to be expanded.)
Data                          \           The data consist of a list of 6 np.arrays 
                                          containing in order:

                                           0. the x-postions (two-theta in degrees),
                                           1. the intensity values (Yobs),
                                           2. the weights for each Yobs value
                                           3. the computed intensity values (Ycalc)
                                           4. the background values
                                           5. Yobs-Ycalc
======================  ===============  ====================================================

Powder Reflection Data Structure
--------------------------------

.. _PowderRefl_table:

.. index::
   single: Powder reflection object description
   single: Data object descriptions; Powder Reflections
   
For every phase in a histogram, the ``Reflection Lists`` value is a dict
one element of which is `'RefList'`, which is a np.array containing
reflections. The columns in that array are documented below.

==========  ====================================================
  index         explanation
==========  ====================================================
 0,1,2       h,k,l (float)
 3           multiplicity
 4           d-space, Angstrom
 5           pos, two-theta
 6           sig, Gaussian width
 7           gam, Lorenzian width
 8           :math:`F_{obs}^2`
 9           :math:`F_{calc}^2`
 10          reflection phase, in degrees
 11          intensity correction for reflection, this times
             :math:`F_{obs}^2` or :math:`F_{calc}^2` gives Iobs or Icalc 
==========  ====================================================

Single Crystal Tree Items
-------------------------

.. _Xtal_table:

.. index::
   single: Single Crystal data object description
   single: Data object descriptions; Single crystal data

Every single crystal diffraction histogram is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "HKLF ". The
diffraction data for that information are directly associated with
that tree item and there are a series of children to that item. The
routines :func:`GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ===============  ====================================================
  key                      sub-key        explanation
======================  ===============  ====================================================
Data                          \           A dict that contains the 
                                          reflection table,
                                          as described in the
                                          :ref:`Single Crystal Reflections
                                          <XtalRefl_table>`
                                          description.

Instrument Parameters         \           A list containing two dicts where the possible
                                          keys in each dict are listed below. The value
                                          for most items is a list containing two values:
                                          the initial value, the current value.
                                          The first and second
                                          values are floats unless otherwise noted.
\                         Lam             Specifies a wavelength in Angstroms (two floats)
\                         Type            Histogram type (two str values): 
                                           * 'SXC' for constant wavelength x-ray
                                           * 'SNC' for constant wavelength neutron
                                           * 'SNT' for time of flight neutron
\                         InstrName       A name for the instrument, used in preparing
                                          a CIF (str).

wtFactor                      \           A weighting factor to increase or decrease
                                          the leverage of data in the histogram (float).
                                          A value of 1.0 weights the data with their
                                          standard uncertainties and a larger value
                                          increases the weighting of the data (equivalent
                                          to decreasing the uncertainties).

hId                           \           The number assigned to the histogram when
                                          the project is loaded or edited (can change)
ranId                         \           A random number id for the histogram
                                          that does not change
======================  ===============  ====================================================

Single Crystal Reflection Data Structure
----------------------------------------

.. _XtalRefl_table:

.. index::
   single: Single Crystal reflection object description
   single: Data object descriptions; Single Crystal Reflections
   
For every single crystal a histogram, the ``'Data'`` item contains
the structure factors as an np.array in item `'RefList'`.
The columns in that array are documented below.

==========  ====================================================
  index         explanation
==========  ====================================================
 0,1,2       h,k,l (float)
 3           multiplicity
 4           d-space, Angstrom
 5           :math:`F_{obs}^2`
 6           :math:`\sigma(F_{obs}^2)`
 7           :math:`F_{calc}^2`
 8           :math:`F_{obs}^2T`
 9           :math:`F_{calc}^2T`
 10          reflection phase, in degrees
 11          intensity correction for reflection, this times
             :math:`F_{obs}^2` or :math:`F_{calc}^2`
             gives Iobs or Icalc
==========  ====================================================

Image Data Structure
--------------------

.. _Image_table:

.. index::
   image: Image data object description
   image: Image object descriptions
   
Every 2-dimensional image is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "IMG ". The
image data are directly associated with that tree item and there 
are a series of children to that item. The routines :func:`GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ======================  ====================================================
  key                      sub-key              explanation
======================  ======================  ====================================================
Comments                       \                Text strings extracted from the original image data 
                                                header or a metafile. These cannot be changed by  
                                                the user; it may be empty.                                                
Image Controls              azmthOff            (float) The offset to be applied to an azimuthal 
                                                value. Accomodates 
                                                detector orientations other than with the detector 
                                                X-axis
                                                horizontal.
\                           background image    (list:str,float) The name of a tree item ("IMG ...") that is to be subtracted
                                                during image integration multiplied by value. It must have the same size/shape as 
                                                the integrated image. NB: value < 0 for subtraction.
\                           calibrant           (str) The material used for determining the position/orientation
                                                of the image. The data is obtained from :func:`ImageCalibrants` 
                                                and UserCalibrants.py (supplied by user).
\                           calibdmin           (float) The minimum d-spacing used during the last calibration run.
\                           calibskip           (int) The number of expected diffraction lines skipped during the last
                                                calibration run.
\                           center              (list:floats) The [X,Y] point in detector coordinates (mm) where the direct beam
                                                strikes the detector plane as determined by calibration. This point 
                                                does not have to be within the limits of the detector boundaries. 
\                           centerAzm           (bool) If True then the azimuth reported for the integrated slice
                                                of the image is at the center line otherwise it is at the leading edge.
\                           color               (str) The name of the colormap used to display the image. Default = 'Paired'.
\                           cutoff              (float) The minimum value of I/Ib for a point selected in a diffraction ring for 
                                                calibration calculations. See pixLimit for details as how point is found.           
\                           DetDepth            (float) Coefficient for penetration correction to distance; accounts for diffraction
                                                ring offset at higher angles. Optionally determined by calibration.
\                           DetDepthRef         (bool) If True then refine DetDepth during calibration/recalibration calculation.
\                           distance            (float) The distance (mm) from sample to detector plane.
\                           ellipses            (list:lists) Each object in ellipses is a list [center,phi,radii,color] where
                                                center (list) is location (mm) of the ellipse center on the detector plane, phi is the 
                                                rotation of the ellipse minor axis from the x-axis, and radii are the minor & major
                                                radii of the ellipse. If radii[0] is negative then parameters describe a hyperbola. Color
                                                is the selected drawing color (one of 'b', 'g' ,'r') for the ellipse/hyperbola.
\                           edgemin             (float) Not used;  parameter in EdgeFinder code.
\                           fullIntegrate       (bool) If True then integrate over full 360 deg azimuthal range.
\                           GonioAngles         (list:floats) The 'Omega','Chi','Phi' goniometer angles used for this image. 
                                                Required for texture calculations.
\                           invert_x            (bool) If True display the image with the x-axis inverted.
\                           invert_y            (bool) If True display the image with the y-axis inverted.
\                           IOtth               (list:floats) The minimum and maximum 2-theta values to be used for integration.
\                           LRazimuth           (list:floats) The minimum and maximum azimuth values to be used for integration.
\                           Oblique             (list:float,bool) If True apply a detector absorption correction using the value to the
                                                intensities obtained during integration.
\                           outAzimuths         (int) The number of azimuth pie slices.
\                           outChannels         (int) The number of 2-theta steps. 
\                           pixelSize           (list:ints) The X,Y dimensions (microns) of each pixel.
\                           pixLimit            (int) A box in the image with 2*pixLimit+1 edges is searched to find the maximum.
                                                This value (I) along with the minimum (Ib) in the box is reported by :func:`GSASIIimage.ImageLocalMax`
                                                and subject to cutoff in :func:`GSASIIimage.makeRing`. 
                                                Locations are used to construct rings of points for calibration calcualtions.
\                           PolaVal             (list:float,bool) If type='SASD' and if True, apply polarization correction to intensities from 
                                                integration using value.
\                           rings               (list:lists) Each entry is [X,Y,dsp] where X & Y are lists of x,y coordinates around a 
                                                diffraction ring with the same d-spacing (dsp)
\                           ring                (list) The x,y coordinates of the >5 points on an inner ring 
                                                selected by the user,
\                           Range               (list) The minimum & maximum values of the image
\                           rotation            (float) The angle between the x-axis and the vector about which the 
                                                detector is tilted. Constrained to -180 to 180 deg.     
\                           SampleShape         (str) Currently only 'Cylinder'. Sample shape for Debye-Scherrer experiments; used for absorption
                                                calculations.
\                           SampleAbs           (list: float,bool) Value of absorption coefficient for Debye-Scherrer experimnents, flag if True
                                                to cause correction to be applied.
\                           setDefault          (bool) If True the use the image controls values for all new images to be read. (might be removed)
\                           setRings            (bool) If True then display all the selected x,y ring positions (vida supra rings) used in the calibration.            
\                           showLines           (bool) If True then isplay the integration limits to be used.
\                           size                (list:int) The number of pixels on the image x & y axes
\                           type                (str) One of 'PWDR', 'SASD' or 'REFL' for powder, small angle or reflectometry data, respectively.
\                           tilt                (float) The angle the detector normal makes with the incident beam; range -90 to 90.
\                           wavelength          (float) Tha radiation wavelength (Angstroms) as entered by the user (or someday obtained from the image header).
                                                
Masks                       Arcs                (list: lists) Each entry [2-theta,[azimuth[0],azimuth[1]],thickness] describes an arc mask
                                                to be excluded from integration
\                           Frames              (list:lists) Each entry describes the x,y points (3 or more - mm) that describe a frame outside
                                                of which is excluded from recalibration and integration. Only one frame is allowed.
\                           Points              (list:lists) Each entry [x,y,radius] (mm) describes an excluded spot on the image to be excluded
                                                from integration.
\                           Polygons            (list:lists) Each entry is a list of 3+ [x,y] points (mm) that describe a polygon on the image
                                                to be excluded from integration.
\                           Rings               (list: lists) Each entry [2-theta,thickness] describes a ring mask
                                                to be excluded from integration.
\                           Thresholds          (list:[tuple,list]) [(Imin,Imax),[Imin,Imax]] This gives lower and upper limits for points on the image to be included
                                                in integrsation. The tuple is the image intensity limits and the list are those set by the user.    
                                                
Stress/Strain               Sample phi          (float) Sample rotation about vertical axis.
\                           Sample z            (float) Sample translation from the calibration sample position (for Sample phi = 0)
                                                These will be restricted by space group symmetry; result of strain fit refinement.
\                           Type                (str) 'True' or 'Conventional': The strain model used for the calculation.
\                           d-zero              (list:dict) Each item is for a diffraction ring on the image; all items are from the same phase 
                                                and are used to determine the strain tensor.
                                                The dictionary items are:
                                                'Dset': (float) True d-spacing for the diffraction ring; entered by the user.
                                                'Dcalc': (float) Average calculated d-spacing determined from strain coeff.
                                                'Emat': (list: float) The strain tensor elements e11, e12 & e22 (e21=e12, rest are 0)
                                                'Esig': (list: float) Esds for Emat from fitting.
                                                'pixLimit': (int) Search range to find highest point on ring for each data point
                                                'cutoff': (float) I/Ib cutoff for searching.
                                                'ImxyObs': (list: lists) [[X],[Y]] observed points to be used for strain calculations.
                                                'ImtaObs': (list: lists) [[d],[azm]] transformed via detector calibration from ImxyObs.
                                                'ImtaCalc': (list: lists [[d],[azm]] calculated d-spacing & azimuth from fit.
                                                
======================  ======================  ====================================================

Parameter Dictionary
-------------------------

.. _parmDict_table:

.. index::
   single: Parameter dictionary

The parameter dictionary contains all of the variable parameters for the refinement.
The dictionary keys are the name of the parameter (<phase>:<hist>:<name>:<atom>). 
It is prepared in two ways. When loaded from the tree
(in :meth:`GSASII.GSASII.MakeLSParmDict` and
:meth:`GSASIIIO.ExportBaseclass.loadParmDict`), 
the values are lists with two elements: ``[value, refine flag]``

When loaded from the GPX file (in
:func:`GSASIIstrMain.Refine` and :func:`GSASIIstrMain.SeqRefine`), the value in the
dict is the actual parameter value (usually a float, but sometimes a 
letter or string flag value (such as I or A for iso/anisotropic). 


*Classes and routines*
----------------------

'''
import re
import imp
import random as ran
import sys
import GSASIIpath
import GSASIImath as G2mth
import numpy as np

GSASIIpath.SetVersionNumber("$Revision$")

DefaultControls = {
    'deriv type':'analytic Hessian',    #default controls
    'min dM/M':0.0001,'shift factor':1.,'max cyc':3,'F**2':True,
    'minF/sig':0,
    'Author':'no name',
    'FreeVar1':'Sample humidity (%)',
    'FreeVar2':'Sample voltage (V)',
    'FreeVar3':'Applied load (MN)',
    }
'''Values to be used as defaults for the initial contents of the ``Controls``
data tree item.
'''

def MakeUniqueLabel(lbl,labellist):
    '''Make sure that every a label is unique against a list by adding
    digits at the end until it is not found in list.

    :param str lbl: the input label
    :param list labellist: the labels that have already been encountered
    :returns: lbl if not found in labellist or lbl with ``_1-9`` (or
      ``_10-99``, etc.) appended at the end
    '''
    lbl = lbl.strip()
    if not lbl: # deal with a blank label
        lbl = '_1'
    if lbl not in labellist:
        labellist.append(lbl)
        return lbl
    i = 1
    prefix = lbl
    if '_' in lbl:
        prefix = lbl[:lbl.rfind('_')]
        suffix = lbl[lbl.rfind('_')+1:]
        try:
            i = int(suffix)+1
        except: # suffix could not be parsed
            i = 1
            prefix = lbl
    while prefix+'_'+str(i) in labellist:
        i += 1
    else:
        lbl = prefix+'_'+str(i)
        labellist.append(lbl)
    return lbl

PhaseIdLookup = {}
'''dict listing phase name and random Id keyed by sequential phase index as a str;
best to access this using :func:`LookupPhaseName`
'''
PhaseRanIdLookup = {}
'''dict listing phase sequential index keyed by phase random Id;
best to access this using :func:`LookupPhaseId`
'''
HistIdLookup = {}
'''dict listing histogram name and random Id, keyed by sequential histogram index as a str;
best to access this using :func:`LookupHistName`
'''
HistRanIdLookup = {}
'''dict listing histogram sequential index keyed by histogram random Id;
best to access this using :func:`LookupHistId`
'''
AtomIdLookup = {}
'''dict listing for each phase index as a str, the atom label and atom random Id,
keyed by atom sequential index as a str;
best to access this using :func:`LookupAtomLabel`
'''
AtomRanIdLookup = {}
'''dict listing for each phase the atom sequential index keyed by atom random Id;
best to access this using :func:`LookupAtomId`
'''
ShortPhaseNames = {}
'''a dict containing a possibly shortened and when non-unique numbered
version of the phase name. Keyed by the phase sequential index.
'''
ShortHistNames = {}
'''a dict containing a possibly shortened and when non-unique numbered
version of the histogram name. Keyed by the histogram sequential index.
'''

VarDesc = {}
''' This dictionary lists descriptions for GSAS-II variables,
as set in :func:`CompileVarDesc`. See that function for a description
for how keys and values are written.
'''

reVarDesc = {}
''' This dictionary lists descriptions for GSAS-II variables with
the same values as :attr:`VarDesc` except that keys have been compiled as
regular expressions. Initialized in :func:`CompileVarDesc`.
'''

def IndexAllIds(Histograms,Phases):
    '''Scan through the used phases & histograms and create an index
    to the random numbers of phases, histograms and atoms. While doing this,
    confirm that assigned random numbers are unique -- just in case lightning
    strikes twice in the same place.

    Note: this code assumes that the atom random Id (ranId) is the last 
    element each atom record.

    This is called in three places (only): :func:`GSASIIstrIO.GetUsedHistogramsAndPhases`
    (which loads the histograms and phases from a GPX file),
    :meth:`~GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree`
    (which loads the histograms and phases from the data tree.) and
    :meth:`GSASIIconstrGUI.UpdateConstraints`
    (which displays & edits the constraints in a GUI)

    TODO: do we need a lookup for rigid body variables?
    '''
    # process phases and atoms
    PhaseIdLookup.clear()
    PhaseRanIdLookup.clear()    
    AtomIdLookup.clear()
    AtomRanIdLookup.clear()
    ShortPhaseNames.clear()
    for ph in Phases:
        cx,ct,cs,cia = Phases[ph]['General']['AtomPtrs']
        ranId = Phases[ph]['ranId'] 
        while ranId in PhaseRanIdLookup:
            # Found duplicate random Id! note and reassign
            print ("\n\n*** Phase "+str(ph)+" has repeated ranId. Fixing.\n")
            Phases[ph]['ranId'] = ranId = ran.randint(0,sys.maxint)
        pId = str(Phases[ph]['pId'])
        PhaseIdLookup[pId] = (ph,ranId)
        PhaseRanIdLookup[ranId] = pId
        shortname = ph[:10]
        while shortname in ShortPhaseNames.values():
            shortname = ph[:8] + ' ('+ pId + ')'
        ShortPhaseNames[pId] = shortname
        AtomIdLookup[pId] = {}
        AtomRanIdLookup[pId] = {}
        for iatm,at in enumerate(Phases[ph]['Atoms']):
            ranId = at[cia+8]
            while ranId in AtomRanIdLookup[pId]: # check for dups
                print ("\n\n*** Phase "+str(ph)+" atom "+str(iatm)+" has repeated ranId. Fixing.\n")
                at[cia+8] = ranId = ran.randint(0,sys.maxint)
            AtomRanIdLookup[pId][ranId] = str(iatm)
            if Phases[ph]['General']['Type'] == 'macromolecular':
                label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
            else:
                label = at[ct-1]
            AtomIdLookup[pId][str(iatm)] = (label,ranId)
    # process histograms
    HistIdLookup.clear()
    HistRanIdLookup.clear()
    ShortHistNames.clear()
    for hist in Histograms:
        ranId = Histograms[hist]['ranId']
        while ranId in HistRanIdLookup:
            # Found duplicate random Id! note and reassign
            print ("\n\n*** Histogram "+str(hist)+" has repeated ranId. Fixing.\n")
            Histograms[hist]['ranId'] = ranId = ran.randint(0,sys.maxint)
        hId = str(Histograms[hist]['hId'])
        HistIdLookup[hId] = (hist,ranId)
        HistRanIdLookup[ranId] = hId
        shortname = hist[:15]
        while shortname in ShortHistNames.values():
            shortname = hist[:11] + ' ('+ hId + ')'
        ShortHistNames[hId] = shortname

def LookupAtomId(pId,ranId):
    '''Get the atom number from a phase and atom random Id

    :param int/str pId: the sequential number of the phase
    :param int ranId: the random Id assigned to an atom

    :returns: the index number of the atom (str)
    '''
    if not AtomRanIdLookup:
        raise Exception,'Error: LookupAtomId called before IndexAllIds was run'
    if pId is None or pId == '':
        raise KeyError,'Error: phase is invalid (None or blank)'
    pId = str(pId)
    if pId not in AtomRanIdLookup:
        raise KeyError,'Error: LookupAtomId does not have phase '+pId
    if ranId not in AtomRanIdLookup[pId]:
        raise KeyError,'Error: LookupAtomId, ranId '+str(ranId)+' not in AtomRanIdLookup['+pId+']'
    return AtomRanIdLookup[pId][ranId]

def LookupAtomLabel(pId,index):
    '''Get the atom label from a phase and atom index number

    :param int/str pId: the sequential number of the phase
    :param int index: the index of the atom in the list of atoms

    :returns: the label for the atom (str) and the random Id of the atom (int)
    '''
    if not AtomIdLookup:
        raise Exception,'Error: LookupAtomLabel called before IndexAllIds was run'
    if pId is None or pId == '':
        raise KeyError,'Error: phase is invalid (None or blank)'
    pId = str(pId)
    if pId not in AtomIdLookup:
        raise KeyError,'Error: LookupAtomLabel does not have phase '+pId
    if index not in AtomIdLookup[pId]:
        raise KeyError,'Error: LookupAtomLabel, ranId '+str(index)+' not in AtomRanIdLookup['+pId+']'
    return AtomIdLookup[pId][index]

def LookupPhaseId(ranId):
    '''Get the phase number and name from a phase random Id

    :param int ranId: the random Id assigned to a phase
    :returns: the sequential Id (pId) number for the phase (str)
    '''
    if not PhaseRanIdLookup:
        raise Exception,'Error: LookupPhaseId called before IndexAllIds was run'
    if ranId not in PhaseRanIdLookup:
        raise KeyError,'Error: LookupPhaseId does not have ranId '+str(ranId)
    return PhaseRanIdLookup[ranId]

def LookupPhaseName(pId):
    '''Get the phase number and name from a phase Id

    :param int/str pId: the sequential assigned to a phase
    :returns:  (phase,ranId) where phase is the name of the phase (str)
      and ranId is the random # id for the phase (int)
    '''
    if not PhaseIdLookup:
        raise Exception,'Error: LookupPhaseName called before IndexAllIds was run'
    if pId is None or pId == '':
        raise KeyError,'Error: phase is invalid (None or blank)'
    pId = str(pId)
    if pId not in PhaseIdLookup:
        raise KeyError,'Error: LookupPhaseName does not have index '+pId
    return PhaseIdLookup[pId]

def LookupHistId(ranId):
    '''Get the histogram number and name from a histogram random Id

    :param int ranId: the random Id assigned to a histogram
    :returns: the sequential Id (hId) number for the histogram (str)
    '''
    if not HistRanIdLookup:
        raise Exception,'Error: LookupHistId called before IndexAllIds was run'
    if ranId not in HistRanIdLookup:
        raise KeyError,'Error: LookupHistId does not have ranId '+str(ranId)
    return HistRanIdLookup[ranId]

def LookupHistName(hId):
    '''Get the histogram number and name from a histogram Id

    :param int/str hId: the sequential assigned to a histogram
    :returns:  (hist,ranId) where hist is the name of the histogram (str)
      and ranId is the random # id for the histogram (int)
    '''
    if not HistIdLookup:
        raise Exception,'Error: LookupHistName called before IndexAllIds was run'
    if hId is None or hId == '':
        raise KeyError,'Error: histogram is invalid (None or blank)'
    hId = str(hId)
    if hId not in HistIdLookup:
        raise KeyError,'Error: LookupHistName does not have index '+hId
    return HistIdLookup[hId]

def fmtVarDescr(varname):
    '''Return a string with a more complete description for a GSAS-II variable 

    :param str varname: A full G2 variable name with 2 or 3 or 4
       colons (<p>:<h>:name[:<a>] or <p>::RBname:<r>:<t>])
       
    :returns: a string with the description
    '''
    s,l = VarDescr(varname)
    return s+": "+l

def VarDescr(varname):
    '''Return two strings with a more complete description for a GSAS-II variable 

    :param str name: A full G2 variable name with 2 or 3 or 4
       colons (<p>:<h>:name[:<a>] or <p>::RBname:<r>:<t>])
       
    :returns: (loc,meaning) where loc describes what item the variable is mapped
      (phase, histogram, etc.) and meaning describes what the variable does.
    '''
    
    # special handling for parameter names without a colons
    # for now, assume self-defining
    if varname.find(':') == -1:
        return "Global",varname
        
    l = getVarDescr(varname)
    if not l:
        return ("invalid variable name ("+str(varname)+")!"),""
#        return "invalid variable name!",""

    if not l[-1]:
        l[-1] = "(variable needs a definition! Set it in CompileVarDesc)"

    if len(l) == 3:         #SASD variable name!
        s = 'component:'+l[1]
        return s,l[-1]
    s = ""
    if l[0] is not None and l[1] is not None: # HAP: keep short
        if l[2] == "Scale": # fix up ambigous name
            l[5] = "Phase fraction"
        if l[0] == '*':
            lbl = 'all'
        else:
            lbl = ShortPhaseNames.get(l[0],'? #'+str(l[0]))
        if l[1] == '*':
            hlbl = 'all'
        else:
            hlbl = ShortHistNames.get(l[1],'? #'+str(l[1]))
        if hlbl[:4] == 'HKLF':
            hlbl = 'Xtl='+hlbl[5:]
        elif hlbl[:4] == 'PWDR':
            hlbl = 'Pwd='+hlbl[5:]
        else:
            hlbl = 'Hist='+hlbl
        s = "Ph="+str(lbl)+" * "+str(hlbl)
    else:
        if l[2] == "Scale": # fix up ambigous name: must be scale factor, since not HAP
            l[5] = "Scale factor"
        if l[2] == 'Back': # background parameters are "special", alas
            s = 'Hist='+ShortHistNames.get(l[1],'? #'+str(l[1]))
            l[-1] += ' #'+str(l[3])
        elif l[4] is not None: # rigid body parameter
            lbl = ShortPhaseNames.get(l[0],'phase?')
            s = "Res #"+str(l[3])+" body #"+str(l[4])+" in "+str(lbl)
        elif l[3] is not None: # atom parameter, 
            lbl = ShortPhaseNames.get(l[0],'phase?')
            try:
                albl = LookupAtomLabel(l[0],l[3])[0]
            except KeyError:
                albl = 'Atom?'
            s = "Atom "+str(albl)+" in "+str(lbl)
        elif l[0] == '*':
            s = "All phases "
        elif l[0] is not None:
            lbl = ShortPhaseNames.get(l[0],'phase?')
            s = "Phase "+str(lbl)
        elif l[1] == '*':
            s = 'All hists'
        elif l[1] is not None:
            hlbl = ShortHistNames.get(l[1],'? #'+str(l[1]))
            if hlbl[:4] == 'HKLF':
                hlbl = 'Xtl='+hlbl[5:]
            elif hlbl[:4] == 'PWDR':
                hlbl = 'Pwd='+hlbl[5:]
            else:
                hlbl = 'Hist='+hlbl
            s = str(hlbl)
    if not s:
        s = 'Global'
    return s,l[-1]

def getVarDescr(varname):
    '''Return a short description for a GSAS-II variable 

    :param str name: A full G2 variable name with 2 or 3 or 4
       colons (<p>:<h>:name[:<a1>][:<a2>])
      
    :returns: a six element list as [`p`,`h`,`name`,`a1`,`a2`,`description`],
      where `p`, `h`, `a1`, `a2` are str values or `None`, for the phase number,
      the histogram number and the atom number; `name` will always be
      a str; and `description` is str or `None`.
      If the variable name is incorrectly formed (for example, wrong
      number of colons), `None` is returned instead of a list.
    '''
    l = varname.split(':')
    if len(l) == 2:     #SASD parameter name
        return varname,l[0],getDescr(l[1])
    if len(l) == 3:
        l += [None,None]
    elif len(l) == 4:
        l += [None]
    elif len(l) != 5:
        return None
    for i in (0,1,3,4):
        if l[i] == "":
            l[i] = None
    l += [getDescr(l[2])]
    return l
    
def CompileVarDesc():
    '''Set the values in the variable description lookup table (:attr:`VarDesc`)
    into :attr:`reVarDesc`. This is called in :func:`getDescr` so the initialization
    is always done before use.

    Note that keys may contain regular expressions, where '[xyz]'
    matches 'x' 'y' or 'z' (equivalently '[x-z]' describes this as range of values).
    '.*' matches any string. For example::

    'AUiso':'Atomic isotropic displacement parameter',

    will match variable ``'p::AUiso:a'``.
    If parentheses are used in the key, the contents of those parentheses can be
    used in the value, such as::

    'AU([123][123])':'Atomic anisotropic displacement parameter U\\1',

    will match ``AU11``, ``AU23``,.. and `U11`, `U23` etc will be displayed
    in the value when used.
    
    '''
    if reVarDesc: return # already done
    for key,value in {
        # derived or other sequential vars
        '([abc])$' : 'Lattice parameter, \\1, from Ai and Djk', # N.B. '$' prevents match if any characters follow
        u'\u03B1' : u'Lattice parameter, \u03B1, from Ai and Djk',
        u'\u03B2' : u'Lattice parameter, \u03B2, from Ai and Djk',
        u'\u03B3' : u'Lattice parameter, \u03B3, from Ai and Djk',
        # ambiguous, alas:
        'Scale' : 'Phase or Histogram scale factor',
        # Phase vars (p::<var>)
        'A([0-5])' : 'Reciprocal metric tensor component \\1',
        'Vol' : 'Unit cell volume',
        # Atom vars (p::<var>:a)
        'dA([xyz])$' : 'change to atomic coordinate, \\1',
        'A([xyz])$' : '\\1 fractional atomic coordinate',
        'AUiso':'Atomic isotropic displacement parameter',
        'AU([123][123])':'Atomic anisotropic displacement parameter U\\1',
        'Afrac': 'Atomic occupancy parameter',
        # Hist & Phase (HAP) vars (p:h:<var>)
        'Back': 'Background term',
        'BkPkint;(.*)':'Background peak #\\1 intensity',
        'BkPkpos;(.*)':'Background peak #\\1 position',
        'BkPksig;(.*)':'Background peak #\\1 Gaussian width',
        'BkPkgam;(.*)':'Background peak #\\1 Cauchy width',
        'Bab([AU])': 'Babinet solvent scattering coef. \\1',
        'D([123][123])' : 'Anisotropic strain coef. \\1',
        'Extinction' : 'Extinction coef.',
        'MD' : 'March-Dollase coef.',
        'Mustrain;.*' : 'Microstrain coef.',
        'Size;.*' : 'Crystallite size value',
        'eA$' : 'Cubic mustrain value',
        'Ep$' : 'Primary extinction',
        'Es$' : 'Secondary type II extinction',
        'Eg$' : 'Secondary type I extinction',
        #Histogram vars (:h:<var>)
        'Absorption' : 'Absorption coef.',
        'Displace([XY])' : 'Debye-Scherrer sample displacement \\1',
        'Lam' : 'Wavelength',
        'Polariz\.' : 'Polarization correction',
        'SH/L' : 'FCJ peak asymmetry correction',
        '([UVW])$' : 'Gaussian instrument broadening \\1',
        '([XY])$' : 'Cauchy instrument broadening \\1',
        'Zero' : 'Debye-Scherrer zero correction',
        'nDebye' : 'Debye model background corr. terms',
        'nPeaks' : 'Fixed peak background corr. terms',
        'RBV.*' : 'Vector rigid body parameter',
        'RBR.*' : 'Residue rigid body parameter',
        'RBRO([aijk])' : 'Residue rigid body orientation parameter',
        'RBRP([xyz])' : 'Residue rigid body position parameter',
        'RBRTr;.*' : 'Residue rigid body torsion parameter',
        'RBR([TLS])([123AB][123AB])' : 'Residue rigid body group disp. param.',
        'constr([0-9]*)' : 'Parameter from constraint',
        # supersymmetry parameters  p::<var>:a:o 'Flen','Fcent'?
        'mV([0-2])$' : 'Modulation vector component \\1',
        'Fsin$'  :   'Sin site fraction modulation',
        'Fcos$'  :   'Cos site fraction modulation',
        '([XYZ])sin'  : 'Sin position wave for \\1',
        '([XYZ])cos'  : 'Cos position wave for \\1',
        'U([123][123])sin$' :  'Sin thermal wave for U\\1',
        'U([123][123])cos$' :  'Cos thermal wave for U\\1',
        'M([XYZ])sin$' :  'Sin mag. moment wave for \\1',
        'M([XYZ])cos$' :  'Cos mag. moment wave for \\1',
        # SASD vars (l:<var>;l = component)
        'Aspect ratio' : 'Particle aspect ratio',
        'Length' : 'Cylinder length',
        'Diameter' : 'Cylinder/disk diameter',
        'Thickness' : 'Disk thickness',
        'Dist' : 'Interparticle distance',
        'VolFr' : 'Dense scatterer volume fraction',
        'epis' : 'Sticky sphere epsilon',
        'Sticky' : 'Stickyness',
        'Depth' : 'Well depth',
        'Width' : 'Well width',
        'Volume' : 'Particle volume',
        'Radius' : 'Sphere/cylinder/disk radius',
        'Mean' : 'Particle mean radius',
        'StdDev' : 'Standard deviation in Mean',
        'G$': 'Guinier prefactor',
        'Rg$': 'Guinier radius of gyration',
        'B$': 'Porod prefactor',
        'P$': 'Porod power',
        'Cutoff': 'Porod cutoff',
        'PkInt': 'Bragg peak intensity',
        'PkPos': 'Bragg peak position',
        'PkSig': 'Bragg peak sigma',
        'PkGam': 'Bragg peak gamma',
        'e([12][12])' : 'strain tensor e\1',   # strain vars e11, e22, e12
        'Dcalc': 'Calc. d-spacing',
        'Back$': 'background parameter',
        'pos$': 'peak position',
        'int$': 'peak intensity',
        }.items():
        VarDesc[key] = value
        reVarDesc[re.compile(key)] = value

def getDescr(name):
    '''Return a short description for a GSAS-II variable 

    :param str name: The descriptive part of the variable name without colons (:)
      
    :returns: a short description or None if not found
    '''

    CompileVarDesc() # compile the regular expressions, if needed
    for key in reVarDesc:
        m = key.match(name)
        if m:
            reVarDesc[key]
            return m.expand(reVarDesc[key])
    return None

def GenWildCard(varlist):
    '''Generate wildcard versions of G2 variables. These introduce '*'
    for a phase, histogram or atom number (but only for one of these
    fields) but only when there is more than one matching variable in the
    input variable list. So if the input is this::
    
      varlist = ['0::AUiso:0', '0::AUiso:1', '1::AUiso:0']

    then the output will be this::
    
       wildList = ['*::AUiso:0', '0::AUiso:*']

    :param list varlist: an input list of GSAS-II variable names
      (such as 0::AUiso:0)

    :returns: wildList, the generated list of wild card variable names.
    '''
    wild = []
    for i in (0,1,3):
        currentL = varlist[:]
        while currentL:
            item1 = currentL.pop(0)
            i1splt = item1.split(':')
            if i >= len(i1splt): continue
            if i1splt[i]:
                nextL = []
                i1splt[i] = '[0-9]+'
                rexp = re.compile(':'.join(i1splt))
                matchlist = [item1]
                for nxtitem in currentL:
                    if rexp.match(nxtitem):
                        matchlist += [nxtitem]
                    else:
                        nextL.append(nxtitem)
                if len(matchlist) > 1:
                    i1splt[i] = '*'
                    wild.append(':'.join(i1splt))
                currentL = nextL
    return wild

def LookupWildCard(varname,varlist):
    '''returns a list of variable names from list varname
    that match wildcard name in varname
    
    :param str varname: a G2 variable name containing a wildcard
      (such as \*::var)
    :param list varlist: the list of all variable names used in
      the current project
    :returns: a list of matching GSAS-II variables (may be empty)  
    '''
    rexp = re.compile(varname.replace('*','[0-9]+'))
    return sorted([var for var in varlist if rexp.match(var)])


def _lookup(dic,key):
    '''Lookup a key in a dictionary, where None returns an empty string
    but an unmatched key returns a question mark. Used in :class:`G2VarObj`
    '''
    if key is None:
        return ""
    elif key == "*":
        return "*"
    else:
        return dic.get(key,'?')

class G2VarObj(object):
    '''Defines a GSAS-II variable either using the phase/atom/histogram
    unique Id numbers or using a character string that specifies
    variables by phase/atom/histogram number (which can change).
    Note that :func:`LoadID` should be used to (re)load the current Ids
    before creating or later using the G2VarObj object.

    This can store rigid body variables, but does not translate the residue # and
    body # to/from random Ids

    A :class:`G2VarObj` object can be created with a single parameter:
    
    :param str/tuple varname: a single value can be used to create a :class:`G2VarObj`
      object. If a string, it must be of form "p:h:var" or "p:h:var:a", where

     * p is the phase number (which may be left blank or may be '*' to indicate all phases); 
     * h is the histogram number (which may be left blank or may be '*' to indicate all histograms); 
     * a is the atom number (which may be left blank in which case the third colon is omitted).
       The atom number can be specified as '*' if a phase number is specified (not as '*').
       For rigid body variables, specify a will be a string of form "residue:body#"

      Alternately a single tuple of form (Phase,Histogram,VarName,AtomID) can be used, where
      Phase, Histogram, and AtomID are None or are ranId values (or one can be '*')
      and VarName is a string. Note that if Phase is '*' then the AtomID is an atom number.
      For a rigid body variables, AtomID is a string of form "residue:body#".

    If four positional arguments are supplied, they are:

    :param str/int phasenum: The number for the phase (or None or '*')
    :param str/int histnum: The number for the histogram (or None or '*')
    :param str varname: a single value can be used to create a :class:`G2VarObj`
    :param str/int atomnum: The number for the atom (or None or '*')
    
    '''
    IDdict = {}
    IDdict['phases'] = {}
    IDdict['hists'] = {}
    IDdict['atoms'] = {}
    def __init__(self,*args):
        self.phase = None
        self.histogram = None
        self.name = ''
        self.atom = None
        if len(args) == 1 and (type(args[0]) is list or type(args[0]) is tuple) and len(args[0]) == 4:
            # single arg with 4 values
            self.phase,self.histogram,self.name,self.atom = args[0]
        elif len(args) == 1 and ':' in args[0]:
            #parse a string 
            lst = args[0].split(':')
            if lst[0] == '*':
                self.phase = '*'
                if len(lst) > 3:
                    self.atom = lst[3]
                self.histogram = HistIdLookup.get(lst[1],[None,None])[1]
            elif lst[1] == '*':            
                self.histogram = '*'
                self.phase = PhaseIdLookup.get(lst[0],[None,None])[1]
            else:
                self.histogram = HistIdLookup.get(lst[1],[None,None])[1]
                self.phase = PhaseIdLookup.get(lst[0],[None,None])[1]
                if len(lst) == 4:
                    if lst[3] == '*':
                        self.atom = '*'
                    else:
                        self.atom = AtomIdLookup[lst[0]].get(lst[3],[None,None])[1]
                elif len(lst) == 5:
                    self.atom = lst[3]+":"+lst[4]
                elif len(lst) == 3:
                    pass
                else:
                    raise Exception,"Too many colons in var name "+str(args[0])
            self.name = lst[2]
        elif len(args) == 4:
            if args[0] == '*':
                self.phase = '*'
                self.atom = args[3]
            else:
                self.phase = PhaseIdLookup.get(str(args[0]),[None,None])[1]
                if args[3] == '*':
                    self.atom = '*'
                elif args[0] is not None:
                    self.atom = AtomIdLookup[args[0]].get(str(args[3]),[None,None])[1]
            if args[1] == '*':
                self.histogram = '*'
            else:
                self.histogram = HistIdLookup.get(str(args[1]),[None,None])[1]
            self.name = args[2]
        else:
            raise Exception,"Incorrectly called GSAS-II parameter name"

        #print "DEBUG: created ",self.phase,self.histogram,self.name,self.atom

    def __str__(self):
        return self.varname()

    def varname(self):
        '''Formats the GSAS-II variable name as a "traditional" GSAS-II variable 
        string (p:h:<var>:a) or (p:h:<var>)

        :returns: the variable name as a str
        '''
        a = ""
        if self.phase == "*":
            ph = "*"
            if self.atom:
                a = ":" + str(self.atom)
        else:
            ph = _lookup(PhaseRanIdLookup,self.phase)
            if self.atom == '*':
                a = ':*'
            elif self.atom:
                if ":" in str(self.atom):
                    a = ":" + str(self.atom)
                elif ph in AtomRanIdLookup:
                    a = ":" + AtomRanIdLookup[ph].get(self.atom,'?')
                else:
                    a = ":?"
        if self.histogram == "*":
            hist = "*"
        else:
            hist = _lookup(HistRanIdLookup,self.histogram)
        s = (ph + ":" + hist + ":" + str(self.name)) + a
        return s
    
    def __repr__(self):
        '''Return the detailed contents of the object
        '''
        s = "<"
        if self.phase == '*':
            s += "Phases: all; "
            if self.atom is not None:
                if ":" in str(self.atom):
                    s += "Rigid body" + str(self.atom) + "; "
                else:
                    s += "Atom #" + str(self.atom) + "; "
        elif self.phase is not None:
            ph =  _lookup(PhaseRanIdLookup,self.phase)
            s += "Phase: rId=" + str(self.phase) + " (#"+ ph + "); "
            if self.atom == '*':
                s += "Atoms: all; "
            elif ":" in self(self.atom):
                s += "Rigid body" + str(self.atom) + "; "
            elif self.atom is not None:
                s += "Atom rId=" + str(self.atom)
                if ph in AtomRanIdLookup:
                    s += " (#" + AtomRanIdLookup[ph].get(self.atom,'?') + "); "
                else:
                    s += " (#? -- not found!); "
        if self.histogram == '*':
            s += "Histograms: all; "
        elif self.histogram is not None:
            hist = _lookup(HistRanIdLookup,self.histogram)
            s += "Histogram: rId=" + str(self.histogram) + " (#"+ hist + "); "
        s += 'Variable name="' + str(self.name) + '">'
        return s+" ("+self.varname()+")"

    def __eq__(self, other):
        if type(other) is type(self):
            return (self.phase == other.phase and
                    self.histogram == other.histogram and
                    self.name == other.name and
                    self.atom == other.atom)
        return False

    def _show(self):
        'For testing, shows the current lookup table'
        print 'phases', self.IDdict['phases']
        print 'hists', self.IDdict['hists']
        print 'atomDict', self.IDdict['atoms']

#==========================================================================
# shortcut routines
exp = np.exp
sind = sin = s = lambda x: np.sin(x*np.pi/180.)
cosd = cos = c = lambda x: np.cos(x*np.pi/180.)
tand = tan = t = lambda x: np.tan(x*np.pi/180.)
sqrt = sq = lambda x: np.sqrt(x)
pi = lambda: np.pi
class ExpressionObj(object):
    '''Defines an object with a user-defined expression, to be used for
    secondary fits or restraints. Object is created null, but is changed
    using :meth:`LoadExpression`. This contains only the minimum
    information that needs to be stored to save and load the expression
    and how it is mapped to GSAS-II variables. 
    '''
    def __init__(self):
        self.expression = ''
        'The expression as a text string'
        self.assgnVars = {}
        '''A dict where keys are label names in the expression mapping to a GSAS-II
        variable. The value a G2 variable name.
        Note that the G2 variable name may contain a wild-card and correspond to
        multiple values.
        '''
        self.freeVars = {}
        '''A dict where keys are label names in the expression mapping to a free
        parameter. The value is a list with:

         * a name assigned to the parameter
         * a value for to the parameter and
         * a flag to determine if the variable is refined.
        ''' 
        self.depVar = None

        self.lastError = ('','')
        '''Shows last encountered error in processing expression
        (list of 1-3 str values)'''

    def LoadExpression(self,expr,exprVarLst,varSelect,varName,varValue,varRefflag):
        '''Load the expression and associated settings into the object. Raises
        an exception if the expression is not parsed, if not all functions
        are defined or if not all needed parameter labels in the expression
        are defined.

        This will not test if the variable referenced in these definitions
        are actually in the parameter dictionary. This is checked when the
        computation for the expression is done in :meth:`SetupCalc`.
        
        :param str expr: the expression
        :param list exprVarLst: parameter labels found in the expression
        :param dict varSelect: this will be 0 for Free parameters
          and non-zero for expression labels linked to G2 variables.
        :param dict varName: Defines a name (str) associated with each free parameter
        :param dict varValue: Defines a value (float) associated with each free parameter
        :param dict varRefflag: Defines a refinement flag (bool)
          associated with each free parameter
        '''
        self.expression = expr
        self.compiledExpr = None
        self.freeVars = {}
        self.assgnVars = {}
        for v in exprVarLst:
            if varSelect[v] == 0:
                self.freeVars[v] = [
                    varName.get(v),
                    varValue.get(v),
                    varRefflag.get(v),
                    ]
            else:
                self.assgnVars[v] = varName[v]
        self.CheckVars()

    def EditExpression(self,exprVarLst,varSelect,varName,varValue,varRefflag):
        '''Load the expression and associated settings from the object into
        arrays used for editing.

        :param list exprVarLst: parameter labels found in the expression
        :param dict varSelect: this will be 0 for Free parameters
          and non-zero for expression labels linked to G2 variables.
        :param dict varName: Defines a name (str) associated with each free parameter
        :param dict varValue: Defines a value (float) associated with each free parameter
        :param dict varRefflag: Defines a refinement flag (bool)
          associated with each free parameter

        :returns: the expression as a str
        '''
        for v in self.freeVars:
            varSelect[v] = 0
            varName[v] = self.freeVars[v][0]
            varValue[v] = self.freeVars[v][1]
            varRefflag[v] = self.freeVars[v][2]
        for v in self.assgnVars:
            varSelect[v] = 1
            varName[v] = self.assgnVars[v]
        return self.expression

    def GetVaried(self):
        'Returns the names of the free parameters that will be refined'
        return ["::"+self.freeVars[v][0] for v in self.freeVars if self.freeVars[v][2]]

    def GetVariedVarVal(self):
        'Returns the names and values of the free parameters that will be refined'
        return [("::"+self.freeVars[v][0],self.freeVars[v][1]) for v in self.freeVars if self.freeVars[v][2]]

    def UpdateVariedVars(self,varyList,values):
        'Updates values for the free parameters (after a refinement); only updates refined vars'
        for v in self.freeVars:
            if not self.freeVars[v][2]: continue
            if "::"+self.freeVars[v][0] not in varyList: continue
            indx = varyList.index("::"+self.freeVars[v][0])
            self.freeVars[v][1] = values[indx]

    def GetIndependentVars(self):
        'Returns the names of the required independent parameters used in expression'
        return [self.assgnVars[v] for v in self.assgnVars]

    def CheckVars(self):
        '''Check that the expression can be parsed, all functions are
        defined and that input loaded into the object is internally
        consistent. If not an Exception is raised.

        :returns: a dict with references to packages needed to
          find functions referenced in the expression.
        '''
        ret = self.ParseExpression(self.expression)
        if not ret:
            raise Exception("Expression parse error")
        exprLblList,fxnpkgdict = ret
        # check each var used in expression is defined
        defined = self.assgnVars.keys() + self.freeVars.keys()
        notfound = []
        for var in exprLblList:
            if var not in defined:
                notfound.append(var)
        if notfound:
            msg = 'Not all variables defined'
            msg1 = 'The following variables were not defined: '
            msg2 = ''
            for var in notfound:
                if msg: msg += ', '
                msg += var
            self.lastError = (msg1,'  '+msg2)
            raise Exception(msg)
        return fxnpkgdict

    def ParseExpression(self,expr):
        '''Parse an expression and return a dict of called functions and
        the variables used in the expression. Returns None in case an error
        is encountered. If packages are referenced in functions, they are loaded 
        and the functions are looked up into the modules global
        workspace. 
        
        Note that no changes are made to the object other than
        saving an error message, so that this can be used for testing prior
        to the save. 

        :returns: a list of used variables
        '''
        self.lastError = ('','')
        import ast
        def FindFunction(f):
            '''Find the object corresponding to function f
            :param str f: a function name such as 'numpy.exp'
            :returns: (pkgdict,pkgobj) where pkgdict contains a dict
              that defines the package location(s) and where pkgobj
              defines the object associated with the function.
              If the function is not found, pkgobj is None.
            '''
            df = f.split('.')
            pkgdict = {}
            # no listed package, try in current namespace
            if len(df) == 1: 
                try:
                    fxnobj = eval(f)
                    return pkgdict,fxnobj
                except (AttributeError, NameError):
                    return None,None
            else:
                try:
                    fxnobj = eval(f)
                    pkgdict[df[0]] = eval(df[0])
                    return pkgdict,fxnobj
                except (AttributeError, NameError):
                    pass
            # includes a package, lets try to load the packages
            pkgname = ''
            path = sys.path
            for pkg in f.split('.')[:-1]: # if needed, descend down the tree
                if pkgname:
                    pkgname += '.' + pkg
                else:
                    pkgname = pkg
                fp = None
                try:
                    fp, fppath,desc = imp.find_module(pkg,path)
                    pkgobj = imp.load_module(pkg,fp,fppath,desc)
                    pkgdict[pkgname] = pkgobj
                    path = [fppath]
                except Exception as msg:
                    print('load of '+pkgname+' failed with error='+str(msg))
                    return {},None
                finally:
                    if fp: fp.close()
                try:
                    #print 'before',pkgdict.keys()
                    fxnobj = eval(f,globals(),pkgdict)
                    #print 'after 1',pkgdict.keys()
                    #fxnobj = eval(f,pkgdict)
                    #print 'after 2',pkgdict.keys()
                    return pkgdict,fxnobj
                except:
                    continue
            return None # not found
        def ASTtransverse(node,fxn=False):
            '''Transverse a AST-parsed expresson, compiling a list of variables
            referenced in the expression. This routine is used recursively.

            :returns: varlist,fxnlist where
              varlist is a list of referenced variable names and
              fxnlist is a list of used functions
            '''
            varlist = []
            fxnlist = []
            if isinstance(node, list):
                for b in node:
                    v,f = ASTtransverse(b,fxn)
                    varlist += v
                    fxnlist += f
            elif isinstance(node, ast.AST):
                for a, b in ast.iter_fields(node):
                    if isinstance(b, ast.AST):
                        if a == 'func': 
                            fxnlist += ['.'.join(ASTtransverse(b,True)[0])]
                            continue
                        v,f = ASTtransverse(b,fxn)
                        varlist += v
                        fxnlist += f
                    elif isinstance(b, list):
                        v,f = ASTtransverse(b,fxn)
                        varlist += v
                        fxnlist += f
                    elif node.__class__.__name__ == "Name":
                        varlist += [b]
                    elif fxn and node.__class__.__name__ == "Attribute":
                        varlist += [b]
            return varlist,fxnlist
        try:
            exprast = ast.parse(expr)
        except SyntaxError as err:
            s = ''
            import traceback
            for i in traceback.format_exc().splitlines()[-3:-1]:
                if s: s += "\n"
                s += str(i)
            self.lastError = ("Error parsing expression:",s)
            return
        # find the variables & functions
        v,f = ASTtransverse(exprast)
        varlist = sorted(list(set(v)))
        fxnlist = list(set(f))
        pkgdict = {}
        # check the functions are defined
        for fxn in fxnlist:
            fxndict,fxnobj = FindFunction(fxn)
            if not fxnobj:
                self.lastError = ("Error: Invalid function",fxn,
                                  "is not defined")
                return
            if not hasattr(fxnobj,'__call__'):
                self.lastError = ("Error: Not a function.",fxn,
                                  "cannot be called as a function")
                return
            pkgdict.update(fxndict)
        return varlist,pkgdict

    def GetDepVar(self):
        'return the dependent variable, or None'
        return self.depVar

    def SetDepVar(self,var):
        'Set the dependent variable, if used'
        self.depVar = var
#==========================================================================
class ExpressionCalcObj(object):
    '''An object used to evaluate an expression from a :class:`ExpressionObj`
    object.
    
    :param ExpressionObj exprObj: a :class:`~ExpressionObj` expression object with
      an expression string and mappings for the parameter labels in that object.
    '''
    def __init__(self,exprObj):
        self.eObj = exprObj
        'The expression and mappings; a :class:`ExpressionObj` object'
        self.compiledExpr = None
        'The expression as compiled byte-code'
        self.exprDict = {}
        '''dict that defines values for labels used in expression and packages
        referenced by functions
        '''
        self.lblLookup = {}
        '''Lookup table that specifies the expression label name that is
        tied to a particular GSAS-II parameters in the parmDict.
        '''
        self.fxnpkgdict = {}
        '''a dict with references to packages needed to
        find functions referenced in the expression.
        '''
        self.varLookup = {}
        '''Lookup table that specifies the GSAS-II variable(s)
        indexed by the expression label name. (Used for only for diagnostics
        not evaluation of expression.)
        '''
        # Patch: for old-style expressions with a (now removed step size)
        for v in self.eObj.assgnVars:
            if not isinstance(self.eObj.assgnVars[v], basestring):
                self.eObj.assgnVars[v] = self.eObj.assgnVars[v][0]

    def SetupCalc(self,parmDict):
        '''Do all preparations to use the expression for computation.
        Adds the free parameter values to the parameter dict (parmDict).
        '''
        self.fxnpkgdict = self.eObj.CheckVars()
        # all is OK, compile the expression
        self.compiledExpr = compile(self.eObj.expression,'','eval')

        # look at first value in parmDict to determine its type
        parmsInList = True
        for key in parmDict:
            val = parmDict[key]
            if isinstance(val, basestring):
                parmsInList = False
                break
            try: # check if values are in lists
                val = parmDict[key][0]
            except (TypeError,IndexError):
                parmsInList = False
            break
            
        # set up the dicts needed to speed computations
        self.exprDict = {}
        self.lblLookup = {}
        self.varLookup = {}
        for v in self.eObj.freeVars:
            varname = self.eObj.freeVars[v][0]
            varname = "::" + varname.lstrip(':').replace(' ','_').replace(':',';')
            self.lblLookup[varname] = v
            self.varLookup[v] = varname
            if parmsInList:
                parmDict[varname] = [self.eObj.freeVars[v][1],self.eObj.freeVars[v][2]]
            else:
                parmDict[varname] = self.eObj.freeVars[v][1]
            self.exprDict[v] = self.eObj.freeVars[v][1]
        for v in self.eObj.assgnVars:
            varname = self.eObj.assgnVars[v]
            if '*' in varname:
                varlist = LookupWildCard(varname,parmDict.keys())
                if len(varlist) == 0:
                    raise Exception,"No variables match "+str(v)
                for var in varlist:
                    self.lblLookup[var] = v
                if parmsInList:
                    self.exprDict[v] = np.array([parmDict[var][0] for var in varlist])
                else:
                    self.exprDict[v] = np.array([parmDict[var] for var in varlist])
                self.varLookup[v] = [var for var in varlist]
            elif varname in parmDict:
                self.lblLookup[varname] = v
                self.varLookup[v] = varname
                if parmsInList:
                    self.exprDict[v] = parmDict[varname][0]
                else:
                    self.exprDict[v] = parmDict[varname]
            else:
                raise Exception,"No value for variable "+str(v)
        self.exprDict.update(self.fxnpkgdict)

    def UpdateVars(self,varList,valList):
        '''Update the dict for the expression with a set of values
        :param list varList: a list of variable names
        :param list valList: a list of corresponding values
        '''
        for var,val in zip(varList,valList):
            self.exprDict[self.lblLookup.get(var,'undefined: '+var)] = val

    def UpdateDict(self,parmDict):
        '''Update the dict for the expression with values in a dict
        :param list parmDict: a dict of values some of which may be in use here
        '''
        for var in parmDict:
            if var in self.lblLookup:
                self.exprDict[self.lblLookup[var]] = parmDict[var]
            
    def EvalExpression(self):
        '''Evaluate an expression. Note that the expression
        and mapping are taken from the :class:`ExpressionObj` expression object
        and the parameter values were specified in :meth:`SetupCalc`.
        :returns: a single value for the expression. If parameter
        values are arrays (for example, from wild-carded variable names),
        the sum of the resulting expression is returned.

        For example, if the expression is ``'A*B'``,
        where A is 2.0 and B maps to ``'1::Afrac:*'``, which evaluates to::

        [0.5, 1, 0.5]

        then the result will be ``4.0``.
        '''
        if self.compiledExpr is None:
            raise Exception,"EvalExpression called before SetupCalc"
        val = eval(self.compiledExpr,globals(),self.exprDict)
        if not np.isscalar(val):
            val = np.sum(val)
        return val


if __name__ == "__main__":
    # test equation evaluation
    def showEQ(calcobj):
        print 50*'='
        print calcobj.eObj.expression,'=',calcobj.EvalExpression()
        for v in sorted(calcobj.varLookup):
            print "  ",v,'=',calcobj.exprDict[v],'=',calcobj.varLookup[v]
        # print '  Derivatives'
        # for v in calcobj.derivStep.keys():
        #     print '    d(Expr)/d('+v+') =',calcobj.EvalDeriv(v)

    obj = ExpressionObj()

    obj.expression = "A*np.exp(B)"
    obj.assgnVars =  {'B': '0::Afrac:1'}
    obj.freeVars =  {'A': [u'A', 0.5, True]}
    #obj.CheckVars()
    parmDict2 = {'0::Afrac:0':[0.0,True], '0::Afrac:1': [1.0,False]}
    calcobj = ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)

    obj.expression = "A*np.exp(B)"
    obj.assgnVars =  {'B': '0::Afrac:*'}
    obj.freeVars =  {'A': [u'Free Prm A', 0.5, True]}
    #obj.CheckVars()
    parmDict1 = {'0::Afrac:0':1.0, '0::Afrac:1': 1.0}
    calcobj = ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict1)
    showEQ(calcobj)

    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)
