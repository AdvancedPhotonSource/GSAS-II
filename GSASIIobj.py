# TODO: change this to assemble the look-up tables of atoms, phases and hists from the tree
# and then save/unsave those values in __init__ & __str__, etc. 

'''
*GSASIIobj: Data objects*
=========================

This module defines and/or documents the data structures used in GSAS-II.


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

Each constraint is defined as a list using a series of terms of form

::

[[<mult1>, <var1>], [<mult2>, <var2>],..., <fixed val>, <vary flag>, <cons type>]

Where the variable pair list item containing two values [<mult>, <var>], 

  * <mult> is a multiplier for the constraint (float)
  * <var> is the name of the variable (str) (or to be implemented a :class:`VarName` object.)

Note that the last three items in the list play a special role:

 * <fixed val> is the fixed value for a constraint equation or is None
 * <vary flag> is True, False or None and is intended to use to indicate if new variables
   should be refined.
 * <cons type> is one of four letters, 'e', 'c', 'h', 'f' that determines the type of constraint.

    * 'e' defines a set of equivalent variables. Only the first variable is refined (if the
      appropriate refine flag is set) and and all other equivalent variables in the list
      are generated from that variable. The vary flag for those variables is ignored.
    * 'c' defines a constraint equation of form, :math:`m_1 \\times var_1 + m_2 \\times var_2 + ... = c`
    * 'h' defines a variable to hold (not vary). Any variable on this list is not varied, even if its refinement
      flag is set. This is of particular value when needing to hold one or more variables in a set such as
      the reciprocal metric tensor or anisotropic displacement parameter. 
    * 'f' defines a relationship to define a new variable according to relationship 
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

Space Group Objects
-------------------

.. _SGData_table:

.. index::
   single: SGData description
   single: Data object descriptions; SGData

Space groups are interpreted by :func:`GSASIIspc.SpcGroup` 
and the information is placed in a SGdata object,
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
are four pointers, ``cx,ct,cs,cia = phasedict['General']['AtomPtrs']``,
which define locations in the atom record, as shown below. Items shown are 
always present; additional ones for macromolecular phases are marked 'mm'

.. tabularcolumns:: |l|p{4.5in}|

==============   ====================================================
location         explanation
==============   ====================================================
cx,cx+1,cx+2      the x,y and z coordinates
cx+3              fractional occupancy (also cs-1)
ct-4              (mm) residue number (str)
ct-3              (mm) residue name (e.g. ALA) (str)
ct-2              (mm) chain label (str)
ct-1              atom label
ct                atom type
ct+1              refinement flags
cs                site symmetry string
cs+1              site multiplicity
cia               ADP flag: Isotropic ('I') or Anisotropic ('A')
cia+1             Uiso
cia+2...cia+6     U11, U22, U33, U12, U13, U23
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
routine :func:`~GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ===============  ====================================================
  key                      sub-key        explanation
======================  ===============  ====================================================
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
                                          (refinable).
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

                                           1. the x-postions (two-theta in degrees),
                                           2. the intensity values (Yobs),
                                           3. the weights for each Yobs value
                                           4. the computed intensity values (Ycalc)
                                           5. the background values
                                           6. Yobs-Ycalc
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
routine :func:`~GSASII.GSASII.GetUsedHistogramsAndPhasesfromTree` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ===============  ====================================================
  key                      sub-key        explanation
======================  ===============  ====================================================
Data
                                          A dict that contains the 
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
======================  ===============  ====================================================

Single Crystal Reflection Data Structure
----------------------------------------

.. _XtalRefl_table:

.. index::
   single: Single Crystal reflection object description
   single: Data object descriptions; Single Crystal Reflections
   
For every simgle crystal a histogram, the ``'Data'`` item contains
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


*Classes and routines*
----------------------

'''
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

def LoadHistogramIDs(histList,idList):
    '''Save the Id values for a series of histograms'''
    VarName.IDdict['hists'] = {}
    for h,i in zip(histList,idList):
        VarName.IDdict['hists'][i] = h

def LoadPhaseIDs(self):
    pass

class VarName(object):
    '''Defines a GSAS-II variable either using the phase/atom/histogram
    unique Id numbers or using a character string that specifies
    variables by phase/atom/histogram number (which can change).
    Note that :func:`LoadID` should be used to (re)load the current Ids
    before creating or later using the VarName object.

    A :class:`VarName` object can be created with a single parameter:
    
    :param str varname: a single value can be used to create a :class:`VarName`
      object. The string must be of form "p:h:var" or "p:h:var:a", where

     * p is the phase number (which may be left blank); 
     * h is the histogram number (which may be left blank); 
     * a is the atom number (which may be left blank in which case the third colon is omitted).

    Alternately, a :class:`VarName` object can be created with exactly four positional parameters:

    :param int phasenum: The number for the phase
    :param int histnum: The number for the histogram
    :param str varname: a single value can be used to create a :class:`VarName`
    :param int atomnum: The number for the atom
    
    '''
    import re
    IDdict = {}
    IDdict['phases'] = {}
    IDdict['hists'] = {}
    IDdict['atoms'] = {}
    # This dictionary lists descriptions for GSAS-II variables.
    # Note that keys may contain regular expressions, where '[xyz]'
    # matches 'x' 'y' or 'z' (equivalently '[x-z]' describes this as range of values).
    # '.*' matches any string
    VarDesc = {
        # Phase vars (p::<var>)
        'A[0-5]' : 'Reciprocal metric tensor component',
        'Vol' : 'Unit cell volume????',
        # Atom vars (p::<var>:a)
        'dA[xyz]' : 'change to atomic position',
        'AUiso':'Atomic isotropic displacement parameter',
        'AU[123][123]':'Atomic anisotropic displacement parameter',
        'AFrac': 'Atomic occupancy parameter',
        # Hist & Phase (HAP) vars (p:h:<var>)
        'Bab[AU]': 'Babinet solvent scattering coef.',
        'D[123][123]' : 'Anisotropic strain coef.',
        'Extinction' : 'Extinction coef.',
        'MD' : 'March-Dollase coef.',
        'Mustrain;.*' : 'Microstrain coef.',
        'Scale' : 'Phase scale factor',
        'Size;.*' : 'Crystallite size value',
        'eA' : '?',
        #Histogram vars (:h:<var>)
        'Absorption' : 'Absorption coef.',
        'Displace[XY]' : 'Debye-Scherrer sample displacement',
        'Lam' : 'Wavelength',
        'Polariz' : 'Polarization correction',
        'SH/L' : 'FCJ peak asymmetry correction',
        'Scale' : 'Histogram scale factor',
        '[UVW]' : 'Gaussian instrument broadening',
        '[XY]' : 'Cauchy instrument broadening',
        'Zero' : 'Debye-Scherrer zero correction',
        'nDebye' : 'Debye model background corr. terms',
        'nPeaks' : 'Fixed peak  background corr. terms',
        # Global vars (::<var>)
        }
    def __init__(self,*args):
        self.phase = None
        self.histogram = None
        self.name = ''
        self.atom = None
        if len(args) == 1:
            lst = args[0].split(':')
            raise Exception, "Need to look up IDs"
            self.phase = lst[0]
            self.histogram = lst[1]
            self.name = lst[2]
            if len(lst) > 3:
                self.atom = lst[3]
        elif len(args) == 4:
            self.phase = args[0]
            self.histogram = args[1]
            self.name = args[2]
            self.atom = args[3]
        else:
            raise Exception,"Incorrectly called GSAS-II parameter name"

    def __str__(self):
        return self.name()

    def name(self):
        '''Formats the GSAS-II variable name as a "traditional" string (p:h:<var>:a)

        :returns: the variable name as a str
        '''
        def _fmt(val):
            if val is None:
                return ""
            return str(val)
        return _fmt(self.phase) + ":" + _fmt(self.histogram) + _fmt(self.name) + _fmt(self.atom)

    def __repr__(self):
        '''Return the detailed contents of the object
        '''
        s = ""
        if self.phase:
            s += "Phase: " + str(self.phase) + "; "

        if self.histogram:
            s += "Histogram: " + str(self.histogram) + "; "
            
        if self.name:
            s += "Variable name: " + str(self.name) + "; "

        if self.atom:
            s += "Atom number: " + str(self.atom) + "; "

        return s+"("+self.name()+")"

    def getDescr(self):
        '''Return a short description for a GSAS-II variable 

        :returns: a short description or 'no definition' if not found
        '''
        # iterating over uncompiled regular expressions is not terribly fast,
        # but this routine should not need to be all that speedy
        for key in self.VarDesc:
            if re.match(key, self.name):
                return self.VarDesc[key]
        return 'no definition'

    def getDescr(self):
        '''Return a short description for a GSAS-II variable 

        :returns: a short description or 'no definition' if not found
        '''
        # iterating over uncompiled regular expressions is not terribly fast,
        # but this routine should not need to be all that speedy
        for key in self.VarDesc:
            if re.match(key, self.name):
                return self.VarDesc[key]
        return 'no definition'

    def fullDescr(self):
        '''Return a longer description for a GSAS-II variable 

        :returns: a short description or 'no definition' if not found
        '''
        # iterating over uncompiled regular expressions is not terribly fast,
        # but this routine should not need to be all that speedy
        str = self.name()
        
        for key in self.VarDesc:
            if re.match(key, self.name):
                return self.VarDesc[key]
        return 'no definition'


    def _show(self):
        'For testing, shows the current lookup table'
        print 'phases', self.IDdict['phases']
        print 'hists', self.IDdict['hists']
        print 'atomDict', self.IDdict['atoms']

