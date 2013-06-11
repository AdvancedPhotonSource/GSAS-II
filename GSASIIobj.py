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

Phase Tree Item
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
  \         AtomPtrs         ? (list)
  \         F000X            x-ray F(000) intensity (float)
  \         F000N            neutron F(000) intensity (float)
  \         Mydir            directory of current .gpx file (str)
  \         MCSA controls    ?
  \         Cell             List with 7 items: cell refinement flag (bool)
                             a, b, c, (Angstrom, float)
                             alpha, beta & gamma (degrees, float)
  \         Type             for now 'nuclear' (str)
  \         Map              dict of map parameters
  \         SH Texture       dict of spherical harmonic preferred orientation 
                             parameters
  \         Isotope          dict of isotopes for each atom type
  \         Isotopes         dict of scattering lengths for each isotope
                             combination for each element in phase  
  \         Name             phase name (str)
  \         SGData           Space group details, 
                             as defined in :mod:`GSASIIspc` as  :ref:`SGData definition <SGData_table>`
  \         Pawley neg wt    Restraint value for negative Pawley intensities
                             (float)
  \         Flip             Charge flip controls dict?
  \         Data plot type   ?
  \         Mass             Mass of unit cell contents in g/mol
  \         POhkl            March-Dollase preferred orientation direction
  \         Z                ?
  \         vdWRadii         ?
  \         Color            Colors for atoms (list of (r,b,g) triplets)
  \         AtomTypes        List of atom types
  \         AtomMass         List of masses for atoms
  \         doPawley         Flag for Pawley intensity extraction (bool)
  \         NoAtoms          Number of atoms per unit cell of each type (dict)
  \         Pawley dmin      maximum Q (as d-space) to use for Pawley 
                             extraction (float)
  \         BondRadii        Radius for each atom used to compute 
                             interatomic distances (list of floats)
  \         AngleRadii       Radius for each atom used to compute 
                             interatomic angles (list of floats)
ranId           \            unique random number Id for phase (int)
pId             \            ? (int)
Atoms           \            Atoms in phase as a list of lists. The outer list
                             is for each atom, the inner list contains 18
                             items:
                             0) atom label, 1) the atom type,
                             2) the refinement flags, 3-6) x, y, z, frac
                             7) site symmetry, 8) site multiplicity,
                             9) 'I' or 'A' for iso/anisotropic,
                             10) Uiso, 10-16) Uij, 16) unique Id #.
                             (list of lists)
Drawing         \            Display parameters (dict)
\           ballScale        Size of spheres in ball-and-stick display (float)
\           bondList         dict with bonds
\           contourLevel     ? (float)
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
\           atomPtrs         ? (list)
\           viewPoint        list of lists. First item in list is [x,y,z]
                             in fractional coordinates for the center of
                             the plot. Second item ?.
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
\           oldxy            ? (list with two floats)
\           cameraPos        Viewing position in A for plot (float)
\           depthFog         ? (bool)
RBModels        \            Rigid body assignments (note Rigid body definitions
                             are stored in their own main top-level tree entry.)
Pawley ref      \            Pawley reflections
Histograms      \            A dict of dicts. The key for the outer dict is
                             the histograms tied to this phase. The inner
                             dict contains the combined phase/histogram
                             parameters for items such as scale factors,
                             size and strain parameters. (dict)
MCSA            \            Monte-Carlo simulated annealing parameters
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
SGOps       symmetry operations as a list in form
            [[M1,T1],[M2,T2,...] where Mn is a 3x3 np.array
            and T is a length 3 np.array.
            Atom coordinates are transformed where the
            Asymmetric unit coordinates [X=(x,y,z)] are
            transformed using ``M*X+T ==> X'``
SGSys       symmetry unit cell: type one of
            'triclinic', 'monoclinic', 'orthorhombic',
            'tetragonal', 'rhombohedral', 'trigonal',
            'hexagonal', 'cubic' (str)
SGPolax     Axes for space group polarity. Will be one of
            '', 'x', 'y', 'x y', 'z', 'x z', 'y z',
            'xyz'. In the case where axes are arbitrary 
            '111' is used (P 1, and ?).
==========  ====================================================

'''

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

