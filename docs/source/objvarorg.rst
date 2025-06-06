=================================================
*GSASII Data object & variable organization*
=================================================

This chapter documents how data is organized within the data
structures used in GSAS-II. 

*Summary/Contents*
============================

.. contents:: Section Contents 

.. index::
   single: Parameter names
   single: GSAS-II variable naming

.. _VarNames_table:

Parameter names in GSAS-II
============================

Parameters in GSAS-II contain values that are used in diffraction
computations. These include atom positions, but also factors that
affect peak shapes or compensate for physical effects such as
absoption. 
Many, but not all, can be optimized. The term variables is intended to
refer to parameters that are being optimized in GSAS-II, but this
usage is not always applied consistently within GSAS-II and in places the term
variables may be applied to unvaried parameters. 

Parameter in GSAS-II are uniquely named using the following pattern,
``p:h:<var>:n``, where ``<var>`` is a variable name, as shown in the following table. Also,
``p`` is the phase number, ``h`` is the histogram number, 
and ``n`` is the atom parameter number
If a parameter does not depend on a histogram, phase or atom, ``h``, ``p`` and/or ``n`` will be omitted, 
so ``p::<var>:n``, ``:h:<var>`` and ``p:h:<var>`` are all valid names.

.. include:: ../source/vars.rst

.. _Constraints_table:

.. index::
   single: Constraints object description
   single: Data object descriptions; Constraints

GSAS-II Data Tree
============================

A GSAS-II project is stored in a data file and is loaded into a
wxPython data tree (wx.TreeCtrl) defined by
:class:`GSASIIctrlGUI.G2TreeCtrl`. Each entry in the tree has a text
label and a data object associated with it. Note that all information
used in a GSAS-II project is stored in the data tree, with the
exception of images (which are too large). For images, a reference to
the file location is saved and images are loaded from the file when
needed.

To save a GSAS-II project, routine :func:`GSASIImiscGUI.ProjFileSave` is
used to convert the tree contents to a "flat" format and write it to a
file. The tree is transversed, and for each first-level tree item, a
list is created, where the first item in that list is a two-element
list containing the label of the tree item and the data object
associated with the label. If there are second-level tree items that
are children of that first-level tree item, additional items are added
to the outermost list with pairs of text labels and data
objects. Finally the outermost list is converted to a binary
representation and written to disk with the Python pickle
function. Note that GSAS-II does not use any data tree items other
than first-level and second-level. 
Routine :func:`GSASIImiscGUI.ProjFileOpen` is used to read a GSAS-II
project file and populate the data tree. GSAS-II project files are
written with the ``.gpx`` extension.

Two pointers are kept for a selected tree entry in the GSAS-II data
tree, saved as class variables in :class:`GSASIIdataGUI.GSASII` (often
referenced as ``G2frame``). These are ``G2frame.PickId``, which points to
the selected data tree item, and ``G2frame.PatternId``, which points to
the parent of the data tree item, when ``G2frame.PickId`` points to a
histogram. The two pointer may be the same when the first-level tree
item for a histogram is selected. 

Constraints Tree Item
----------------------

Constraints are stored in a dict, separated into groups.
Note that parameter are named in the following pattern,
p:h:<var>:n, where p is the phase number, h is the histogram number
<var> is a variable name and n is the parameter number.
If a parameter does not depend on a histogram or phase or is unnumbered, that
number is omitted.
Note that the contents of each dict item is a List where each element in the
list is a :ref:`constraint definition objects <Constraint_definitions_table>`.
The constraints in this form are converted in
:func:`GSASIImapvars.ProcessConstraints` to the form used in :mod:`GSASIImapvars`

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
  * <var> a :class:`G2VarObj` object. (Note that in very old .gpx files this might be a str with a variable name of form 'p:h:name[:at]')

Note that the last three items in the list play a special role:

 * <fixedval> is the fixed value for a `constant equation` (``constype=c``)
   constraint or is None. For a `New variable` (``constype=f``) constraint,
   a variable name can be specified as a str (used for externally
   generated constraints)
 * <varyflag> is True or False for `New variable` (``constype=f``) constraints
   or is None. This indicates if this variable should be refined.
 * <constype> is one of four letters, 'e', 'c', 'h', 'f' that determines the type of constraint:

    * 'e' defines a set of equivalent variables. Only the first variable is refined (if the
      appropriate refine flag is set) and and all other equivalent variables in the list
      are generated from that variable, using the appropriate multipliers.
    * 'c' defines a constraint equation of form,
      :math:`m_1 \times var_1 + m_2 \times var_2 + ... = c`
    * 'h' defines a variable to hold (not vary). Any variable on this list is not varied,
      even if its refinement flag is set. Only one [mult,var] pair is allowed in a hold
      constraint and the mult value is ignored.
      This is of particular value when needing to hold one or more variables where a
      single flag controls a set of variables such as, coordinates,
      the reciprocal metric tensor or anisotropic displacement parameter.
    * 'f' defines a new variable (function) according to relationship
      :math:`newvar = m_1 \times var_1 + m_2 \times var_2 + ...`

.. _Covariance_table:

.. index::
   single: Covariance description
   single: Data object descriptions; Covariance

Covariance Tree Item
--------------------

The Covariance tree item has results from the last least-squares run. They
are stored in a dict with these keys:

.. tabularcolumns:: |l|l|p{4in}|

=============  ===============  ===========================================================================
  key            sub-key        explanation
=============  ===============  ===========================================================================
newCellDict    \                (dict) lattice parameters computed by
                                :func:`GSASIIstrMath.GetNewCellParms`
title          \                (str) Name of gpx file
variables      \                (list) Values for refined variables
                                (list of float values, length N,
                                ordered to match varyList)
sig            \                (list) Standard uncertainty values for refined variables
                                (list of float values, length N,
                                ordered to match varyList)
varyList       \                (list of str values, length N) List of directly refined variables
varyListStart  \                (list) initial refined variables before dependent vars are removed
newAtomDict    \                (dict) atom position values computed in
                                :func:`GSASIIstrMath.ApplyXYZshifts`
Lastshft       \                (list) The shifts applied to each variable in the last refinement
                                run. (list of float values, length N,
                                ordered to match varyList)
depSigDict     \                (dict) Values along with standard uncertainty values for
                                dependent variables
covMatrix      \                (np.array) The (NxN) covVariance matrix
freshCOV       \                (bool) indicates if the covMatrix has been freshly computed
msg            \                Warning/error messages from the last refinement run

Rvals          \                (dict) R-factors, GOF, Marquardt value for last
                                refinement cycle 
\              Nobs             (int) Number of observed data points 
\              Nvars            (int) Number of refined parameters
\              Rwp              (float) overall weighted profile R-factor (%)
\              chisq            (float) :math:`\sum w*(I_{obs}-I_{calc})^2`                                
                                for all data.
                                Note: this is what GSAS-II calls :math:`\chi^2`,
                                which is not the same thing as the reduced :math:`\chi^2`. 
\              lamMax           (float) Marquardt value applied to Hessian diagonal
\              GOF              (float) The goodness-of-fit, aka square root of
                                the reduced :math:`\chi^2` squared, after refinement. 
\              GOF0             (float) The goodness-of-fit, aka square root of
                                the reduced :math:`\chi^2` square, before refinement.
\              lastShifts       (dict) values of the shifts applied in the last
                                refinement cycle (note: differs from `Lastshft`,
                                which has values from the last run).
\              SVD0             (int) number of singular value decomposition
                                (SVD) singularities
\              converged        (bool) True if last refinement run converged
\              DelChi2          (float) change in :math:`\chi^2` in last refinement cycle
\              RestraintSum     (float) sum of restraints
\              RestraintTerms   (float) total number of restraints 
\              Max shft/sig     (float) maximum shift/s.u. for shifts applied in last refinement run
=============  ===============  ===========================================================================

.. _Phase_table:

.. index::
   single: Phase object description
   single: Data object descriptions; Phase

Phase Tree Items
----------------

Phase information is stored in the GSAS-II data tree as children of the
Phases item in a dict with keys:

.. tabularcolumns:: |l|l|p{4in}|

==========  ===============     =====================================================================================================
  key         sub-key           explanation
==========  ===============     =====================================================================================================
General         \               (dict) Overall information for the phase 
  \         3Dproj              (list of str) projections for 3D pole distribution plots
  \         AngleRadii          (list of floats) Default radius for each atom used to compute
                                interatomic angles 
  \         AtomMass            (list of floats) Masses for atoms
  \         AtomPtrs            (list of int) four locations (cx,ct,cs & cu) to use to pull info
                                from the atom records 
  \         AtomTypes           (llist of str) Atom types
  \         BondRadii           (list of floats) Default radius for each atom used to compute
                                interatomic distances 
  \         Cell                Unit cell parameters & ref. flag
                                (list with 8 items. All but first item are float.)

                                 0: cell refinement flag (True/False),

                                 1-3: a, b, c, (:math:`\AA`)

                                 4-6: alpha, beta & gamma, (degrees)

                                 7: volume (:math:`\AA^3`)

  \         Color               (list of (r,b,g) triplets) Colors for atoms 
  \         Compare             (dict) Polygon comparison parameters
  \         Data plot type      (str) data plot type ('Mustrain', 'Size' or
                                'Preferred orientation') for powder data 
  \         DisAglCtls          (dDict) with distance/angle search controls,
                                which has keys 'Name', 'AtomTypes',
                                'BondRadii', 'AngleRadii' which are as above
                                except are possibly edited. Also contains
                                'Factors', which is a 2 element list with
                                a multiplier for bond and angle search range
                                [typically (0.85,0.85)].
  \         F000X               (float) x-ray F(000) intensity 
  \         F000N               (float) neutron F(000) intensity 
  \         Flip                (dict) Charge flip controls
  \         HydIds              (dict) geometrically generated hydrogen atoms
  \         Isotope             (dict) Isotopes for each atom type
  \         Isotopes            (dict) Scattering lengths for each isotope
                                combination for each element in phase
  \         MCSA controls       (dict) Monte Carlo-Simulated Annealing controls 
  \         Map                 (dict) Map parameters
  \         Mass                (float) Mass of unit cell contents in g/mol
  \         Modulated           (bool) True if phase modulated
  \         Mydir               (str) Directory of current .gpx file 
  \         Name                (str) Phase name 
  \         NoAtoms             (dict) Number of atoms per unit cell of each type 
  \         POhkl               (list) March-Dollase preferred orientation direction
  \         Pawley dmin         (float) maximum Q (as d-space) to use for Pawley extraction 
  \         Pawley dmax         (float) minimum Q (as d-space) to use for Pawley extraction 
  \         Pawley neg wt       (float) Restraint value for negative Pawley intensities
  \         SGData              (object) Space group details as a 
                                :ref:`space group (SGData) <SGData_table>` 
                                object, as defined in :func:`GSASIIspc.SpcGroup`.
  \         SH Texture          (dict) Spherical harmonic preferred orientation parameters
  \         Super               (int) dimension of super group (0,1 only)
  \         Type                (str) phase type (e.g. 'nuclear')
  \         Z                   (dict) Atomic numbers for each atom type
  \         doDysnomia          (bool) flag for max ent map modification via Dysnomia
  \         doPawley            (bool) Flag for Pawley intensity extraction
  \         vdWRadii            (dict) Van der Waals radii for each atom type
ranId           \               (int) unique random number Id for phase 
pId             \               (int) Phase Id number for current project.
Atoms           \               (list of lists) Atoms in phase as a list of lists. The outer list
                                is for each atom, the inner list contains varying
                                items depending on the type of phase, see
                                the :ref:`Atom Records <Atoms_table>` description.
Drawing         \               (dict) Display parameters 
\           Atoms               (list of lists) with an entry for each atom that is drawn
\           Plane               (list) Controls for contour density plane display
\           Quaternion          (4 element np.array) Viewing quaternion 
\           Zclip               (float) clipping distance in :math:`\AA`
\           Zstep               (float) Step to de/increase Z-clip 
\           atomPtrs            (list) positions of x, type, site sym, ADP flag in Draw Atoms 
\           backColor           (list) background for plot as and R,G,B triplet
                                (default = [0, 0, 0], black).
\           ballScale           (float) Radius of spheres in ball-and-stick display 
\           bondList            (dict) Bonds
\           bondRadius          (float) Radius of binds in :math:`\AA` 
\           cameraPos           (float) Viewing position in :math:`\AA` for plot 
\           contourLevel        (float) map contour level in :math:`e/\AA^3` 
\           contourMax          (float) map contour maximum
\           depthFog            (bool) True if use depthFog on plot - set currently as False 
\           ellipseProb         (float) Probability limit for display of thermal
                                ellipsoids in % .
\           magMult             (float) multiplier for magnetic moment arrows
\           mapSize             (float) x & y dimensions of contourmap (fixed internally)
\           modelView           (4,4 array) from openGL drawing transofmation matrix
\           oldxy               (list with two floats) previous view point 
\           radiusFactor        (float) Distance ratio for searching for bonds. Bonds
                                are located that are within r(Ra+Rb) and (Ra+Rb)/r
                                where Ra and Rb are the atomic radii.
\           selectedAtoms       (list of int values) List of selected atoms 
\           showABC             (bool) Flag to show view point triplet. True=show.
\           showHydrogen        (bool) Flag to control plotting of H atoms.
\           showRigidBodies     (bool) Flag to highlight rigid body placement
\           showSlice           (bool) flag to show contour map
\           sizeH               (float) Size ratio for H atoms 
\           unitCellBox         (bool) Flag to control display of the unit cell.
\           vdwScale            (float) Multiplier of van der Waals radius for display of vdW spheres.
\           viewDir             (np.array with three floats) cartesian viewing direction 
\           viewPoint           (list of lists) First item in list is [x,y,z]
                                in fractional coordinates for the center of
                                the plot. Second item list of previous & current
                                atom number viewed (may be [0,0])
ISODISTORT      \               (dict) contains controls for running ISODISTORT and results from it
\           ISOmethod           (int) ISODISTORT method (currently 1 or 4; 2 & 3 not implemented in GSAS-II)
\           ParentCIF           (str) parent cif file name for ISODISTORT method 4
\           ChildCIF            (str) child cif file name for ISODISTORT method 4
\           SGselect            (dict) selection list for lattice types in radio result from ISODISTORT method 1
\           selection           (int) chosen selection from radio
\           radio               (list) results from ISODISTORT method 1
\           ChildMatrix         (3x3 array) transformation matrix for method 3 (not currently used)
\           ChildSprGp          (str) child space group for method 3 (not currently used)
\           ChildCell           (str) cell ordering for nonstandard orthorhombic ChildSprGrp in method 3 (not currently used)
\           G2ModeList          (list) ISODISTORT mode names
\           modeDispl           (list) distortion mode values; refinable parameters
\           ISOmodeDispl        (list) distortion mode values as determined in method 4 by ISODISTORT
\           NormList            (list) ISODISTORT normalization values; to convert mode value to fractional coordinate dsplacement
\           G2parentCoords      (list) full set of parent structure coordinates transformed to child structure; starting basis for mode displacements
\           G2VarList           (list) 
\           IsoVarList          (list)
\           G2coordOffset       (list) only adjustible set of parent structure coordinates
\           G2OccVarList        (list) 
\           Var2ModeMatrix      (array) atom variable to distortion mode transformation 
\           Mode2VarMatrix      (array) distortion mode to atom variable transformation
\           rundata             (dict) saved input information for use by ISODISTORT method 1

RBModels        \               Rigid body assignments (note Rigid body definitions
                                are stored in their own main top-level tree entry.)
RMC             \               (dict) RMCProfile, PDFfit & fullrmc controls
Pawley ref      \               (list) Pawley reflections
Histograms      \               (dict of dicts) The key for the outer dict is
                                the histograms tied to this phase. The inner
                                dict contains the combined phase/histogram
                                parameters for items such as scale factors,
                                size and strain parameters. The following are the
                                keys to the inner dict. (dict)
\           Babinet             (dict) For protein crystallography. Dictionary with two
                                entries, 'BabA', 'BabU'
\           Extinction          (list of float, bool) Extinction parameter 
\           Flack               (list of [float, bool]) Flack parameter & refine flag
\           HStrain             (list of two lists) Hydrostatic strain. The first is
                                a list of the HStrain parameters (1, 2, 3, 4, or 6
                                depending on unit cell), the second is a list of boolean
                                refinement parameters (same length)
\           Histogram           (str) The name of the associated histogram 
\           Layer Disp          (list of [float, bool]) Layer displacement in beam direction & refine flag
\           LeBail              (bool) Flag for LeBail extraction 
\           Mustrain            (list) Microstrain parameters, in order:
    
                                0. Type, one of  u'isotropic', u'uniaxial', u'generalized'
                                1. Isotropic/uniaxial parameters - list of 3 floats
                                2. Refinement flags - list of 3 bools
                                3. Microstrain axis - list of 3 ints, [h, k, l]
                                4. Generalized mustrain parameters - list of 2-6 floats, depending on space group
                                5. Generalized refinement flags - list of bools, corresponding to the parameters of (4)
\           Pref.Ori.           (list) Preferred Orientation. List of eight parameters.
                                Items marked SH are only used for Spherical Harmonics.
                                
                                0. (str) Type, 'MD' for March-Dollase or 'SH' for Spherical Harmonics
                                1. (float) Value 
                                2. (bool) Refinement flag 
                                3. (list) Preferred direction, list of ints, [h, k, l]
                                4. (int) SH - number of terms
                                5. (dict) SH -  
                                6. (list) SH
                                7. (float) SH 
\           Scale               (list of [float, bool]) Phase fraction & refine flag
\           Size                List of crystallite size parameters, in order:

                                0. (str) Type, one of  u'isotropic', u'uniaxial', u'ellipsoidal'
                                1. (list) Isotropic/uniaxial parameters - list of 3 floats
                                2. (list) Refinement flags - list of 3 bools
                                3. (list) Size axis - list of 3 ints, [h, k, l]
                                4. (list) Ellipsoidal size parameters - list of 6 floats
                                5. (list) Ellipsoidal refinement flags - list of bools, corresponding to the parameters of (4)
\           Use                 (bool) True if this histogram is to be used in refinement
MCSA            \               (dict) Monte-Carlo simulated annealing parameters 
==========  ===============     =====================================================================================================

.. _RBData_table:

.. index::
   single: Rigid Body Data description
   single: Data object descriptions; Rigid Body Data

Rigid Body Objects
------------------

Rigid body descriptions are available for two types of rigid bodies: 'Vector'
and 'Residue'. Vector rigid bodies are developed by a sequence of translations each
with a refinable magnitude and Residue rigid bodies are described as Cartesian coordinates
with defined refinable torsion angles.

.. tabularcolumns:: |l|l|p{4in}|

==========  ===============     ====================================================
  key         sub-key           explanation
==========  ===============     ====================================================
Vector      RBId                (dict of dict) vector rigid bodies 
\           AtInfo              (dict) Drad, Color: atom drawing radius & color for each atom type 
\           RBname              (str) Name assigned by user to rigid body 
\           VectMag             (list) vector magnitudes in :math:`\AA` 
\           rbXYZ               (list of 3 float Cartesian coordinates for Vector rigid body )
\           rbRef               (list of 3 int & 1 bool) 3 assigned reference atom nos. in rigid body for origin
                                definition, use center of atoms flag 
\           VectRef             (list of bool refinement flags for VectMag values )
\           rbTypes             (list of str) Atom types for each atom in rigid body 
\           rbVect              (list of lists) Cartesian vectors for each translation used to build rigid body 
\           useCount            (int) Number of times rigid body is used in any structure 
Residue     RBId                (dict of dict) residue rigid bodies 
\           AtInfo              (dict) Drad, Color: atom drawing radius & color for each atom type
\           RBname              (str) Name assigned by user to rigid body 
\           rbXYZ               (list of 3 float) Cartesian coordinates for Residue rigid body 
\           rbTypes             (list of str) Atom types for each atom in rigid body 
\           atNames             (list of str) Names of each atom in rigid body (e.g. C1,N2...) 
\           rbRef               (list of 3 int & 1 bool) 3 assigned reference atom nos. in rigid body for origin
                                definition, use center of atoms flag 
\           rbSeq               (list) Orig,Piv,angle,Riding : definition of internal rigid body
                                torsion; origin atom (int), pivot atom (int), torsion angle (float),
                                riding atoms (list of int)
\           SelSeq              (int,int) used by SeqSizer to identify objects
\           useCount            (int)Number of times rigid body is used in any structure 
RBIds           \               (dict) unique Ids generated upon creation of each rigid body 
\           Vector              (list) Ids for each Vector rigid body 
\           Residue             (list) Ids for each Residue rigid body 
==========  ===============     ====================================================

.. _SGData_table:

.. index::
   single: Space Group Data description
   single: Data object descriptions; Space Group Data

Space Group Objects
-------------------

Space groups are interpreted by :func:`GSASIIspc.SpcGroup`
and the information is placed in a SGdata object
which is a dict with these keys. Magnetic ones are marked "mag"

.. tabularcolumns:: |l|p{4.5in}|

==========  ========================================================================================
  key         explanation
==========  ========================================================================================
BNSlattsym  mag - (str) BNS magnetic space group symbol and centering vector
GenFlg      mag - (list) symmetry generators indices
GenSym      mag - (list) names for each generator
MagMom      mag - (list) "time reversals" for each magnetic operator
MagPtGp     mag - (str) Magnetic point group symbol
MagSpGrp    mag - (str) Magnetic space group symbol
OprNames    mag - (list) names for each space group operation
SGCen       (np.array) Symmetry cell centering vectors. A (n,3) np.array
            of centers. Will always have at least one row: ``np.array([[0, 0, 0]])``
SGFixed     (bool) Only True if phase mported from a magnetic cif file
            then the space group can not be changed by the user because 
            operator set from cif may be nonstandard
SGGen       (list) generators
SGGray      (bool) True if space group is a gray group (incommensurate magnetic structures)
SGInv       (bool) True if centrosymmetric, False if not
SGLatt      (str)Lattice centering type. Will be one of
            P, A, B, C, I, F, R 
SGLaue      (str) one of the following 14 Laue classes:
            -1, 2/m, mmm, 4/m, 4/mmm, 3R,
            3mR, 3, 3m1, 31m, 6/m, 6/mmm, m3, m3m
SGOps       (list) symmetry operations as a list of form
            ``[[M1,T1], [M2,T2],...]``
            where :math:`M_n` is a 3x3 np.array
            and :math:`T_n` is a length 3 np.array.
            Atom coordinates are transformed where the
            Asymmetric unit coordinates [X is (x,y,z)]
            are transformed using
            :math:`X^\prime = M_n*X+T_n`
SGPolax     (str) Axes for space group polarity. Will be one of
            '', 'x', 'y', 'x y', 'z', 'x z', 'y z',
            'xyz'. In the case where axes are arbitrary
            '111' is used (P 1, and ?).
SGPtGrp     (str) Point group of the space group
SGUniq      unique axis if monoclinic. Will be
            a, b, or c for monoclinic space groups.
            Will be blank for non-monoclinic.
SGSpin      mag - (list) of spin flip operatiors (+1 or -1) for the space group operations
SGSys       (str) symmetry unit cell: type one of
            'triclinic', 'monoclinic', 'orthorhombic',
            'tetragonal', 'rhombohedral', 'trigonal',
            'hexagonal', 'cubic' 
SSGK1       (list) Superspace multipliers
SpGrp       (str) space group symbol 
SpnFlp      mag - (list) Magnetic spin flips for every magnetic space group operator
==========  ========================================================================================

.. _SSGData_table:

.. index::
   single: Superspace Group Data description
   single: Data object descriptions; Superspace Group Data

Superspace groups [3+1] are interpreted by :func:`GSASIIspc.SSpcGroup`
and the information is placed in a SSGdata object
which is a dict with these keys:

.. tabularcolumns:: |l|p{4.5in}|

==========  ====================================================
  key         explanation
==========  ====================================================
SSGCen      (list) 4D cell centering vectors [0,0,0,0] at least
SSGK1       (list) Superspace multipliers
SSGOps      (list) 4D symmetry operations as [M,T] so that M*x+T = x'
SSpGrp      (str) superspace group symbol extension to space group
            symbol, accidental spaces removed
modQ        (list) modulation/propagation vector
modSymb     (list of str) Modulation symbols
==========  ====================================================

Phase Information
--------------------

.. _Phase_Information:

.. index::
   single: Phase information record description

Phase information is placed in one of the following keys:

.. tabularcolumns:: |l|p{4.5in}|

==========  ==============================================================
  key         explanation
==========  ==============================================================
General       Overall information about a phase
Histograms    Information about each histogram linked to the 
              current phase as well as parameters that 
              are defined for each histogram and phase
              (such as sample peak widths and preferred 
              orientation parameters. 
Atoms         Contains a list of atoms, as described in the 
              :ref:`Atom Records <Atoms_table>` description.
Drawing       Parameters that determine how the phase is 
              displayed, including a list of atoms to be 
              included, as described in the 
              :ref:`Drawing Atom Records <Drawing_atoms_table>`
              description
MCSA          Monte-Carlo simulated annealing parameters
pId           The index of each phase in the project, numbered
              starting at 0
ranId         An int value with a unique value for each phase 
RBModels      A list of dicts with parameters for each 
              rigid body inserted into the current phase, 
              as defined in the 
              :ref:`Rigid Body Insertions <Rigid_Body_Insertions>`.
              Note that the rigid bodies are defined as 
              :ref:`Rigid Body Objects <RBData_table>` 
RMC           PDF modeling parameters
Pawley ref    Pawley refinement parameters 

==========  ==============================================================

.. _Atoms_table:

.. index::
   single: Atoms record description
   single: Data object descriptions; Atoms record

--------------------
Atom Records
--------------------

If ``phasedict`` points to the phase information in the data tree, then
atoms are contained in a list of atom records (list) in
``phasedict['Atoms']``. Also needed to read atom information
are four pointers, ``cx,ct,cs,cia = phasedict['General']['AtomPtrs']``,
which define locations in the atom record, as shown below. Items shown are
always present; additional ones for macromolecular phases are marked 'mm', 
and those for magnetic structures are marked 'mg'

.. tabularcolumns:: |l|p{4.5in}|

==============      ====================================================
location            explanation
==============      ====================================================
ct-4                mm - (str) residue number
ct-3                mm - (str) residue name (e.g. ALA) 
ct-2                mm - (str) chain label 
ct-1                (str) atom label 
ct                  (str) atom type 
ct+1                (str) refinement flags; combination of 'F', 'X', 'U', 'M' 
cx,cx+1,cx+2        (3 floats) the x,y and z coordinates 
cx+3                (float) site occupancy 
cx+4,cx+5,cx+6      mg - (list) atom magnetic moment along a,b,c in Bohr magnetons 
cs                  (str) site symmetry 
cs+1                (int) site multiplicity 
cia                 (str) ADP flag: Isotropic ('I') or Anisotropic ('A')
cia+1               (float) Uiso 
cia+2...cia+7       (6 floats) U11, U22, U33, U12, U13, U23 
atom[cia+8]         (int) unique atom identifier 

==============      ====================================================

.. _Drawing_atoms_table:

.. index::
   single: Drawing atoms record description
   single: Data object descriptions; Drawing atoms record

----------------------------
Drawing Atom Records
----------------------------

If ``phasedict`` points to the phase information in the data tree, then
drawing atoms are contained in a list of drawing atom records (list) in
``phasedict['Drawing']['Atoms']``. Also needed to read atom information
are four pointers, ``cx,ct,cs,ci = phasedict['Drawing']['AtomPtrs']``,
which define locations in the atom record, as shown below. Items shown are
always present; additional ones for macromolecular phases are marked 'mm',
and those for magnetic structures are marked 'mg'

.. tabularcolumns:: |l|p{4.5in}|

==============   ===================================================================================
location            explanation
==============   ===================================================================================
ct-4                mm - (str) residue number 
ct-3                mm - (str) residue name (e.g. ALA) 
ct-2                mm - (str) chain label 
ct-1                (str) atom label
ct                  (str) atom type 
cx,cx+1,cx+2        (3 floats) the x,y and z coordinates 
cx+3,cx+4,cx+5      mg - (3 floats) atom magnetic moment along a,b,c in Bohr magnetons 
cs-1                (str) Sym Op symbol; sym. op number + unit cell id (e.g. '1,0,-1') 
cs                  (str) atom drawing style; e.g. 'balls & sticks' 
cs+1                (str) atom label style (e.g. 'name') 
cs+2                (int) atom color (RBG triplet) 
cs+3                (str) ADP flag: Isotropic ('I') or Anisotropic ('A')
cs+4                (float) Uiso 
cs+5...cs+11        (6 floats) U11, U22, U33, U12, U13, U23 
ci                  (int) unique atom identifier; matches source atom Id in Atom Records 
==============   ===================================================================================

.. _Rigid_Body_Insertions:

----------------------------
Rigid Body Insertions
----------------------------

If ``phasedict`` points to the phase information in the data tree, then
rigid body information is contained in list(s) in
``phasedict['RBModels']['Residue']`` and/or ``phasedict['RBModels']['Vector']``
for each rigid body inserted into the current phase. 

.. tabularcolumns:: |l|p{4.5in}|

==============   ===================================================================================
key              explanation
==============   ===================================================================================
fixOrig           Should the origin be fixed (when editing, not the refinement flag)
Ids               Ids for assignment of atoms in the rigid body
numChain          Chain number for macromolecular fits 
Orient            Orientation of the RB as a quaternion and a refinement flag (' ', 'A' or 'AV')
OrientVec         Orientation of the RB expressed as a vector and azimuthal rotation angle
Orig              Origin of the RB in fractional coordinates and refinement flag (bool)
RBId              References the unique ID of a rigid body in the 
                  :ref:`Rigid Body Objects <RBData_table>`
RBname            The name for the rigid body (str)
AtomFrac          The atom fractions for the rigid body
ThermalMotion     The thermal motion description for the rigid body, which includes a choice for 
                  the model and can include TLS parameters or an overall Uiso value. 
Torsions          Defines the torsion angle and refinement flag for each torsion defined in 
                  the :ref:`Rigid Body Object <RBData_table>`
==============   ===================================================================================

.. _Powder_table:

.. index::
   single: Powder data object description
   single: Data object descriptions; Powder Data

Powder Diffraction Tree Items
-----------------------------

Every powder diffraction histogram is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "PWDR ". The
diffraction data for that information are directly associated with
that tree item and there are a series of children to that item. The
routines :func:`GSASIIdataGUI.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |p{1in}|p{1in}|p{4in}|

======================     ===============  ===========================================================
  key                       sub-key          explanation
======================     ===============  ===========================================================
Comments                    \               (list of str) Text strings extracted from the original powder
                                            data header. These cannot be changed by the user;
                                            it may be empty.
Limits                      \               (list) two two element lists, as [[Ld,Hd],[L,H]]
                                            where L and Ld are the current and default lowest
                                            two-theta value to be used and
                                            where H and Hd are the current and default highest
                                            two-theta value to be used.
Reflection Lists            \               (dict of dicts) with an entry for each phase in the
                                            histogram. The contents of each dict item
                                            is a dict containing reflections, as described in
                                            the :ref:`Powder Reflections <PowderRefl_table>`
                                            description.
Instrument Parameters       \               (dict) The instrument parameters uses different dicts 
                                            for the constant wavelength (CW) and time-of-flight (TOF)
                                            cases. See below for the descriptions of each. 
wtFactor                    \               (float) A weighting factor to increase or decrease
                                            the leverage of data in the histogram .
                                            A value of 1.0 weights the data with their
                                            standard uncertainties and a larger value
                                            increases the weighting of the data (equivalent
                                            to decreasing the uncertainties).
Sample Parameters           \               (dict) Parameters that describe how
                                            the data were collected, as listed
                                            below. Refinable parameters are a list containing
                                            a float and a bool, where the second value
                                            specifies if the value is refined, otherwise
                                            the value is a float unless otherwise noted.
\                           Scale           The histogram scale factor (refinable)
\                           Absorption      The sample absorption coefficient as
                                            :math:`\mu r` where r is the radius
                                            (refinable). Only valid for Debye-Scherrer geometry.
\                           SurfaceRoughA   Surface roughness parameter A as defined by
                                            Surotti, *J. Appl. Cryst*, **5**, 325-331, 1972.
                                            (refinable - only valid for Bragg-Brentano geometry)
\                           SurfaceRoughB   Surface roughness parameter B (refinable -
                                            only valid for Bragg-Brentano geometry)
\                           DisplaceX,      Sample displacement from goniometer center
                            DisplaceY       where Y is along the beam direction and
                                            X is perpendicular. Units are :math:`\mu m`
                                            (refinable).
\                           Phi, Chi,       Goniometer sample setting angles, in degrees.
                            Omega
\                           Gonio. radius   Radius of the diffractometer in mm
\                           InstrName       (str) A name for the instrument, used in preparing
                                            a CIF .
\                           Force,          Variables that describe how the measurement
                            Temperature,    was performed. Not used directly in
                            Humidity,       any computations.
                            Pressure,
                            Voltage
\                           ranId           (int) The random-number Id for the histogram
                                            (same value as where top-level key is ranId)
\                           Type            (str) Type of diffraction data, may be 'Debye-Scherrer'
                                            or 'Bragg-Brentano' .
hId                         \               (int) The number assigned to the histogram when
                                            the project is loaded or edited (can change)
ranId                       \               (int) A random number id for the histogram
                                            that does not change
Background                  \               (list) The background is stored as a list with where
                                            the first item in the list is list and the second
                                            item is a dict. The list contains the background
                                            function and its coefficients; the dict contains
                                            Debye diffuse terms and background peaks.
                                            (TODO: this needs to be expanded.)
Data                        \               (list) The data consist of a list of 6 np.arrays
                                            containing in order:

                                            0. the x-postions (two-theta in degrees),
                                            1. the intensity values (Yobs),
                                            2. the weights for each Yobs value
                                            3. the computed intensity values (Ycalc)
                                            4. the background values
                                            5. Yobs-Ycalc
======================     ===============  ===========================================================

.. _CWPowder_table:

.. index::
   single: Powder data CW Instrument Parameters

-----------------------------
CW Instrument Parameters
-----------------------------

Instrument Parameters are placed in a list of two dicts, 
where the keys in the first dict are listed below. Note that the dict contents are different for 
constant wavelength (CW) vs. time-of-flight (TOF) histograms. 
The value for each item is a list containing three values: the initial value, the current value
and a refinement flag which can have a value of True, False or 0 where 0 indicates a value that
cannot be refined. The first and second values are floats unless otherwise noted.
Items not refined are noted as [*]

.. tabularcolumns:: |l|p{1in}|p{4in}|

========================    ===============  ===========================================================
  key                       sub-key           explanation
========================    ===============  ===========================================================
Instrument Parameters[0]    Type [*]            (str) Histogram type:
                                                * 'PXC' for constant wavelength x-ray
                                                * 'PNC' for constant wavelength neutron
\                           Bank [*]            (int) Data set number in a multidata file (usually 1)
\                           Lam                 (float) Specifies a wavelength in :math:`\AA`
\                           Lam1 [*]            (float) Specifies the primary wavelength in
                                                :math:`\AA`, used in place of Lam 
                                                when an :math:`\alpha_1, \alpha_2`
                                                source is used.
\                           Lam2 [*]            (float) Specifies the secondary wavelength in
                                                :math:`\AA`, used with Lam1
\                           I(L2)/I(L1)         (float) Ratio of Lam2 to Lam1, used with Lam1
\                           Zero                (float) Two-theta zero correction in *degrees*
\                           Azimuth [*]         (float) Azimuthal setting angle for data recorded with differing setting angles
\                           U, V, W             (float) Cagliotti profile coefficients
                                                for Gaussian instrumental broadening, where the
                                                FWHM goes as
                                                :math:`U \tan^2\theta + V \tan\theta + W`
\                           X, Y, Z             (float) Cauchy (Lorentzian) instrumental broadening coefficients
\                           SH/L                (float) Variant of the Finger-Cox-Jephcoat asymmetric
                                                peak broadening ratio. Note that this is the
                                                sum of S/L and H/L where S is
                                                sample height, H is the slit height and
                                                L is the goniometer diameter.
\                           Polariz.            (float) Polarization coefficient. 
Instrument Parameters[1]                        (empty dict)
========================    ===============  ===========================================================

.. _TOFPowder_table:

.. index::
   single: Powder data TOF Instrument Parameters

-----------------------------
TOF Instrument Parameters
-----------------------------

Instrument Parameters are also placed in a list of two dicts, 
where the keys in each dict listed below, but here for 
time-of-flight (TOF) histograms. 
The value for each item is a list containing three values: the initial value, the current value
and a refinement flag which can have a value of True, False or 0 where 0 indicates a value that
cannot be refined. The first and second values are floats unless otherwise noted.
Items not refined are noted as [*]

.. tabularcolumns:: |l|p{1.5in}|p{4in}|

========================    ===============  ===========================================================
  key                        sub-key          explanation
========================    ===============  ===========================================================
Instrument Parameters[0]    Type [*]            (str) Histogram type:
                                                * 'PNT' for time of flight neutron
\                           Bank                (int) Data set number in a multidata file
\                           2-theta [*]         (float) Nominal scattering angle for the detector
\                           fltPath [*]         (float) Total flight path source-sample-detector
\                           Azimuth [*]         (float) Azimuth angle for detector right hand rotation 
                                                from horizontal away from source
\                           difC,difA,          (float) Diffractometer constants for conversion of d-spacing to TOF
                            difB                in microseconds
\                           Zero                (float) Zero point offset (microseconds)
\                           alpha               (float) Exponential rise profile coefficients
\                           beta-0              (float) Exponential decay profile coefficients
                            beta-1
                            beta-q
\                           sig-0               (float) Gaussian profile coefficients
                            sig-1
                            sig-2
                            sig-q    
\                           X,Y,Z               (float) Lorentzian profile coefficients
Instrument Parameters[1]    Pdabc               (list of 4 float lists) Originally created for use in gsas as optional tables 
                                                of d, alp, bet, d-true; for a reflection alpha & beta are obtained via interpolation
                                                from the d-spacing and these tables. The d-true column is apparently unused.
========================    ===============  ===========================================================


.. _PowderRefl_table:

.. index::
   single: Powder reflection object description
   single: Data object descriptions; Powder Reflections

Powder Reflection Data Structure
--------------------------------

The data tree entry for powder diffraction histograms contains an
entry labeled ``Reflection Lists`` containing a dict keyed by phase
name, for every phase linked to the histogram. Each entry is itself a
dict with four entries, with keys:

==========  ====================================================
  key         explanation
==========  ====================================================
RefList      This contains the reflection list, as described 
             below.
FF           Contains a dict with two entries, 
             ``El`` which contains a list of ``n`` 
             element types and
             ``FF`` which contains a 55 x ``n`` np.array of 
             of form factor values.
Type         Contains a string specifying the type of 
             histogram, such as 'PXC'
Super        Contains a bool value, which is True when
             the phase has a superspace spacegroup (3+1 
             dimension).
==========  ====================================================

one element of which is `'RefList'`, which is a np.array containing
reflections. The columns in that array are documented below.

==========  ====================================================
  index         explanation
==========  ====================================================
 0,1,2           h,k,l 
 3               multiplicity
 4               d-space, :math:`\AA`
 5               pos, two-theta
 6               sig, Gaussian width
 7               gam, Lorenzian width
 8               :math:`F_{obs}^2`
 9               :math:`F_{calc}^2`
 10              reflection phase, in degrees
 11              intensity correction for reflection, this 
                 times :math:`F_{obs}^2` or :math:`F_{calc}^2` 
                 gives Iobs or Icalc
 12              Preferred orientation correction
 13              Transmission (absorption correction)
 14              Extinction correction
==========  ====================================================

Note that when the ``Super`` entry in the phase's main dict is True,
indicating that the phase is a 3+1 super-space group, the columns are: 

==========  ====================================================
  index         explanation
==========  ====================================================
 0,1,2,3         h,k,l,m
 4               multiplicity
 5               d-space, :math:`\AA`
 6               pos, two-theta
 7               sig, Gaussian width
 8               gam, Lorenzian width
 9               :math:`F_{obs}^2`
 10              :math:`F_{calc}^2`
 11              reflection phase, in degrees
 12              intensity correction for reflection, this
                 times :math:`F_{obs}^2` or :math:`F_{calc}^2`
                 gives Iobs or Icalc
 13              Preferred orientation correction
 14              Transmission (absorption correction)
 15              Extinction correction
==========  ====================================================


.. _Xtal_table:

.. index::
   single: Single Crystal data object description
   single: Data object descriptions; Single crystal data

Single Crystal Tree Items
-------------------------

Every single crystal diffraction histogram is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "HKLF ". The
diffraction data for that information are directly associated with
that tree item and there are a series of children to that item. The
routines :func:`GSASIIdataGUI.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ===============     ====================================================
  key                      sub-key          explanation
======================  ===============     ====================================================
Data                        \               (dict) that contains the
                                            reflection table,
                                            as described in the
                                            :ref:`Single Crystal Reflections
                                            <XtalRefl_table>`
                                            description.

Instrument Parameters       \               (list) containing two dicts where the possible
                                            keys in each dict are listed below. The value
                                            for most items is a list containing two values:
                                            the initial value, the current value.
                                            The first and second
                                            values are floats unless otherwise noted.
\                           Lam             (two floats) Specifies a wavelength in :math:`\AA` 
\                           Type            (two str values) Histogram type :
                                            * 'SXC' for constant wavelength x-ray
                                            * 'SNC' for constant wavelength neutron
                                            * 'SNT' for time of flight neutron
                                            * 'SEC' for constant wavelength electrons (e.g. micro-ED)
\                           InstrName       (str) A name for the instrument, used in preparing a CIF
wtFactor                    \               (float) A weighting factor to increase or decrease
                                            the leverage of data in the histogram.
                                            A value of 1.0 weights the data with their
                                            standard uncertainties and a larger value
                                            increases the weighting of the data (equivalent
                                            to decreasing the uncertainties).

hId                         \               (int) The number assigned to the histogram when
                                            the project is loaded or edited (can change)
ranId                       \               (int) A random number id for the histogram
                                            that does not change
======================  ===============     ====================================================

.. _XtalRefl_table:

.. index::
   single: Single Crystal reflection object description
   single: Data object descriptions; Single Crystal Reflections

Single Crystal Reflection Data Structure
----------------------------------------

For every single crystal a histogram, the ``'Data'`` item contains
the structure factors as an np.array in item `'RefList'`.
The columns in that array are documented below for
non-superspace phases. 

.. tabularcolumns:: |l|l|p{4in}|

==========  ==========  ====================================================
  index     3+1 index    explanation
==========  ==========  ====================================================
 0,1,2        0,1,2      reflection indices, h,k,l 
 \              3        3+1 superspace index, m
 3              4        flag (0 absent, 1 observed)
 4              5        d-space, :math:`\AA`
 5              6        :math:`F_{obs}^2`
 6              7        :math:`\sigma(F_{obs}^2)`
 7              8        :math:`F_{calc}^2`
 8              9        :math:`F_{obs}^2(T)`
 9             10        :math:`F_{calc}^2(T)`
 10            11        reflection phase, in degrees
 11            12        intensity correction for reflection, this times
                         :math:`F_{obs}^2` or :math:`F_{calc}^2`
                         gives Iobs or Icalc
==========  ==========  ====================================================

Notes:

  * The annotation "(T)" in the second set of :math:`F^2(T)` values 
    stands for "true," where the values are on an absolute scale 
    through application of the scale factor. 
  * The left-most column gives the entry index for three dimensional 
    spacegroups, the column to the right of that has the index for 
    3+1 superspace phases, where there are four reflection indices 
    h, k, l, m. 

.. _Image_table:

.. index::
   image: Image data object description
   image: Image object descriptions

Image Data Structure
--------------------

Every 2-dimensional image is stored in the GSAS-II data tree
with a top-level entry named beginning with the string "IMG ". The
image data are directly associated with that tree item and there
are a series of children to that item. The routines :func:`GSASIIdataGUI.GSASII.GetUsedHistogramsAndPhasesfromTree`
and :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` will
load this information into a dictionary where the child tree name is
used as a key, and the information in the main entry is assigned
a key of ``Data``, as outlined below.

.. tabularcolumns:: |l|l|p{4in}|

======================  ======================  ====================================================
  key                      sub-key              explanation
======================  ======================  ====================================================
Comments                    \                   (list of str) Text strings extracted from the original image data
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
\                           wavelength          (float) The radiation wavelength (:math:`\AA`) as entered by the user 
                                                (or someday obtained from the image header).
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
\                           SpotMask            (dict: int & array)
                                                'esdMul'(int) number of standard deviations above mean ring intensity to mask
                                                'spotMask' (bool array) the spot mask for every pixel in image         

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

.. _parmDict_table:

.. index::
   single: Parameter dictionary

Parameter Dictionary
=========================

The parameter dictionary contains all of the variable parameters for the refinement.
The dictionary keys are the name of the parameter (<phase>:<hist>:<name>:<atom>).
It is prepared in two ways. When loaded from the tree
(in :meth:`GSASIIdataGUI.GSASII.MakeLSParmDict` and
:meth:`GSASIIfiles.ExportBaseclass.loadParmDict`),
the values are lists with two elements: ``[value, refine flag]``

When loaded from the GPX file (in
:func:`GSASIIstrMain.Refine` and :func:`GSASIIstrMain.SeqRefine`), the value in the
dict is the actual parameter value (usually a float, but sometimes a
letter or string flag value (such as I or A for iso/anisotropic).

Texture implementation
==============================

There are two different places where texture can be treated in GSAS-II. 
One is for mitigating the effects of texture in a structural refinement.
The other is for texture characterization. 

For reducing the effect of texture in a structural refinement
there are entries labeled preferred orientation in each phase's
data tab. Two different approaches can be used for this, the March-Dollase 
model and spherical harmonics.

For the March-Dollase model, one axis in reciprocal space is designated as 
unique (defaulting to the 001 axis) and reflections are corrected 
according to the angle they make with this axis depending on 
the March-Dollase ratio. (If unity, no correction is made). 
The ratio can be greater than one or less than one depending on if 
crystallites oriented along the designated axis are 
overrepresented or underrepresented. For most crystal systems there is an 
obvious choice for the direction of the unique axis and then only a single
term needs to be refined. If the number is close to 1, then the correction 
is not needed. 

The second method for reducing the effect of texture in a structural 
refinement is to create a crystallite orientation probability surface as an 
expansion in terms spherical harmonic functions. Only functions consistent with 
cylindrical diffraction suymmetry and having texture symmetry 
consistent with the Laue class of phase are used and are allowed, 
so the higher the symmetry the fewer terms that are available for a given spherical harmonics order. 
To use this correction, select the lowest order that provides 
refinable terms and perform a refinement. If the texture index remains close to 
one, then the correction is not needed. If a significant improvement is 
noted in the profile Rwp, one may wish to see if a higher order expansion
gives an even larger improvement. 

To characterize texture in a material, generally one needs data collected with the 
sample at multiple orientations or, for TOF, with detectors at multiple 
locations around the sample. In this case the detector orientation is given in 
each histogram's Sample Parameters and the sample's orientation is described 
with the Euler angles specifed on the phase's Texture tab, which is also 
where the texture type (cylindrical, rolling,...) and the spherical 
harmonic order is selected. This should not be used with a single dataset and 
should not be used if the preferred orientations corrections are used. 

The coordinate system used for texture characterization is defined where 
the sample coordinates (Psi, gamma) are defined with an instrument coordinate
system (I, J, K) such that K is normal to the diffraction plane and J is coincident with the
direction of the incident radiation beam toward the source. We further define
a standard set of right-handed goniometer eulerian angles (Omega, Chi, Phi) so that Omega and Phi are
rotations about K and Chi is a rotation about J when Omega = 0. Finally, as the sample
may be mounted so that the sample coordinate system (Is, Js, Ks) does not coincide with
the instrument coordinate system (I, J, K), we define three eulerian sample rotation angles
(Omega-s, Chi-s, Phi-s) that describe the rotation from (Is, Js, Ks) to (I, J, K). The sample rotation
angles are defined so that with the goniometer angles at zero Omega-s and Phi-s are rotations
about K and Chi-s is a rotation about J.

Three typical examples:

    1) Bragg-Brentano laboratory diffractometer: Chi=0
    2) Debye-Scherrer counter detector; sample capillary axis perpendicular to diffraction plane: Chi=90
    3) Debye-Scherrer 2D area detector positioned directly behind sample; sample capillary axis horizontal; Chi=0

NB: The area detector azimuthal angle will equal 0 in horizontal plane to right as viewed from x-ray source and will equal 
90 at vertical "up" direction.
            
ISODISTORT implementation
==============================

CIFs prepared with the ISODISTORT web site 
https://stokes.byu.edu/iso/isodistort_version5.6.1/isodistort.php
[B. J. Campbell, H. T. Stokes, D. E. Tanner, and D. M. Hatch, "ISODISPLACE: An Internet Tool for Exploring Structural Distortions." 
J. Appl. Cryst. 39, 607-614 (2006).] can be read into GSAS-II using import CIF. This will cause constraints to be established for 
structural distortion modes read from the CIF. At present, of the five types of modes  only displacive(``_iso_displacivemode``...) 
and occupancy (``_iso_occupancymode``...) are processed. Not yet processed: ``_iso_magneticmode``..., 
``_iso_rotationalmode``... & ``_iso_strainmode``...

The CIF importer :mod:`G2phase_CIF` implements class :class:`G2phase_CIF.CIFPhaseReader` which offers two methods associated 
with ISODISTORT (ID) input. Method :meth:`G2phase_CIF.CIFPhaseReader.ISODISTORT_test` checks to see if a CIF block contains 
the loops with ``_iso_displacivemode_label`` or  ``_iso_occupancymode_label`` items. If so, method 
:meth:`G2phase_CIF.CIFPhaseReader.ISODISTORT_proc` is called to read and interpret them. The results are placed into the 
reader object's ``.Phase`` class variable as a dict item with key ``'ISODISTORT'``. 

Note that each mode ID has a long label with a name such as  Pm-3m[1/2,1/2,1/2]R5+(a,a,0)[La:b:dsp]T1u(a). Function 
:func:`G2phase_CIF.ISODISTORT_shortLbl` is used to create a short name for this, such as R5_T1u(a) which is made unique 
by addition of _n if the short name is duplicated. As each mode is processed, a constraint corresponding to that mode is 
created and is added to list in the reader object's ``.Constraints`` class variable. Items placed into that list can either 
be a list, which corresponds to a function (new var) type :ref:`constraint definition <Constraints_table>` entry, or an item 
can be a dict, which provides help information for each constraint.

Displacive modes
------------------------------

The coordinate variables, as named by ISODISTORT, are placed in ``.Phase['ISODISTORT']['IsoVarList']`` and the 
corresponding :class:`GSASIIobj.G2VarObj` objects for each are placed in ``.Phase['ISODISTORT']['G2VarList']``. 
The mode variables, as named by ISODISTORT, are placed in ``.Phase['ISODISTORT']['IsoModeList']`` and the 
corresponding :class:`GSASIIobj.G2VarObj` objects for each are placed in ``.Phase['ISODISTORT']['G2ModeList']``.
[Use ``str(G2VarObj)`` to get the variable name from the G2VarObj object, but note that the phase number, *n*, for the prefix 
"*n*::" cannot be determined as the phase number is not yet assigned.]

Displacive modes are a bit complex in that they relate to delta displacements, relative to an offset value for each coordinate, 
and because the modes are normalized. While GSAS-II also uses displacements,  these are added to the coordinates after 
each refinement cycle and then the delta values are set to zero. 
ISODISTORT uses fixed offsets (subtracted from the actual position 
to obtain the delta values) that are taken from the parent structure coordinate and the initial offset value 
(in ``_iso_deltacoordinate_value``) and these are placed in 
``.Phase['ISODISTORT']['G2coordOffset']`` in the same order as ``.Phase['ISODISTORT']['G2ModeList']``, 
``.Phase['ISODISTORT']['IsoVarList']`` and ''.Phase[ISODISTORT']['G2parentCoords']''.' 

The normalization factors (which the delta values are divided by) 
are taken from ``_iso_displacivemodenorm_value`` and are placed in ``.Phase['ISODISTORT']['NormList']`` in the same 
order as as ``...['IsoModeList']`` and ``...['G2ModeList']``.

The CIF contains a sparse matrix, from the ``loop_`` containing ``_iso_displacivemodematrix_value`` which provides the equations 
for determining the mode values from the coordinates, that matrix is placed in ``.Phase['ISODISTORT']['Mode2VarMatrix']``. 
The matrix is inverted to produce ``.Phase['ISODISTORT']['Var2ModeMatrix']``, which determines how to compute the
mode values from the delta coordinate values. These values are used for the in :func:`GSASIIconstrGUI.ShowIsoDistortCalc`,
which shows coordinate and mode values, the latter with s.u. values. 

Occupancy modes
------------------------------


The delta occupancy variables, as named by ISODISTORT, are placed in 
``.Phase['ISODISTORT']['OccVarList']`` and the corresponding :class:`GSASIIobj.G2VarObj` objects for each are placed 
in ``.Phase['ISODISTORT']['G2OccVarList']``. The mode variables, as named by ISODISTORT, are placed in 
``.Phase['ISODISTORT']['OccModeList']`` and the corresponding :class:`GSASIIobj.G2VarObj` objects for each are placed 
in ``.Phase['ISODISTORT']['G2OccModeList']``.

Occupancy modes, like Displacive modes, are also refined as delta values.  However, GSAS-II directly refines the fractional 
occupancies. Offset values for each atom, are taken from ``_iso_occupancy_formula`` and are placed in 
``.Phase['ISODISTORT']['ParentOcc]``. (Offset values are subtracted from the actual position to obtain the delta values.) 
Modes are normalized (where the mode values are divided by the normalization factor) are taken from ``_iso_occupancymodenorm_value`` 
and are placed in ``.Phase['ISODISTORT']['OccNormList']`` in the same order as as ``...['OccModeList']`` and
``...['G2OccModeList']``.

The CIF contains a sparse matrix, from the ``loop_`` containing ``_iso_occupancymodematrix_value``, which provides the 
equations for determining the mode values from the coordinates. That matrix is placed in ``.Phase['ISODISTORT']['Occ2VarMatrix']``. 
The matrix is inverted to produce ``.Phase['ISODISTORT']['Var2OccMatrix']``, which determines how to compute the
mode values from the delta coordinate values.


Mode Computations
------------------------------

Constraints are processed after the CIF has been read in :meth:`GSASIIdataGUI.GSASII.OnImportPhase` or  
:meth:`GSASIIscriptable.G2Project.add_phase` by moving them from the reader object's ``.Constraints`` 
class variable to the Constraints tree entry's ['Phase'] list (for list items defining constraints) or
the Constraints tree entry's ['_Explain'] dict (for dict items defining constraint help information)

The information in ``.Phase['ISODISTORT']`` is used in :func:`GSASIIconstrGUI.ShowIsoDistortCalc` which shows coordinate and mode
values, the latter with s.u. values. This can be called from the Constraints and Phase/Atoms tree items. 

Before each refinement, constraints are processed as :ref:`described elsewhere <Constraints_processing>`. After a refinement
is complete, :func:`GSASIIstrIO.PrintIndependentVars` shows the shifts and s.u.'s on the refined modes, 
using GSAS-II values, but :func:`GSASIIstrIO.PrintISOmodes` prints the ISODISTORT modes as computed in the web site.


.. _ParameterLimits:

.. index::
   single: Parameter limits

Parameter Limits
==============================

One of the most often requested "enhancements" for GSAS-II would be the inclusion 
of constraints to force parameters such as occupancies or Uiso values to stay within 
expected ranges. While it is possible for users to supply their own restraints that would
perform this by supplying an appropriate expression with the "General" restraints, the 
GSAS-II authors do not feel that use of restraints or constraints are a good solution for
this common problem where parameters refine to non-physical values. This is because when 
this occurs, most likely one of the following cases is occurring: 

#. there is a significant problem 
   with the model, for example for an x-ray fit if an O atom is placed where a S is actually 
   present, the Uiso will refine artificially small or the occupancy much larger than unity 
   to try to compensate for the missing electrons; or
 
#. the data are simply insensitive 
   to the parameter or combination of parameters, for example unless very high-Q data 
   are included, the effects of a occupancy and Uiso value can have compensating effects,
   so an assumption must be made; likewise, with neutron data natural-abundance V atoms 
   are nearly invisible due to weak coherent scattering. No parameters can be fit for a 
   V atom with neutrons.

#. the parameter is non-physical (such as a negative Uiso value) but within 
   two sigma (sigma = standard uncertainty, aka e.s.d.) of a reasonable value, 
   in which case the
   value is not problematic as it is experimentally indistinguishable from an 
   expected value.  

#. there is a systematic problem with the data (experimental error)

In all these cases, this situation needs to be reviewed by a crystallographer to decide 
how to best determine a structural model for these data. An implementation with a constraint
or restraint is likely to simply hide the problem from the user, making it more probable 
that a poor model choice is obtained. 

What GSAS-II does implement is to allow users to specify ranges for parameters 
that works by disabling 
refinement of parameters that refine beyond either a lower limit or an upper limit, where 
either or both may be optionally specified. Parameters limits are specified in the Controls 
tree entry in dicts named as ``Controls['parmMaxDict']`` and ``Controls['parmMinDict']``, where 
the keys are :class:`G2VarObj` objects corresponding to standard GSAS-II variable 
(see :func:`getVarDescr` and :func:`CompileVarDesc`) names, where a 
wildcard ('*') may optionally be used for histogram number or atom number 
(phase number is intentionally not  allowed as a wildcard as it makes little sense 
to group the same parameter together different phases). Note
that :func:`prmLookup` is used to see if a name matches a wildcard. The upper or lower limit
is placed into these dicts as a float value. These values can be edited using the window 
created by the Calculate/"View LS parms" menu command or in scripting with the 
:meth:`GSASIIscriptable.G2Project.set_Controls` function. 
In the GUI, a checkbox labeled "match all histograms/atoms" is used to insert a wildcard
into the appropriate part of the variable name.

When a refinement is conducted, routine :func:`GSASIIstrMain.dropOOBvars` is used to 
find parameters that have refined to values outside their limits. If this occurs, the parameter
is set to the limiting value and the variable name is added to a list of frozen variables 
(as a :class:`G2VarObj` objects) kept in a list in the
``Controls['parmFrozen']`` dict. In a sequential refinement, this is kept separate for 
each histogram as a list in 
``Controls['parmFrozen'][histogram]`` (where the key is the histogram name) or as a list in 
``Controls['parmFrozen']['FrozenList']`` for a non-sequential fit. 
This allows different variables
to be frozen in each section of a sequential fit. 
Frozen parameters are not included in refinements through removal from the 
list of parameters to be refined (``varyList``) in :func:`GSASIIstrMain.Refine` or 
:func:`GSASIIstrMain.SeqRefine`. 
The data window for the Controls tree item shows the number of Frozen variables and
the individual variables can be viewed with the Calculate/"View LS parms" menu window or 
obtained with :meth:`GSASIIscriptable.G2Project.get_Frozen`.
Once a variable is frozen, it will not be refined in any 
future refinements unless the the variable is removed (manually) from the list. This can also 
be done with the Calculate/"View LS parms" menu window or 
:meth:`GSASIIscriptable.G2Project.set_Frozen`.


.. seealso::
  :class:`G2VarObj`
  :func:`getVarDescr` 
  :func:`CompileVarDesc`
  :func:`prmLookup`
  :class:`GSASIIctrlGUI.ShowLSParms`
  :class:`GSASIIctrlGUI.VirtualVarBox`
  :func:`GSASIIstrIO.SaveUsedHistogramsAndPhases`
  :func:`GSASIIstrIO.SaveUpdatedHistogramsAndPhases`
  :func:`GSASIIstrIO.SetSeqResult`
  :func:`GSASIIstrMain.dropOOBvars`
  :meth:`GSASIIscriptable.G2Project.set_Controls`
  :meth:`GSASIIscriptable.G2Project.get_Frozen`
  :meth:`GSASIIscriptable.G2Project.set_Frozen`
