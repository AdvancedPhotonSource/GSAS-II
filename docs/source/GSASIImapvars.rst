*GSASIImapvars: Param Constraints*
=================================================

*Summary/Contents*
----------------------------

Module to implements algebraic contraints, parameter redefinition
and parameter simplification contraints. 

.. contents:: Section Contents 

*Externally-Accessible Routines*
---------------------------------

To define a set of constrained and unconstrained relations, one
defines a list of dictionary defining constraint parameters and their
values, a list of fixed values for each constraint and a list of
parameters to be varied. In addition, one uses
:func:`StoreEquivalence` to define parameters that are equivalent. 
See the :ref:`Constraints Processing section<Constraints_processing>` for details on how 
processing of constraints is done. 

.. tabularcolumns:: |l|p{4in}|

=============================  ===================================================================
  routine                      explanation
=============================  ===================================================================
:func:`InitVars`               This is used to clear out all defined previously 
                               defined constraint information
  
:func:`StoreEquivalence`       Implements parameter redefinition. 
                               This should be called for every set of equivalence relationships.
                               Use :func:`StoreEquivalence` before calling 
                               :func:`GenerateConstraints` 

:func:`ProcessConstraints`     Initially constraints of all types are maintained in lists of 
                               dict entries that are stored in the data tree, 
                               with parameters are stored as 
                               :class:`~GSASIIobj.G2VarObj` objects so that they can 
                               be resolved if the phase/histogram order changes. 
                               :func:`ProcessConstraints` processes this list of dict entries,
                               separating the "Equivalence", "Hold", “Const” and “New Var” 
                               entries for subsequent use.
                               See the :ref:`Constraint Reorganization <ProcessConstraints>`
                               section for more details. 

:func:`EvaluateMultipliers`    Convert any string-specified (formula-based) multipliers to 
                               numbers. Call this before using :func:`GenerateConstraints`. 
                               At present only values in dict for phase (atom/cell) parameters
                               are used to evaluate multipliers containint formulae,
                               but this could be changed if needed. 

:func:`GenerateConstraints`    Generates the internally-used tables from constraints and 
                               equivalences. Checks for internal consistency and repairs 
                               problems where possible. See the 
                               :ref:`Constraint Checking and Grouping <GenerateConstraints>`
                               and :ref:`Equivalence Checking<CheckEquivalences>`
                               sections for more details. 

:func:`Map2Dict`               To determine values for any parameters created in this module,
                               call Map2Dict. This will not apply contraints.

:func:`Dict2Map`               To apply the constraints and equivalences, call this. 
                               It takes values from the new independent parameters and 
                               constraints, and applies them to the parameter dict. 

:func:`Dict2Deriv`             This determines derivatives on independent parameters
                               from those on dependent ones.

:func:`ComputeDepESD`          Use ComputeDepESD to compute uncertainties on dependent variables.

:func:`VarRemapShow`           Use this to show a summary of the parameter remapping.
                               Call after :func:`GenerateConstraints`. 
=============================  ===================================================================

Types of constraints
--------------------

There are four ways to specify constraints, as listed below. 
Note that constraints are initially stored in the 
main section of the GSAS-II data tree under heading ``Constraints``. 
This dict has four keys, 'Hist', 'HAP', 'Global', and 'Phase', 
each containing a list of constraints. An additional set of constraints 
are generated for each phase based on symmetry considerations by calling 
:func:`GSASIIstrIO.GetPhaseData`. 

Note that in the constraints, as stored in the GSAS-II data tree, parameters 
are stored as :class:`GSASIIobj.G2VarObj` objects, as these objects allow for 
changes in numbering of phases, histograms and atoms since :class:`GSASIIobj.G2VarObj` objects 
use random Id's for references.
When constraints are interpreted (in :func:`ProcessConstraints`), 
these references are resolved to the numbered objects by looking up random Id's 
so that the parameter object is converted to a string of form ``ph:hst:VARNAM:at``.

Constraints are initially stored as described in the 
:ref:`constraint definitions <Constraint_definitions_table>`, where the last value in the 
list determines which type of constraint is defined.

Alternate parameters (New Var)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameter redefinition ("New Var" constraints) 
is done by creating an expression that relates several 
parameters: ::

   Mx1 * Px + My1 * Py +... = ::newvar1
   Mx2 * Px + Mz2 * Pz + ... = ::newvar2

where Pj is a GSAS-II parameter name and Mjk is a constant (float) multiplier. 
Alternately, multipliers Mjk can contain a formula (str) that will be evaluated prior 
to the start of the refinement. In a formula, GSAS-II parameters will be replaced by the 
value of the parameter before the formula is evaluated, so ``'np.cos(0::Ax:2)'`` is a valid 
multiplier. At present, only phase (atom/cell) parameters are available for use in 
a formula, from the GUI but this can be expanded if needed. 

This type of constraint describes an alternate 
degree of freedom where parameter Px and Py, etc. are varied to keep 
their ratio 
fixed according the expression. A new variable parameter is assigned to each degree of 
freedom when refined. An example where this can be valuable is when 
two parameters, P1 and P2, have similar values and are highly correlated. It is often better to create a new variable, Ps = P1 + P2, and refine Ps. 
In the later stages of refinement, a second 
variable, Pd = P1 - P2, can be defined and it can be seen if refining Pd is
supported by the data. Another use will be to define parameters that 
express "irrep modes" in terms of the fundamental structural parameters. 

The "New Var" constraints are stored as a type "f"
:ref:`constraint (see definitions)<Constraint_definitions_table>`.

Constrained parameters (Const)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A constraint on a set of variables can be supplied in the form of a
linear algebraic equation: ::

  Nx1 * Px + Ny1 * Py +... = C1

where Cn is a constant (float), where Pj is a GSAS-II parameter name, 
and where Njk is a constant multiplier (float) or a formula (str) that will be evaluated prior 
to the start of the refinement. In a formula, GSAS-II parameters will be replaced by the 
value of the parameter before the formula is evaluated, so ``'np.cos(0::Ax:2)'`` is a valid 
multiplier. At present, only phase (atom/cell) parameters are available for use in 
a formula, but this can be expanded if needed. 

These equations set an interdependence between parameters. 
Common uses of parameter constraints are to set rules that decrease the number of parameters, 
such as restricting the sum of fractional occupancies for atoms that share 
a site to sum to unity, thus reducing the effective number of variables by one. 
Likewise, the Uiso value for a H atom "riding" on a C, N or O atom 
can be related by a fixed offset (the so called B+1 "rule"). 

The "Const" constraints are stored as a type "c" 
:ref:`constraint (see definitions)<Constraint_definitions_table>`.

Equivalenced parameters (Equiv)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simplifed way to set up a constraint equation is to define an equivalence, 
which can be of form: ::

  C1 * P1 = C2 * Py 

or:: 

  C1 * P1 = C2 * P2 = C3 * P3 = ...

where Cn is a constant (float), where Pj is a GSAS-II parameter name. This 
means that parameters Py (or P2 and P3) are determined from (or "slaved" to)
parameter P1. Alternately, equivalences can be created with :func:`StoreEquivalence`
where the multipliers can be a formula (str) that will be evaluated prior to the start of
the refinement. In a formula, GSAS-II parameters will be replaced by the value of the
parameter before the formula is evaluated, so a multiplier can be specifed as ``'2*np.cos(0::Ax:2)'``. 
At present, only phase (atom/cell) parameters are available for use in 
such a formula, but this can be expanded if needed. 

The first parameter (P1 above) 
is considered the independent variable 
and the remaining parameters are dependent variables. The dependent variables 
are then set from the independent variable. 

Note that a constraint expression is conceptually identical to 
defining constraint equations. 
The previous set of equalities could also be written as a set of constraint 
equations in this way: ::

  C1 * P1 - C2 * P2 = 0
  C1 * P1 - C3 * P3 = 0
  ...

In practice, however, 
equivalenced parameters are processed in a 
different and more direct manner than constraint equations. 

A parameter can be used in multiple 
equivalences where it is an independent variable,
but if a parameter were used both as a dependent and independent variable then the order that 
shifts are applied becomes potentially significant. As an example, in this case these two 
equivalences are "chained":: 

  C1 * P1 = C2 * P2
  C2 * P2 = C3 * P3

where P2 is both a dependent and independent variable. Likewise, if a parameter is used both in equivalences 
and in "New Var" or "Const" constraints, it also becomes unclear how this should be processed. It is 
possible to specify equivalences that conflict with constraints. 
Should parameter be used as both a dependent and an independent variable or if a parameter is used both in 
an the equivalence and in a "New Var" or "Const" constraints, the equivalence 
is converted to a constraint (Const) to 
avoid conflicts. The equivalences that require this are addressed in ::func:`GenerateConstraints` where
:func:`CheckEquivalences` is used to locate problematic variables in equivalences 
and then change these equivalences to "Const" equations. Also, unneeded equivalences are removed.

For an example of how equivalences may be used, consider
a material that has **N** O atoms in the asymmetric unit, all in fairly similar bonding environments
and where the diffraction data are sparse. One may wish to reduce the complexity of the model fit to 
these data by defining Uiso for all O atoms to be the same. This is done by selecting Uiso for any one O atom 
as an independent variable in a equivalence and setting the remaining **N-1** other O atom Uiso 
variables as dependent parameters with multipliers of 1. This will require that all O atom Uiso values 
be identical. 
The results of this refinement will be simpler to understand than if a set of
constraint equations is used, because the refined parameter (named as ``ph::Uiso:n``) will be the 
independent variable, corresponding to the first O atom and all other variables would be 
expressed in terms of that variable with a single Equivalence expression. 
The alternate would require **N-1** constraint equations, leaving one degree of freedom with a 
variable would that is likely only indirectly related to the Uiso values.

Equivalenced parameters ("EQUIV" constraints), when defined by users, 
or when created to relate phases, are stored as a type "e" 
:ref:`constraint (see definitions)<Constraint_definitions_table>`.
Symmetry-generated equivalences are generated prior 
to display or refinement in :func:`GSASIIstrIO.GetPhaseData`.
These are not stored in the data tree. 

Hold parameters (Fixed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When parameters are refined where a single refinement flag determines that several variables 
are refined at the same time (examples are: cell parameters, atom positions, anisotropic
displacement parameters, magnetic moments,...) it can be useful to specify that a 
specific parameter should not be varied. These will most commonly be generated due to symmetry, 
but under specific conditions, there may be other good reasons to constrain a parameter. 

The "Hold" constraints are stored as a type "h" 
:ref:`constraint (see definitions)<Constraint_definitions_table>`.

.. _Constraints_processing:

Constraint Processing
---------------------

When constraints will be used or edited, they are processed using a series of
calls. This is done in GSAS-II from several locations:

* For error checking from the tree in :mod:`GSASIIconstrGUI`,
  :func:`GSASIIconstrGUI.CheckConstraints` loads constraints from 
  the data tree.

* When the table of refined parameters is shown, constraints are also 
  processed in function :func:`GSASIIdataGUI.GSASII.OnShowLSParms` using 
  :func:`GSASIIconstrGUI.CheckConstraints`

* To write parameters in the Export sections of the program, 
  :func:`GSASIIfiles.ExportBaseclass.loadParmDict` loads results as well as constraints 
  from the tree. This works a bit differently from the above, so it
  makes direct calls to the constraints routines. 

* For error checking from a GPX file 
  :func:`GSASIIstrIO.ReadCheckConstraints` loads constraints 
  (called in :mod:`GSASIIdataGUI` and :mod:`GSASIIscriptable`), 
  which is similar to :func:`GSASIIconstrGUI.CheckConstraints`. 
  :func:`~GSASIIstrIO.ReadCheckConstraints` is called by 
  :meth:`GSASIIdataGUI.GSASII.OnRefine` and 
  :meth:`GSASIIdataGUI.GSASII.OnSeqRefine` 
  before constraints are generated for use in refinements so they can 
  be shown in the GUI. This is also called to check for errors in
  :class:`GSASIIscriptable.G2Project`. 

* To create the constraints for use in a refinement, in 
  :mod:`GSASIIstrMain`, functions :func:`GSASIIstrMain.Refine` and 
  :func:`GSASIIstrMain.SeqRefine` load and process the constraints again.
  This is repeated here because :func:`~GSASIIstrMain.Refine` and 
  :func:`~GSASIIstrMain.SeqRefine` are intended to operate as stand-alone
  routines that may be called directly.

* After sequential fits have been completed, the previously processed 
  constraint info is read from the sequential results section of the 
  data tree. Function 
  :func:`GSASIIseqGUI.UpdateSeqResults` displays the sequential results
  table also processes constraints. 

TODO: Note that G2stIO.makeTwinFrConstr is called only in one place. It probably needs to be included in all of the above.

When constraints are processed, the following steps are used:  

#. Constraints are stored in separate lists in the data tree to 
   simplify their creation and their GUI display.
   In the initial processing, all of the stored constraints are appended
   into a single list.

#. Then :func:`InitVars` is used to initialize the global variables in 
   this module (:mod:`GSASIImapvars`). This may be done before the previous 
   step.

#. Then :func:`ProcessConstraints` is used to initially process the 
   constraints user-supplied constraints (from the data tree), 
   as described in :ref:`Constraint Reorganization <ProcessConstraints>`.
   When constraints are read from a GPX file, rather than the data tree, use 
   :func:`GSASIIstrIO.ReadConstraints` (which calls :func:`ProcessConstraints`).

#. Symmetry-generated equivalences are then created in 
   :func:`GSASIIstrIO.GetPhaseData`, which also calls 
   :func:`GSASIIstrIO.cellVary` and for Pawley refinements 
   :func:`GSASIIstrIO.GetPawleyConstr`. These are entered directly into this
   module's globals using :func:`StoreEquivalence`.

#. Constraints/equivalences are then checked for possible conflicts with
   :func:`GenerateConstraints` which in turn calls :func:`CheckEquivalences`. 
   These routines group the constraints
   and possibly reorganize them, as discussed below for
   :func:`GenerateConstraints` (:ref:`discussed here <GenerateConstraints>`)
   and for :func:`CheckEquivalences` (:ref:`discussed here <CheckEquivalences>`).
	
#. Note that for debugging, :func:`VarRemapShow` can be called at any point
   after :func:`GenerateConstraints` has been called. This will display the
   generated constraints. 

.. _ProcessConstraints:

Constraint Reorganization (:func:`ProcessConstraints`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`ProcessConstraints` is used to initially process the 
constraints from the list of dict entries. The "Const" and "New Var" are placed into two 
lists (:data:`constrDict` and :data:`fixedList`) that are later used for parameter 
grouping (in :func:`GenerateConstraints`). "Hold" and "Equivalence" constraints are 
separated into separate storage.
 
   For "**Const**" entries, 
     a dict with multiple entries is placed in :data:`constrDict` where 
     each dict key is the parameter name and the value is the multiplier for the parameter, 
     while :data:`fixedList` gets a string value corresponding to the constant value for 
     the expression. 

   For "**New Var**" entries,
     a dict with multiple entries defined identically to 
     to that used in "Const" entries. The differences between "New Var" and "Const" entries is 
     that for "Const" entries, a constant value (as a string) is placed in :data:`fixedList` while
     for "New Var" entries corresponding entry in :data:`fixedList` is None. 
     Also, one or two additional entries are created in the dict for "New Var" constraints:
     an entry with key "_vary" is given the value of True or False
     depending on the refinement flag setting; 
     an entry with key "_name" will be created if the "New Var" parameter has a supplied name.

   For "**Hold**" entries, 
     User-supplied “Hold” constraints are stored in global variable :data:`holdParmList`.
     Initialized in :func:`InitVars`; set in :func:`StoreHold`. Type of hold is stored in
     :data:`holdParmType`.

   **Equivalences** are stored using :func:`StoreEquivalence` into this module's globals 
     (:data:`dependentParmList`, :data:`arrayList`, :data:`invarrayList`, :data:`indParmList`,
     and :data:`symGenList`).
     For each equivalence:

     * a list with one entry, the name of the independent parameter is placed in :data:`indParmList`;
     * a list with one or more parameter name is placed in :data:`dependentParmList`;
     * the value None is added to :data:`arrayList`;
     * a list of multipliers for each dependent variable is placed in :data:`invarrayList`
     * an entry of either True or False is placed in :data:`symGenList`, where True indicates that the entry has been generated from symmetry.

The output from :func:`ProcessConstraints` will have the form as below, 
where the first entry is a "Const" and the second is a "New Var". 

  .. code-block:: python

    constrDict = [
         {'0:12:Scale': 2.0, '0:14:Scale': 4.0, '0:13:Scale': 3.0, '0:0:Scale': 0.5},
         {'2::C(10,6,1)': 1.0, '1::C(10,6,1)': 1.0, '_vary':True}]
    fixedList = ['5.0', None]

.. _GenerateConstraints:

Constraint Checking and Grouping (:func:`GenerateConstraints`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Function :func:`GenerateConstraints` is used to
process the parameter equivalences and constraint lists created in 
:func:`ProcessConstraints` (``constrDict`` and ``fixedList``). :func:`GenerateConstraints` 
is used to generate error/warning messages, to set up lists that are used to show this 
information for the GUI (using :func:`getConstrError`) and to 
generate the information stored in :ref:`global arrays <GlobalVariables>` that are used later to 
apply the constraints. 

When a sequential refinement is in progress, the constraints are scanned for parameters
that have a wildcard (*) for the histogram number, such as 1:*:Scale which would refer 
to the phase fraction for Phase ` in every histogram. The "*" will be replaced with 
the number of the current histogram.

Equivalences are checked with :func:`CheckEquivalences` (described in detail 
:ref:`below <CheckEquivalences>`). This may result in the creation of additional "Hold"
and "Constr" constraints being added to the ``constrDict`` and ``fixedList`` lists.

The "Const" and "New Var" constraint expressions are then scanned for problems:

Constraints cannot be processed without changes if any of the terms within have the following:

* **Undefined parameters** or **Multiplier of zero**

  If any parameters in a constraint are undefined or have a parameter multiplier of zero
  the constraint group is not used.

  If some, but not all, parameters in a constraint are undefined or have a parameter 
  multiplier of zero and remaining valid parameters will be set as "Hold". 
  One exception: atom position constraints (p::dA[xyz]:#) will be assumed as zero. 

* **Hold (Fixed) parameters** and **Unvaried parameters**: New Var constraints

  If any parameters in a new var constraint are either not refined, or are marked as "Hold"
  the constraint can not be varied. Any parameters in that group will be set as "Hold"

* **Hold (Fixed) parameters** and **Unvaried parameters**: Constraint Equations

  If any parameters in a constraint equation are either not refined, or are marked as "Hold"
  those parameters can be removed from the constraint, with an adjustment of the equation 
  sum.  

Constraint expressions ("Const" and "New Var") are sorted by routine :func:`GroupConstraints` into 
groups so that each group contains the minimum number of entries that 
ensures each parameter is referenced in only one group.
This is done by scanning the 
list of dicts in :data:`constrDict` one by one and making a list 
of parameters used in that constraint expression. Any expression that contains
a parameter in that list is added to the current group and those 
parameters are added to this list of parameters. The list of ungrouped 
expressions is then scanned again until no more expressions are added to the 
current group. This process is repeated until every expression has been 
placed in a group. Function :func:`GroupConstraints` returns two lists of lists.
The first has, for each group, a list of the indices in :data:`constrDict` 
that comprise the group (there can be only one). The second list contains, 
for each group, the unique parameter names in that group. 

Each constraint group is then processed. First, wildcard parameters are
renamed (in a sequential refinement). Any held parameters that are used 
in constraints are noted as errors. The number of refined parameters and
the number of parameters that are not defined in the current refinement are 
also noted. It is fine if all parameters in a group are not defined or all are 
not varied, but if some are defined and others not or some are varied and 
others not, this constitutes an error.

The contents of each group is then examined. Groups with a single 
parameter (holds) are ignored. Then for each group, the number 
of parameters in the group (Np) and the number of expressions in the
group (Nc) are counted and for each expression. If Nc > Np, then the constraint 
is overdetermined, which also constitutes an error. 

The parameter multipliers for each expression are then assembled: 

::

   M1a * P1 + M2a * P2 +... Mka * Pk
   M1b * P1 + M2b * P2 +... Mkb * Pk
   ...
   M1j * P1 + M2j * P2 +... Mkj * Pk


From this it becomes possible to create a Nc x Np matrix, which is 
called the constraint matrix:

 .. math::

   \left( \begin{matrix}
   M_{1a}  & M_{2a} &... & M_{ka} \\
   M_{1b}  & M_{2b} &... & M_{kb} \\
   ...  \\
   M_{1j}  & M_{2j}  &... & M_{kj}
   \end{matrix}\right)

When Nc<Np, then additional rows need to be added to the matrix and to 
the vector that contains the value for each row (:data:`fixedList`) where 
values are ``None`` for New Vars and a constant for fixed values. 
This then can describe a system of Np simultaneous equations: 

 .. math::

   \left( \begin{matrix}
   M_{1a}  & M_{2a} &... & M_{ka} \\
   M_{1b}  & M_{2b} &... & M_{kb} \\
   ...  \\
   M_{1j}  & M_{2j}  &... & M_{kj}
   \end{matrix}\right)
   \left( \begin{matrix}
   P_{1} \\
   P_{2} \\
   ...  \\
   P_{k}
   \end{matrix}\right)
   = 
   \left( \begin{matrix}
   C_{1} & C_{2} &  ... & C_{k}
   \end{matrix}\right)

This is done by creating a square matrix from the group using 
:func:`_FillArray`. The top Nc rows in the matrix are filled
as described above. Then :func:`_RowEchelon` is used to see if
those entries in the matrix can be coverted to row-echelon form. This 
will raise an Exception there is linear dependence between the initial Nc rows 
(which means that no matter what values are used for any remaining rows, that the matrix 
will be singular). If that is not the case and Nc<Np then any remaining rows that
were not specified are filled in. For each of these rows, first only the 
diagonal element in that row of the matrix is set to 1 
and the upper portion of the matrix is again tested with :func:`_RowEchelon` 
to check for linear independence. This is likely to be non-singular, 
but should :func:`_RowEchelon` fail, 
:func:`_FillArray` will then try setting each other element in that row to either 
1 or -1. One of those options should be linearly independent from every other 
row of the matrix. 

The  
`Gram-Schmidt process <http://en.wikipedia.org/wiki/Gram-Schmidt>`_, 
implemented  in :func:`GramSchmidtOrtho`, is used to find orthonormal unit 
vectors which are used to replace the remaining Np-Nc rows of the matrix. This will fail with 
a ``ConstraintException`` if this is not possible (singular matrix), but that would be 
unexpected since the matrix could be converted to row-echelon form. The 
Gram-Schmidt result is placed in :data:`constrArr` as a numpy array. 

Rows in the matrix corresponding to "New Var" constraints and those that 
were generated by the Gram-Schmidt process are provided with parameter names.
These names are generated using :data:`paramPrefix`, which is set to ``"::constr"``, 
plus a number to make the new parameter name unique, 
unless a name was specified for the 
"New Var" entry by using a ``"_name"`` element in the constraint dict.

Finally the parameters used as input to the constraint are placed in 
this module's globals
:data:`dependentParmList` and the constraint matrix is 
placed in into  :data:`arrayList`. This can be used to compute 
the initial values for "New Var" parameters. The inverse of the 
constraint matrix is placed in :data:`invarrayList` and a list 
of the "New Var" parameters and a list of the fixed values (as str's) 
is placed in :data:`indParmList`. 
Finally the appropriate entry in :data:`symGenList` is set to 
False to indicate that this is not a symmetry generated constraint. 

.. _CheckEquivalences:

Equivalence Checking and Reorganization (:func:`CheckEquivalences`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equivalences need to be checked for usages that could be internally conflicted 
or have possible conflicts with other constraints. Function :func:`CheckEquivalences`
is called within :func:`GenerateConstraints` to diagnose and where
possible resolve such uses, as discussed below. 

**Mixed parameter use:**

 **Note** that multiple passes, cycling through the equivalences may
 be needed to find all mixed-use parameters, as will be discussed
 further, below.

* A parameter should not show up as a dependent variable in two equivalence expressions, 
  such as:: 

      ::x1 -> ::x3
      ::x2 -> ::x3

  This will be processed by turning the equivalences into two constraint equations::

      ::x1 - ::x3 = 0
      ::x2 - ::x3 = 0

  which can be satisfied when ``::x1 = ::x2 = ::x3``. If  ``::x1`` and ``::x2`` had been 
  intended to be independent parameters, then the above equivalences would be conflict and
  cannot be statisfied.

* If a parameter is used both as an independent and as a dependent variable (*chaining*), 
  as is in these two equivalence expressions::

      ::x1 -> ::x2 & ::x4
      ::x2 -> ::x3

  This can also be addressed by turning these equivalences into three constraint equations::

   ::x1 - ::x2 = 0
   ::x1 - ::x4 = 0
   ::x2 - ::x3 = 0

  which can be satisfied when ``::x1 = ::x2 = ::x3 = ::x4``

*  Use of parameters in both equivalences and "Const" or "New Var" constraint expressions makes
   logical sense::

   ::x1 -> ::x2 & ::x4
   ::x2 + ::x3 = 0

   This can also be addressed by turning the equivalence into two constraint equations::

   ::x1 - ::x2 = 0
   ::x1 - ::x4 = 0

   With the addition of the "Const" equation (``::x2 + ::x3 = 0``), the solution will require 
   ``::x1 = ::x2 = -1.0*::x3 = ::x4``

* Cycling is needed to find all equivalences that must be converted.
  Consider this set of constraints::

   ::x2 + ::x3 = 0
   ::x1 -> ::x2
   ::x1 -> ::x4

 In the first pass the equivalence with ``::x2`` would be converted to a "Const" constraint 
 and in the second pass
 the other equivalence with ``::x1`` would be converted. 


**Mixing Hold (Fixed) parameters in equivalences**

* If one parameter (or more) is designated as a "Hold" in an equivalence, then all parameters in that 
  equivalence cannot be varied. Considering this equivalence::

      ::x1 -> ::x2 & ::x4

  If any of the three parameters (``::x1``, ``::x2``, or `::x4`) are marked as Hold, then 
  the other two parameters may not be varied and will also be set with a "Hold".


**Unvaried parameters in equivalences**

* If no parameters in an equivalence are varied, then the equivalence is ignored.

* If only some parameters are marked as varied then 
  *none of the parameters can be varied*; any varied parameters will be set with a "Hold".


**Undefined parameters in equivalences**

  Parameters may be placed in equivalences that are not actually defined in a project. 
  This can occur in two ways. If an equivalence is created in the GUI for a parameter that
  is later supplanted with a different model (for example, changing from isotropic size 
  broadening to uniaxial broadening replaces the isotropic broadening term with two different
  uniaxial terms) or symmetry may require restrictions on anisotropic ADPs that are not 
  in use).

* If the independent parameter is undefined, then any dependent parameters that are defined 
  are set as "Hold" and the equivalence is ignored. 

* If all dependent parameters are undefined, then the equivalence is ignored. 

* If a dependent parameter is undefined, then that parameter is dropped from the equivalence. 


**Multiplier of zero in equivalences**

  Any dependent parameter that has a multiplier of zero will be dropped from the equivalence.
  If no terms remain, then the equivalence is ignored. (Independent parameters do not 
  have a multiplier). 


.. _GlobalVariables:

*Global Variables*
------------------

This module uses a number of global variables. One set is used to store the 
constraints and equivalences after processing by :func:`StoreEquivalence` and 
:func:`GenerateConstraints`.  
These globals are expected to be used only by this module's (:mod:`GSASIImapvars`) internal routines.

Lists with information from Constraint Equation and New Var constraints. Each entry
in these variables is related to a group of constraints. 

.. tabularcolumns:: |l|p{4.5in}|

=============================  ===================================================================
  variable                      explanation
=============================  ===================================================================
:data:`dependentParmList`        a list containing group of lists of
                                 parameters used in the group. 
                                 The columns of the matrices in :data:`arrayList` match 
                                 the order of parameters here.
                                 Note that parameters listed in
                                 dependentParmList will not be included in the Hessian as their
                                 derivatives will not affect the model

:data:`indParmList`              a list containing groups of variables or constants matching 
                                 the columns of the matrices in :data:`invarrayList`. 

:data:`arrayList`                a list containing group of relationship matrices to relate
                                 parameters in dependentParmList to those in indParmList. 

:data:`invarrayList`             a list containing group of relationship matrices to relate
                                 parameters in indParmList to those in dependentParmList. 
                                 Unlikely to be used externally.

:data:`symGenList`               a list of boolean values that will be True to indicate 
                                 that an equivalence was generated internally GSAS-II 
                                 meaning it is generated based on symmetry, twining 
                                 or Pawley overlap.

=============================  ===================================================================

Lists with information from Hold and Equivalence constraints. Each entry
in these variables is related to a group of constraints. 

.. tabularcolumns:: |l|p{4.5in}|

=============================  ===================================================================
  variable                       explanation
=============================  ===================================================================
:data:`holdParmList`             a list of parameters that have been marked as "Hold". 
                                 Unlikely to be accessed outside this module.
                                 Initialized in :func:`InitVars`; set in :func:`StoreHold`.

:data:`holdParmType`             The reason why a parameter has been marked as "Hold". 
                                 Unlikely to be accessed outside this module.
                                 Initialized in :func:`InitVars`; set in :func:`StoreHold`.

:data:`constrParms`              dict with lists of variables in equivalences, 
                                 constraint equations and new var expressions.
                                 Used within :func:`GetIndependentVars`,
                                 and :func:`GetDependentVars`. 
                                 Best if not referenced outside this module. 
                                 Contains elements:

                                 * 'dep-equiv': dependent parameters set by equivalences
                                 * 'dep-constr': dependent parameters set by 
                                   constraint equations or new var expressions
                                 * 'indep-equiv': dependent parameters used in equivalences
                                 * 'indep-constr': dependent parameters created from 
                                   constraint equations or new var expressions

:data:`saveVaryList`             a list of the varied parameters used when constraints 
                                 were last processed.
=============================  ===================================================================


A second set of global variables are set in :func:`GenerateConstraints` with lists of parameter
names from equivalences and constraints. Used in :func:`CheckEquivalences` and 
:func:`getConstrError`.

.. tabularcolumns:: |l|p{4.5in}|

=============================  ===================================================================
  variable                      explanation
=============================  ===================================================================
:data:`depVarList`              a list of the parameters used in equivalences as dependent 
                                parameters for all equivalences initially specified (including 
                                those to be reclassified as "Constr" constraints.)

:data:`indepVarList`            a list of the parameters used in equivalences as independent 
                                parameters for all equivalences initially specified (including 
                                those to be reclassified as "Constr" constraints.)

:data:`constrVarList`           a list of the parameters that are used in "Constr" or 
                                "New Var" constraints. Does not include those in equivalences
                                to be reclassified as "Constr" constraints.)
=============================  ===================================================================


A third set of global variables to store equivalence warning information. 
Set in :func:`CheckEquivalences` and :func:`GenerateConstraints`.
Used in :func:`getConstrError` to display warning messages.

.. tabularcolumns:: |l|p{4.5in}|

=============================  ===================================================================
  variable                      explanation
=============================  ===================================================================
:data:`convVarList`             parameters in equivalences that were converted to "Const"
                                constraints

:data:`multdepVarList`          parameters used as dependent parameters in equivalences 
                                multiple times

:data:`unvariedParmsList`       parameters used in equivalences and constraints 
                                that are not varied

:data:`undefinedVars`           parameters used in equivalences 
                                that are not defined in the parameter dict (parmDict)

:data:`groupErrors`             parameters in constraints that cause grouping errors
=============================  ===================================================================



*GSASIImapvars Routines/variables*
---------------------------------------

.. automodule:: GSASIImapvars
    :members: 
    :private-members:
    :special-members:
