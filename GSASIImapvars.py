# TODO: revisit SeqRefine and :meth:`GSASIIdataGUI.GSASII.OnSeqRefine` and :func:`GSASIIseqGUI.UpdateSeqResults`

# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
"""
*GSASIImapvars: Parameter constraints*
======================================

Module to implements algebraic contraints, parameter redefinition
and parameter simplification contraints. 


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
  :func:`GSASIIIO.loadParmDict` loads results as well as constraints 
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
   :func:`GenerateConstraints`, this requires grouping the constraints, 
   as described below.

#. :func:`GenerateConstraints` is then called to 
   create the constraints that will be used, 
   :ref:`see below <GenerateConstraints>` for more details. 

#. For debugging constraints, :func:`VarRemapShow` can be called after 
   :func:`GenerateConstraints` to display the generated constraints. 

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

     * a list with one entry, the name of the independent parameter is placed in :data:`~GSASIImapvars.indParmList`;
     * a list with one or more parameter name is placed in :data:`~GSASIImapvars.dependentParmList`;
     * the value None is added to :data:`~GSASIImapvars.arrayList`;
     * a list of multipliers for each dependent variable is placed in :data:`~GSASIImapvars.invarrayList`
     * an entry of either True or False is placed in :data:`~GSASIImapvars.symGenList`, where True indicates that the entry has been generated from symmetry.

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

   \\left( \\begin{matrix}
   M_{1a}  & M_{2a} &... & M_{ka} \\\\
   M_{1b}  & M_{2b} &... & M_{kb} \\\\
   ...  \\\\
   M_{1j}  & M_{2j}  &... & M_{kj}
   \\end{matrix}\\right)

When Nc<Np, then additional rows need to be added to the matrix and to 
the vector that contains the value for each row (:data:`fixedList`) where 
values are ``None`` for New Vars and a constant for fixed values. 
This then can describe a system of Np simultaneous equations: 

 .. math::

   \\left( \\begin{matrix}
   M_{1a}  & M_{2a} &... & M_{ka} \\\\
   M_{1b}  & M_{2b} &... & M_{kb} \\\\
   ...  \\\\
   M_{1j}  & M_{2j}  &... & M_{kj}
   \\end{matrix}\\right)
   \\left( \\begin{matrix}
   P_{1} \\\\
   P_{2} \\\\
   ...  \\\\
   P_{k}
   \\end{matrix}\\right)
   = 
   \\left( \\begin{matrix}
   C_{1} & C_{2} &  ... & C_{k}
   \\end{matrix}\\right)

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
:data:`~GSASIImapvars.dependentParmList` and the constraint matrix is 
placed in into  :data:`~GSASIImapvars.arrayList`. This can be used to compute 
the initial values for "New Var" parameters. The inverse of the 
constraint matrix is placed in :data:`~GSASIImapvars.invarrayList` and a list 
of the "New Var" parameters and a list of the fixed values (as str's) 
is placed in :data:`~GSASIImapvars.indParmList`. 
Finally the appropriate entry in :data:`~GSASIImapvars.symGenList` is set to 
False to indicate that this is not a symmetry generated constraint. 

.. _CheckEquivalences:

Equivalence Checking and Reorganization (:func:`CheckEquivalences`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equivalences need to be checked for usages that could be internally conflicted 
or have possible conflicts with other constraints. 

**Mixed parameter use:**

 **Note** that cycling through the equivalences may be needed to find all mixed-use 
 parameters, see below.

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

:data:`dependentVars`            a list of dependent variables in equivalences, compiled 
                                 from (:data:`dependentParmList`).
                                 Used within :func:`GetDependentVars`.

:data:`independentVars`          a list of dependent variables in equivalences.
                                 Used within :func:`GetIndependentVars`.

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

Note that parameter names in GSAS-II are strings of form ``<ph#>:<hst#>:<nam>`` or ``<ph#>::<nam>:<at#>``
where ``<ph#>`` is a phase number, ``<hst#>`` is a histogram number and ``<at#>`` is an atom number. 
``<nam>`` is a name that determines the parameter type (see :func:`GSASIIobj.CompileVarDesc`). When 
stored in the data tree, parameters are saved as :class:`GSASIIobj.G2VarObj` objects 
so that they can be resolved if the phase/histogram order changes. 

"""

from __future__ import division, print_function
import copy
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIobj as G2obj 
# data used for constraints; 
debug = False # turns on printing as constraint input is processed

#------------------------------------------------------------------------------------
# Global vars used for storing constraints and equivalences after processing
#   note that new var and constraint equations are stored together in groups,
#   where each constraint group contains those parameters that must be handled 
#   together. Equivalences are also stored in these 

dependentParmList = []
'''a list of lists where each item contains a list of parameters in each constraint group. 
note that parameters listed in dependentParmList should not be refined directly.'''
indParmList = [] # a list of names for the new parameters
'''a list of lists where each item contains a list for each constraint group with 
fixed values for constraint equations and names of generated/New Var parameters.
In the case of equivalences, the name of a single independent parameter is stored.
'''
arrayList = []
'''a list of of relationship matrices that map model parameters in each 
constraint group (in :data:`dependentParmList`) to 
generated (New Var) parameters.
'''
invarrayList = []
'''a list of of inverse-relationship matrices that map constrained values and 
generated (New Var) parameters (in :data:`indParmList`) to model parameters
(in :data:`dependentParmList`). 
'''
symGenList = []
'''A list of flags that if True indicates a constraint was generated by symmetry
'''
holdParmList = []
'''List of parameters that should not be refined ("Hold"s).
Set in :func:`StoreHold`. Initialized in :func:`InitVars`.
'''
holdParmType = {}
'''The reason why a parameter has been marked as "Hold". 
Initialized in :func:`InitVars`; set in :func:`StoreHold`.
'''
dependentVars = []
'''A list of dependent variables in equivalences, compiled from (:data:`dependentParmList`).
Used within :func:`GetDependentVars`
'''
independentVars = []
'''A list of dependent variables in equivalences, compiled from (:data:`indParmList`).
Used within :func:`GetIndependentVars`
'''
saveVaryList = []
'''A list of the varied parameters that was last supplied when constraints were
processed. This is set in :func:`GenerateConstraints` and updated in 
:func:`Map2Dict`. Used in :func:`VarRemapShow`
'''
#------------------------------------------------------------------------------------
# Global vars set in :func:`GenerateConstraints`. Used for intermediate processing of
# constraints.

constrVarList = []
'List of parameters used in "Constr" and "New Var" constraints'
indepVarList = []
'A list of all independent parameters in equivalences'
depVarList = []
'A list of all dependent parameters in equivalences'

#------------------------------------------------------------------------------------
# global variables used in :func:`getConstrError` to store error and warning information.
# set in CheckEquivalences and in GenerateConstraints
convVarList = []
'parameters in equivalences that were converted to "Const" constraints'
multdepVarList = []
'parameters used as dependents multiple times in equivalences'
unvariedParmsList = []
'parameters used in equivalences that are not varied'
undefinedVars = []
'parameters used in equivalences that are not defined in the parameter dict'
groupErrors = []
'parameters in constraints where parameter grouping and matrix inversion fails'

#------------------------------------------------------------------------------------
paramPrefix = "::constr"
'A prefix for generated parameter names'
consNum = 0
'The number to be assigned to the next constraint to be created'

#------------------------------------------------------------------------------------
class ConstraintException(Exception):
    '''Defines an Exception that is used when an exception is raised processing constraints.
    Raised in :func:`GenerateConstraints` during sequential fits. Possible (but highly unlikely) 
    to be raised in :func:`CheckEquivalences` (called by :func:`GenerateConstraints`) if an 
    infinite loop is detected.
    Also raised in :func:`GramSchmidtOrtho` and :func:`_SwapColumns` but caught 
    within :func:`GenerateConstraints`.
    '''
    pass

def InitVars():
    '''Initializes all constraint information'''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum,symGenList
    dependentParmList = [] # contains a list of parameters in each group
    arrayList = [] # a list of of relationship matrices 
    invarrayList = [] # a list of inverse relationship matrices 
    indParmList = [] # a list of names for the new parameters
    consNum = 0 # number of the next constraint to be created
    symGenList = [] # Flag if constraint is generated by symmetry
    global holdParmList,holdParmType
    holdParmList = []
    holdParmType = {}

def VarKeys(constr):
    """Finds the keys in a constraint that represent parameters
    e.g. eliminates any that start with '_'

    :param dict constr: a single constraint entry of form::

        {'var1': mult1, 'var2': mult2,... '_notVar': val,...}

        (see :func:`GroupConstraints`)
    :returns: a list of keys where any keys beginning with '_' are
      removed.
    """
    return [i for i in constr.keys() if not i.startswith('_')]


def GroupConstraints(constrDict):
    """Divide the constraints into groups that share no parameters.

    :param dict constrDict: a list of dicts defining relationships/constraints

    ::
    
       constrDict = [{<constr1>}, {<constr2>}, ...]

    where {<constr1>} is {'var1': mult1, 'var2': mult2,... }

    :returns: two lists of lists:
    
      * a list of grouped contraints where each constraint grouped containts a list
        of indices for constraint constrDict entries
      * a list containing lists of parameter names contained in each group
      
    """
    assignedlist = [] # relationships that have been used
    groups = [] # contains a list of grouplists
    ParmList = []
    for i,constrI in enumerate(constrDict):
        if i in assignedlist: continue # already in a group, skip
        # starting a new group
        grouplist = [i,]
        assignedlist.append(i)
        groupset = set(VarKeys(constrI))
        changes = True # always loop at least once
        while(changes): # loop until we can't find anything to add to the current group
            changes = False # but don't loop again unless we find something
            for j,constrJ in enumerate(constrDict):
                if j in assignedlist: continue # already in a group, skip
                if len(set(VarKeys(constrJ)) & groupset) > 0: # true if this needs to be added
                    changes = True
                    grouplist.append(j)
                    assignedlist.append(j)
                    groupset = groupset | set(VarKeys(constrJ))
        group = sorted(grouplist)
        varlist = sorted(list(groupset))
        groups.append(group)
        ParmList.append(varlist)
    return groups,ParmList

def GenerateConstraints(varyList,constrDict,fixedList,parmDict=None,SeqHist=None):
    '''Takes a list of relationship entries that have been stored by 
    :func:`ProcessConstraints` into lists ``constrDict`` and ``fixedList``

    This routine then calls :func:`CheckEquivalences` for internal 
    consistency. This includes converting equivalenced variables into 
    constraints when a variable is used in both. 

    Once checked, parameters are grouped so that any parameters that are used in 
    more than one constraint are grouped together. This allows checking for incompatible
    logic (for example, when four constraints are specified for three variables). 

    If parmDict is not None, the parameter groups are checked for constraints where 
    some parameters are varied, but not others. If so, the value for that unvaried 
    parameter is subtracted from the constant in the constraint. 

    Once all checks are complete, the constraints are then
    converted to the form used to apply them, saving them as global variables within 
    this module. 

    :param list varyList: a list of parameters names (strings of form
      ``<ph>:<hst>:<nam>``) that will be varied. Note that this is changed here. 
    
    :param dict constrDict: a list of dicts defining relationships/constraints
      (as described in :func:`GroupConstraints`)

    :param list fixedList: a list of values specifying a fixed value for each
      dict in constrDict. Values are either strings that can be converted to
      floats, float values or None if the constraint defines a new parameter.
      
    :param dict parmDict: a dict containing all parameters defined in current
      refinement.

    :param int SeqHist: number of current histogram, when used in a sequential
      refinement. None (default) otherwise. Wildcard parameter names are
      set to the current histogram, when found if not None.

    :returns: errmsg,warning,groups,parmlist

      **errmsg**
        Is an error message or empty if no errors were found
      **warning**
        Is a warning message about constraints that have been ignored or changed
      **groups**
        Lists parameter groups
      **parmlist**
        Lists parameters in each parameter groups
    '''
    warninfo = {'msg':'', 'shown':{}}
    def warn(msg,cdict=None,val=None):
        if cdict is not None and cdict != warninfo['shown']:
            warninfo['shown'] = cdict
            if warninfo['msg']: warninfo['msg'] += '\n'
            if '_vary' in cdict:
                warninfo['msg'] += '\nProblem with new var expression: ' + _FormatConstraint(cdict,cdict.get('_name','New Var'))
            else:
                warninfo['msg'] += '\nProblem with constraint equation: ' + _FormatConstraint(cdict,val)
        if warninfo['msg']: warninfo['msg'] += '\n'
        warninfo['msg'] += '  ' + msg
    
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    # lists of parameters used for error reporting
    global undefinedVars # parameters that are used in equivalences but are not defined
    undefinedVars = []
    global groupErrors # parameters in constraints that cause grouping errors
    groupErrors = []
    global saveVaryList
    saveVaryList = copy.copy(varyList)
    
    global dependentVars # List of dependent and independent variables for all constraints
    global independentVars

    errmsg = ''  # save error messages here. If non-blank, constraints cannot be used.
    warning = '' # save informational text messages here.

    # Process the equivalences; If there are conflicting parameters, move them into constraints
    warning = CheckEquivalences(constrDict,varyList,fixedList,parmDict)
    
    # find parameters used in constraint equations & new var assignments (all are dependent)
    global constrVarList 
    constrVarList = []
    for cdict in constrDict:
        constrVarList += [i for i in cdict if i not in constrVarList and not i.startswith('_')]

    # look through "Constr" and "New Var" constraints looking for zero multipliers and
    # Hold, Unvaried & Undefined parameters
    skipList = []
    invalidParms = []
    for cnum,(cdict,fixVal) in enumerate(zip(constrDict,fixedList)):
        valid = 0          # count of good parameters
        # error reporting
        zeroList = []      # parameters with zero multipliers
        holdList = []      # parameters with "Hold"'s
        noVaryList = []    # parameters not varied
        noWildcardList = [] # wildcard parameters in non-sequential fit
        notDefList = []    # parameters not defined
        # processing to be done
        problem = False    # constraint must be dropped 
        dropList = []      # parameters to remove
        for var in VarKeys(cdict):   # assemble warning info
            if cdict[var] == 0:    # zero multiplier
                if var not in zeroList: zeroList.append(var)
                if var not in dropList: dropList.append(var)
            elif var in holdParmList:   # hold invalid in New Var, drop from constraint eqn
                holdList.append(var)
                if fixVal is None:
                    problem = True
                else:
                    if var not in dropList: dropList.append(var)
            elif ':*:' in var :  # wildcard still present should be treated as undefined
                if var not in undefinedVars: undefinedVars.append(var)
                noWildcardList.append(var)
                problem = True
            elif parmDict is not None and var not in parmDict: # not defined, constraint will not be used
                if var not in undefinedVars: undefinedVars.append(var)
                notDefList.append(var)
                if ':dAx:' in var or ':dAy:' in var or ':dAz:' in var: # coordinates from undefined atoms 
                    if fixVal is None:
                        problem = True  # invalid in New Var
                    else:
                        if var not in dropList: dropList.append(var) # ignore in constraint eqn
                else:
                    problem = True
            elif var not in varyList and fixVal is not None:  # unvaried, constraint eq. only
                if var not in unvariedParmsList: unvariedParmsList.append(var)
                noVaryList.append(var)
                dropList.append(var)
            else:
                valid += 1
        for l,m in ((zeroList,"have zero multipliers"), # show warning
                      (holdList,'set as "Hold"'),
                      (noVaryList,"not varied"),
                      (noWildcardList,"wildcard in non-sequential fit"),
                      (notDefList,"not defined")):
            if l and cdict.get('_vary',True): # true for constr eq & varied New Var
                msg = "parameter(s) " + m + ': '
                for i,v in enumerate(l):
                    if i != 0: msg += ', '
                    msg += v
                warn(msg,cdict,fixVal)
        if valid == 0: # no valid entries
            warn('Ignoring constraint, no valid entries',cdict)
            skipList.append(cnum)
        elif problem: # mix of valid & refined and undefined items, cannot use this
            if cdict.get('_vary',True): # true for constr eq & varied New Var
                warn('New Var constraint will be ignored',cdict)
                skipList.append(cnum)
            invalidParms += VarKeys(cdict)
        elif len(dropList) > 0: # mix of valid and problematic items, drop problem vars, but keep rest
            if GSASIIpath.GetConfigValue('debug'): 
                msg = ''
                for v in dropList:
                    if msg: msg += ' ,'
                    msg += v
                warn('removing: '+msg,cdict)
            value = fixedList[cnum]
            for var in dropList:   # do cleanup
                # NB expressions in constraint multipliers have already been evaluated
                if ':dAx:' in var or ':dAy:' in var or ':dAz:' in var: 
                    pass # treat delta coords as 0; no change in total is needed
                elif cdict[var] != 0: 
                    value = float(value) - cdict[var]*parmDict[var]
                del cdict[var]
            if float(value) != float(fixedList[cnum]): fixedList[cnum] = float(np.round(value,12))
            if GSASIIpath.GetConfigValue('debug'):
                warn('revised as: '+_FormatConstraint(constrDict[cnum],fixedList[cnum]))
    for i in list(range(len(constrDict)-1,-1,-1)): # remove the dropped constraints
        if i in skipList:
            del constrDict[i]
            del fixedList[i]
            
    for i in invalidParms: StoreHold(i,"Used in invalid constraint")
    if warning: warning += '\n'
    warning += warninfo['msg']
            
    groups,parmlist = GroupConstraints(constrDict)

    # now process each group and create the relations that are needed to form
    # a non-singular square matrix
    # Now check that all parameters are varied (probably do not need to do this
    # any more). For constraint equations, if all are varied, set VaryFree to True
    # and all newly created relationships will be varied. For NewVar constraints,
    # vary if the vary flag was set. 
    for group,depPrmList in zip(groups,parmlist):
        if len(depPrmList) < len(group): # too many relationships -- no can do
            if errmsg: errmsg += '\n'
            errmsg += "Over-constrained input. "
            errmsg += "There are more constraints (" + str(len(group))
            errmsg += ") than parameters (" + str(len(depPrmList)) + ")\nin these constraints:"
            for rel in group:
                errmsg += '\n\t'+ _FormatConstraint(constrDict[rel],fixedList[rel])
            groupErrors += depPrmList
            continue # go on to next group

        try:
            constrArr = _FillArray(group,constrDict,depPrmList)
        except Exception as err:
            if errmsg: errmsg += '\n'
            if 'Initial' in str(err): 
                errmsg += "\nSingular input. "
                errmsg += "There are internal inconsistencies in these constraints:"
            else:
                errmsg += "\nError expanding matrix with these constraints:"
                for rel in group:
                    errmsg += '\n\t' + _FormatConstraint(constrDict[rel],fixedList[rel])
            groupErrors += depPrmList
            continue

        try:
            GramSchmidtOrtho(constrArr,len(group))
        except:
            if errmsg: errmsg += '\n'
            errmsg += "\nUnexpected singularity with constraints group (in Gram-Schmidt)"
            for rel in group:
                errmsg += '\n\t' + _FormatConstraint(constrDict[rel],fixedList[rel])
            groupErrors += depPrmList
            continue
        
        try: 
            invConstrArr = np.linalg.inv(constrArr)
        except:
            if errmsg: errmsg += '\n'
            errmsg += "\nSingular input. "
            errmsg += "The following constraints are not "
            errmsg += "linearly independent\nor do not "
            errmsg += "allow for generation of a non-singular set.\n"
            errmsg += 'This is unexpected. Please report this (toby@anl.gov)'
            for rel in group:
                errmsg += '\n\t' + _FormatConstraint(constrDict[rel],fixedList[rel])
            groupErrors += depPrmList
            continue

        # scan through current group looking for new var assignments
        hasNewVar = False
        for i,rel in enumerate(group): 
            if fixedList[rel] is None: 
                hasNewVar = True # there a New Var relationship in this group
                break
        else: # only constraint equations, check for unvaried parameters in one
            unvaried = False
            # this should not happen as they should have been removed
            for var in depPrmList:
                if var not in varyList: 
                    unvaried = True
                    break
            if unvaried: # something is not varied: skip group & remove all parameters from varyList
                for var in depPrmList: StoreHold(var,'Unexpected: mixed use')
                if GSASIIpath.GetConfigValue('debug'): 
                    print('Unexpected: Constraint group ignored (some parameters unvaried)')
                    for rel in group:
                        print ('   '+_FormatConstraint(constrDict[rel],fixedList[rel]))                    
                continue
        maplist = []   # value or param name mapped to each row in expression matrix            
        if not hasNewVar: # constraint equations; all remaining entries varied, vary generated 
            for i in range(len(depPrmList)):
                if i >= len(group): # tag generated degrees of freedom with names and vary them
                    varname = paramPrefix + str(consNum) # assign a unique name
                    consNum += 1
                    maplist.append(varname)
                    varyList.append(varname)
                else:
                    maplist.append(fixedList[rel])
        else:   # ------------------------- groups with new var assignments, vary only NV's w/flags set
            for i,rel in enumerate(group):
                if fixedList[rel] is None:
                    varname = constrDict[rel].get('_name','::?????')
                    maplist.append(varname)
                    if  constrDict[rel].get('_vary',False):
                        varyList.append(varname)
                else:
                   maplist.append(fixedList[rel])   # constraint equation
            for i in range(len(depPrmList)):
                if i >= len(group): # name generated degrees of freedom
                    varname = paramPrefix + str(consNum) # assign a unique name
                    consNum += 1
                    maplist.append(varname)
            for var in depPrmList: StoreHold(var,'New Var use')
        # keep this group
        dependentParmList.append(depPrmList)
        arrayList.append(constrArr)
        invarrayList.append(invConstrArr)
        indParmList.append(maplist)
        symGenList.append(False)
            
    if errmsg and SeqHist is not None:
        print (' *** ERROR in constraint definitions! ***')
        print (errmsg)
        if warning:
            print (' also note warnings in constraint processing:')
            print (warning)
        raise ConstraintException
    elif errmsg:
        return errmsg,warning,None,None

    # Make list of dependent and independent variables for all constraints
    dependentVars = []
    independentVars = []
    for varlist,mapvars in zip(dependentParmList,indParmList):  # process all constraints
        for mv in mapvars:
            if type(mv) is float: continue
            if mv not in independentVars: independentVars.append(mv)
        for mv in varlist:
            if mv not in dependentVars: dependentVars.append(mv)
            StoreHold(mv,'dependent param')
    saveVaryList = copy.copy(varyList)  # save varyList so it can be used within module

    # if equivMoved:
    #     print(60*'=')
    #     print('Constraints were reclassified to avoid conflicts, as below:')
    #     print(mvMsg)
    #     print('New constraints are:')
    #     print (VarRemapShow(varyList,True))
    #     print(60*'=')
    return errmsg,warning,groups,parmlist # saved for sequential fits
    
def CheckEquivalences(constrDict,varyList,fixedList,parmDict=None):
    '''Process equivalence constraints, looking for conflicts such as 
    where a parameter is used in both an equivalence and a constraint expression
    or where chaining is done (A->B and B->C). 

    Removes equivalences or parameters from equivalences or converts equivalences to
    constraints as described for :ref:`Equivalence Checking and Reorganization <CheckEquivalences>`.

    :param dict constrDict: a list of dicts defining relationships/constraints
    :param list varyList: list of varied parameters (defined during refinements only)
    :param list fixedList: a list of values specifying a fixed value for each
       dict in constrDict. Values are either strings that can be converted to
       floats or ``None`` if the constraint defines a new parameter rather
       than a constant.
    :param dict parmDict: a dict containing defined parameters and their values. Used to find 
       equivalences where a parameter is has been removed from a refinement. 

    :returns: warning messages about changes that need to be made to equivalences 
    '''
    
    warninfo = {'msg':'', 'shown':-1}
    def warnEqv(msg,cnum=None):
        if cnum is not None and cnum != warninfo['shown']:
            warninfo['shown'] = cnum
            if warninfo['msg']: warninfo['msg'] += '\n'
            warninfo['msg'] += '\nProblem with equivalence: ' + _showEquiv(
                dependentParmList[cnum],indParmList[cnum],invarrayList[cnum])
        if warninfo['msg']: warninfo['msg'] += '\n'
        warninfo['msg'] += '  ' + msg

    global depVarList # parameters used in equivalences as dependent parameters
    global indepVarList # parameters used in equivalences as independent parameters
    global constrVarList  # parameters used in other constraints
             
    # lists of parameters used for error reporting
    global undefinedVars # parameters that are used in equivalences but are not defined
    global convVarList # parameters in equivalences that will be converted to constraints
    convVarList = [] # parameters in equivalences to be made into constraints
    global multdepVarList
    multdepVarList = [] # list of dependent parameters used in more than one equivalence

    # local vars
    dropVarList = []  # parameters that can be removed from equivalences 
    removeList = []  # equivalences that are not needed
    convertList = [] # equivalences that should be converted to "Const" constraints
    
    # tabulate parameters used in equivalences by type
    depVarList = []  # list of all dependent parameters in equivalences
    indepVarList = []  # list of all independent parameters in equivalences
    for cnum,(varlist,mapvars,multarr,invmultarr) in enumerate(zip(
            dependentParmList,indParmList,arrayList,invarrayList)):
        #if multarr is not None: continue # equivalence
        indepVarList += [mv for mv in mapvars if mv not in indepVarList]
        depVarList += [v for v in varlist if v not in depVarList]

    # process equivalences: make a list of dependent and independent vars
    #    and check for repeated uses (repetition of a parameter as an
    #    independent var is OK)
    # look for parameters in equivalences that are used more than once as dependent parameters
    seenOnce = []
    for cnum,(varlist,multarr) in enumerate(zip(dependentParmList,arrayList)):
        if multarr is not None: continue # equivalences only
        for v in varlist:
            if v not in seenOnce:
                seenOnce.append(v)
            elif v not in multdepVarList:
                multdepVarList.append(v)                

    # scan through equivalences looking for other "dual uses". Stop when no new ones are found
    changed = True
    count = 0
    while changed:
        changed = False
        count += 1
        if count > 1000:
            raise ConstraintException("Too many loops in CheckEquivalences")
        
        # look for repeated dependent vars
        convVarList = [] # parameters in equivalences to be made into constraints
        for cnum,(varlist,mapvars,multarr,invmultarr) in enumerate(zip(
            dependentParmList,indParmList,arrayList,invarrayList)):
            if multarr is not None: continue # equivalences only
            if cnum in convertList:
                convVarList += [v for v in mapvars+varlist if v not in convVarList and type(v) is not float]
                continue
            
            # identify equivalences that need to be converted to constraints.
            #  Where parameters:
            #     are used both in equivalences & constraints,
            #     are used as dependent multiple times or
            #     where are used as both dependent and independent (chained)
            msg = False
            for v in mapvars:
                if v in constrVarList+convVarList:
                    changed = True
                    msg = True
                    warnEqv("Independent parameter "+str(v)+' used in constraint',cnum)
                    if cnum not in convertList: convertList.append(cnum)
            for v in varlist:
                if v in multdepVarList:
                    changed = True
                    msg = True
                    warnEqv("Dependent parameter "+str(v)+' repeated',cnum)
                    if cnum not in convertList: convertList.append(cnum)
                elif v in indepVarList:
                    changed = True
                    msg = True
                    warnEqv("Dependent parameter "+str(v)+' used elsewhere as independent',cnum)
                    if cnum not in convertList: convertList.append(cnum)
                elif v in constrVarList+convVarList:
                    changed = True
                    msg = True
                    warnEqv("Dependent parameter "+str(v)+' used in constraint',cnum)
                    if cnum not in convertList: convertList.append(cnum)
            if msg:
                warnEqv('Converting to "Constr"',cnum)
            
    global unvariedParmsList
    unvariedParmsList = []  # parameters in equivalences that are not varied
    # scan equivalences: look for holds
    for cnum,(varlist,mapvars,multarr,invmultarr) in enumerate(zip(
        dependentParmList,indParmList,arrayList,invarrayList)):
        if multarr is not None: continue # not an equivalence
        if cnum in convertList: continue

        # look for holds
        gotHold = False
        holdList = []
        for v in varlist+mapvars:
            if v in holdParmList:
                gotHold = True
            elif type(v) is not float:
                holdList.append(v)
        if gotHold:
            if holdList:
                msg = '  Some parameters set as "Hold"; setting remainder as "Hold": '
                for i,var in enumerate(holdList):
                    if i != 0: msg += ", "
                    msg += var
                    StoreHold(var,'Equiv fixed')
            else:
                msg = '  All parameters set as "Hold" '
            msg += "\n   and ignoring equivalence"
            warnEqv(msg,cnum)
            removeList.append(cnum)
            continue
        
        # look for unvaried parameters
        gotVary = False
        gotNotVary = False
        holdList = []
        for v in varlist+mapvars:
            if v in varyList:
                gotVary = True
                holdList.append(v)
            elif type(v) is not float:
                gotNotVary = True
                if v not in unvariedParmsList: unvariedParmsList.append(v)
        if gotNotVary:  # at least some unvaried parameters
            if gotVary:  # mix of varied and unvaried parameters
                msg = '\nSome parameters not varied; setting remainder as "Hold": '
                for i,var in enumerate(holdList):
                    if i != 0: msg += ", "
                    msg += var
                    StoreHold(var,'Equiv fixed')
            else:
                msg = 'No parameters varied '
            msg += "\nand ignoring equivalence"
            warnEqv(msg,cnum)
            removeList.append(cnum)
            continue
            
        # look for undefined or zero multipliers
        holdList = []
        drop = 0
        for v,m in zip(varlist,invmultarr):    
            if parmDict is not None and v not in parmDict:
                if v not in undefinedVars: undefinedVars.append(v)
                if v not in dropVarList: dropVarList.append(v)
                drop += 1
            elif m == 0:
                warnEqv("Parameter "+str(v)+" has a zero multiplier, dropping",cnum)
                if v not in dropVarList: dropVarList.append(v)
                drop += 1
            else:
                holdList.append(v)
        if drop == len(varlist):
            warnEqv("No dependent parameters defined, ignoring equivalence",cnum)
            removeList.append(cnum)
            continue
        for mv in mapvars:
            if type(mv) is float: continue
            if parmDict is not None and mv not in parmDict:
                # independent parameter is undefined, but some dependent parameters are defined
                # hold them
                if mv not in undefinedVars: undefinedVars.append(mv)
                msg = "Parameter(s) "+str(mv)
                for v in varlist:
                    if v in dropVarList:
                        msg += ', ' + v
                msg += "not defined in this refinement\n"
                msg = "Setting holds for: "
                for i,var in enumerate(holdList):
                    if i != 0: msg += ", "
                    msg += var
                warnEqv(msg,cnum)
                drop += 1
        if drop: # independent var and at least one dependent variable is defined
            msg = "Dropping undefined parameter(s) "
            i = 0
            for v in varlist:
                if v in dropVarList:
                    if i != 0: msg += ', '
                    i += 1
                    msg += v
            warnEqv(msg,cnum)
            msg = "Some parameters not defined. Setting holds for: "
            for i,var in enumerate(holdList):
                if i != 0: msg += ", "
                msg += var
                StoreHold(var,'Equiv fixed')
            warnEqv(msg,cnum)
    
    # Convert equivalences where noted
    for cnum,varlist in enumerate(dependentParmList):
        if cnum not in convertList: continue
        indvar = indParmList[cnum][0]
        # msg = '\nChanging equivalence:\n    ' + _showEquiv(
        #     dependentParmList[cnum],indParmList[cnum],invarrayList[cnum])
        for dep,mult in zip(dependentParmList[cnum],invarrayList[cnum]):
            constrDict += [{indvar:-1.,dep:mult[0]}]
            fixedList += ['0.0']
        #msg += '\n  to constraint(s):'
        #msg += '\n    ' + _FormatConstraint(constrDict[-1],fixedList[-1])
        removeList.append(cnum)
    # Drop equivalences where noted
    if removeList:
        for i in sorted(set(removeList),reverse=True):
            del dependentParmList[i],indParmList[i],arrayList[i],invarrayList[i],symGenList[i]
    # Drop variables from remaining equivalences
    for cnum,varlist in enumerate(dependentParmList):
        for j,v in enumerate(varlist):
            drop = []
            if v in dropVarList:
                drop.append(j)
        if drop:
            for j in sorted(drop,reverse=True):
                del indParmList[cnum][j]
    return warninfo['msg']

def ProcessConstraints(constList,seqmode='use-all',seqhst=None):
    """Interpret the constraints in the constList input into a dictionary, etc.
    All :class:`GSASIIobj.G2VarObj` objects are mapped to the appropriate
    phase/hist/atoms based on the object internals (random Ids). If this can't be
    done (if a phase has been deleted, etc.), the variable is ignored.
    If the constraint cannot be used due to too many dropped variables,
    it is counted as ignored. In the case of sequential refinements, 
    the current histogram number is substituted for a histogram number of "*".

    NB: this processing does not include symmetry imposed constraints
    
    :param list constList: a list of lists where each item in the outer list
      specifies a constraint of some form, as described in the :mod:`GSASIIobj`
      :ref:`Constraint definitions <Constraint_definitions_table>`.
    :param str seqmode: one of 'use-all', 'wildcards-only' or 'auto-wildcard'. 
       When seqmode=='wildcards-only' then any constraint with a numerical 
       histogram number is skipped. With seqmode=='auto-wildcard',
       any non-null constraint number is set to the selected histogram.
    :param int seqhst: number for current histogram (used for 
      'wildcards-only' or 'auto-wildcard' only). Should be None for 
      non-sequential fits.

    :returns:  a tuple of (constrDict,fixedList,ignored) where:
      
      * constrDict (list of dicts) contains the constraint relationships
      * fixedList (list) contains the fixed values for each type
        of constraint.
      * ignored (int) counts the number of invalid constraint items
        (should always be zero!)
    """
    constrDict = []
    fixedList = []
    ignored = 0
    namedVarList = []
    for constr in constList:
        terms = copy.deepcopy(constr[:-3]) # don't change the tree contents
        # deal with wildcards in sequential fits
        if seqmode == 'wildcards-only' and seqhst is not None:
            skip = False
            for term in terms:
                if term[1].histogram == '*':
                    term[1] = term[1].varname(seqhst)
                elif term[1].histogram:
                    skip = True
            if skip: continue
        elif seqmode == 'auto-wildcard' and seqhst is not None:
            for term in terms:
                term[1] = term[1].varname(seqhst)
        elif seqhst is not None:
            for term in terms:
                if term[1].histogram == '*':
                    term[1] = term[1].varname(seqhst)
                # else:
                #     term[1] = term[1].varname()   # does this change anything???
        # separate processing by constraint type
        if constr[-1] == 'h':
            # process a hold
            var = str(terms[0][1])
            if '?' not in var:
                StoreHold(var,'User supplied')
            else:
                ignored += 1
        elif constr[-1] == 'f':
            # process a new variable
            fixedList.append(None)
            D = {}
            varyFlag = constr[-2]
            varname = constr[-3]
            for term in terms:
                var = str(term[1])
                if '?' not in var:
                    D[var] = term[0]
            # add extra dict terms for input variable name and vary flag
            if varname is None: # no assigned name, create one
                global consNum
                varname = str(consNum) 
                consNum += 1
            else:
                varname = str(varname) # in case this is a G2VarObj
            if '::' in varname:
                D['_name'] = varname.replace('::','::nv-')
            else:
                D['_name'] = '::nv-' + varname
            D['_name'] = G2obj.MakeUniqueLabel(D['_name'],namedVarList)
            D['_vary'] = varyFlag == True # force to bool
            constrDict.append(D)
        elif constr[-1] == 'c': 
            # process a constraint equation
            D = {}
            for term in terms:
                var = str(term[1])
                if '?' not in var:
                    D[var] = term[0]
            if len(D) >= 1:
                fixedList.append(float(constr[-3]))
                constrDict.append(D)
            else:
                ignored += 1
        elif constr[-1] == 'e':
            # process an equivalence
            firstmult = None
            eqlist = []
            for term in terms:
                if term[0] == 0: term[0] = 1.0
                var = str(term[1])
                if '?' in var: continue
                if firstmult is None:
                    firstmult = term[0]
                    firstvar = var
                else:
                    eqlist.append([var,firstmult/term[0]])
            if len(eqlist) > 0:
                StoreEquivalence(firstvar,eqlist,False)
            else:
                ignored += 1
        else:
            ignored += 1
    return constrDict,fixedList,ignored

def StoreHold(var,holdType=None):
    '''Takes a variable name and prepares it to be removed from the 
    refined variables.

    Called with user-supplied constraints by :func:`ProcessConstraints`.
    At present symGen is not used, but could be set up to track Holds generated
    by symmetry.
    '''
    global holdParmList,holdParmType
    if var not in holdParmList:
        holdParmList.append(var)
        if holdType: holdParmType[var] = holdType
    
def StoreEquivalence(independentVar,dependentList,symGen=True):
    '''Takes a list of dependent parameter(s) and stores their
    relationship to a single independent parameter (independentVar).

    Called with user-supplied constraints by :func:`ProcessConstraints`,
    with Pawley constraints from :func:`GSASIIstrIO.GetPawleyConstr`, 
    with Unit Cell constraints from :func:`GSASIIstrIO.cellVary`
    with symmetry-generated atom constraints from :func:`GSASIIstrIO.GetPhaseData`

      There is no harm in using StoreEquivalence with the same independent variable::

       StoreEquivalence('x',('y',))
       StoreEquivalence('x',('z',))

      but the same outcome can be obtained with a single call::

       StoreEquivalence('x',('y','z'))

      The latter will run more efficiently. 


      Note that mixing independent and dependent variables, such as:: 

        StoreEquivalence('x',('y',))
        StoreEquivalence('y',('z',))

      is a poor choice. The module will attempt to fix this by transforming the equivalence to a 
      "Const" constraint.

    :param str independentVar: name of master parameter that will be used to determine the value
      to set the dependent variables

    :param list dependentList: a list of parameters that will set from
         independentVar. Each item in the list can be a string with the parameter
         name or a tuple containing a name and multiplier:
         ``['::parm1',('::parm2',.5),]``

    '''
    
    global dependentParmList,arrayList,invarrayList,indParmList,symGenList
    mapList = []
    multlist = []
    allfloat = True
    for var in dependentList:
        if isinstance(var, str):
            mult = 1.0
        elif len(var) == 2:
            var,mult = var
        else:
            raise Exception("Cannot parse "+repr(var) + " as var or (var,multiplier)")
        mapList.append(var)
        try:            
            multlist.append(tuple((float(mult),)))
        except:
            allfloat = False
            multlist.append(tuple((mult,)))
    # added relationships to stored values
    arrayList.append(None)
    if allfloat:
        invarrayList.append(np.array(multlist))
    else:
        invarrayList.append(multlist)
    indParmList.append(list((independentVar,)))
    dependentParmList.append(mapList)
    symGenList.append(symGen)
    return

def SubfromParmDict(s,prmDict):
    '''Process a string as a multiplier and convert it to a float value. This
    is done by subsituting any GSAS-II parameter names that appear in the 
    string that have associated values in the parameter dict with the value 
    for that parameter. 

    :param str s: a string to be converted to a value
    :param dict prmDict: a dictionary with keys as GSAS-II parameter names
      and values the corresponding parameter value.
    :returns: the evaluated expression as a float.
    '''
    # TODO: perhaps SubfromParmDict should be called to convert the
    # fixed-val in constraint equations from strings to values.
    for key in prmDict:
        if key in s:
            s = s.replace(key,str(prmDict[key]))
    return eval(s)

def EvaluateMultipliers(constList,*dicts):
    '''Convert multipliers for constraints and equivalences that are specified
    as strings into values. The strings can specify values in the parameter dicts as 
    well as normal Python functions, such as "2*np.cos(0::Ax:2/2.)"
    
    :param list constList: a list of dicts containing constraint expressions
    :param \*dicts: one or more dicts containing GSAS-II parameters and their values 
       can be specified
    :returns: an empty string if there were no errors, or an error message listing
       the strings that could not be converted.
    '''
    prmDict = {}
    for d in dicts: prmDict.update(d) # combine all passed parameter dicts
    problemList = ""
    # loop through multipliers in contraint expressions
    for const in constList:
        for key in const:
            if key.startswith('_'): continue
            try: # is this already a float, etc? 
                1+const[key]
                continue
            except:
                pass
            try:
                newval = SubfromParmDict(const[key][:],prmDict)
                if GSASIIpath.GetConfigValue('debug'):
                    print('Changing ',const[key],'to',newval)
                const[key] = newval                
            except:
                if problemList: problemList += ", "
                problemList += const[key]
    # loop through multipliers in equivalences
    global arrayList,invarrayList
    for i,(a,valList) in enumerate(zip(arrayList,invarrayList)):
        if a is not None: continue # ignore if not equiv
        try:
            valList.shape
            continue # ignore if already a numpy array
        except:
            pass
        repList = []
        for v in valList:
            try:
                1+v[0]
                repList.append(tuple((v[0],)))
                continue
            except:
                pass
            try:
                newval = SubfromParmDict(v[0][:],prmDict)
                if GSASIIpath.GetConfigValue('debug'):
                    print('Changing ',v[0],'to',newval)
                repList.append(tuple((newval,)))
            except:
                if problemList: problemList += ", "
                problemList += v[0]
                repList.append(tuple(('error',)))
        invarrayList[i] = np.array(repList)
    return problemList

def GetDependentVars():
    '''Return a list of dependent variables: e.g. parameters that are
    constrained in terms of other parameters

    :returns: a list of parameter names

    '''
    global dependentVars
    return dependentVars

def GetIndependentVars():
    '''Return a list of independent variables: e.g. parameters that are
    slaved to other parameters by constraints

    :returns: a list of parameter names

    '''
    global independentVars
    return independentVars

def PrintIndependentVars(parmDict,varyList,sigDict,PrintAll=False,pFile=None):
    '''Print the values and uncertainties on the independent parameters'''
    global dependentParmList,arrayList,invarrayList,indParmList
    printlist = []
    mvs = GetIndependentVars()
    for i,name in enumerate(mvs):
        if PrintAll or name in varyList:
            sig = sigDict.get(name)
            printlist.append([name,parmDict[name],sig])
    if len(printlist) == 0: return
    s1 = ''
    s2 = ''
    s3 = ''
    pFile.write(130*'-'+'\n')
    pFile.write("Parameters generated by constraints\n")
    printlist.append(3*[None])
    for name,val,esd in printlist:
        if len(s1) > 120 or name is None:
            pFile.write(''+'\n')
            pFile.write(s1+'\n')
            pFile.write(s2+'\n')
            pFile.write(s3+'\n')
            s1 = ''
            if name is None: break
        if s1 == "":
            s1 = ' name  :'
            s2 = ' value :'
            s3 = ' sig   :'
        s1 += '%15s' % (name)
        s2 += '%15.5f' % (val)
        if esd is None:
            s3 += '%15s' % ('n/a')
        else:
            s3 += '%15.5f' % (esd)

def getConstrError(constrLst,seqmode,seqhst):
    '''This is used to display error messages for constraints and 
    equivalence relations

    :parm list constrLst: a single constraint or equivalence as saved 
      in the data tree
      (see :ref:`constraint definitions <Constraint_definitions_table>`).
    :param str seqmode: one of 'use-all', 'wildcards-only' or 'auto-wildcard'
    :param int seqhst: number for current histogram (used for 
      'wildcards-only' or 'auto-wildcard' only). Should be None for 
      non-sequential fits.

    :returns: error, msg where error (bool) is True if the 
      constraint/equivalence creates an error, msg (str) can be a warning 
      or an error
    '''
    msg = ''
    note = ''
    terms = copy.deepcopy(constrLst[:-3])
    if seqmode == 'wildcards-only' and seqhst is not None:
        if constrLst[-1] == 'e':
            msg = 'equivalence'
        else:
            msg = 'constraint'
        for term in terms:
            if term[1].histogram == '*':
                term[1] = term[1].varname(seqhst)
            elif term[1].histogram:
                return False,"Ignoring non-wildcard "+msg, "Ignore"
    elif seqmode == 'auto-wildcard' and seqhst is not None:
        for term in terms:
            term[1] = term[1].varname(seqhst)
    else:
        for term in terms:
            if term[1].histogram == '*':
                if seqhst is None:
                    msg = "Parameter "+str(terms[0][1])+" contains a wildcard, which are used only sequential refinements. Constraint ignored."
                    return False,msg, "Ignored"
                else:
                    term[1] = term[1].varname(seqhst)
            else:
                term[1] = term[1].varname()
    if constrLst[-1] == 'e':
        # conflicting uses
        if terms[0][1] in constrVarList+convVarList:
            msg = "Parameter "+str(terms[0][1])+" used in constraint. To be recast as constraint."
            return False,msg, "Recast as constraint."
        varList = []
        for m,v in terms:
            if v in constrVarList+convVarList:
                varList.append(str(v))
        if varList:
            msg = "Parameter(s) used in constraint: "
            for i,v in enumerate(varList):
                if i != 0: msg += ', '
                msg += v
            return False,msg,"Recast as constraint."
        varList = []
        for m,v in terms[1:]:
            if v in indepVarList:
                varList.append(str(v))
        if varList:
            msg = "Parameter(s) used elsewhere as independent: "
            for i,v in enumerate(varList):
                if i != 0: msg += ', '
                msg += v
            return False,msg,"Recast as constraint."
        varList = []
        for m,v in terms[1:]:
            if v in multdepVarList:
                varList.append(str(v))
        if varList:
            msg += "Parameter(s) repeated as dependent: "
            for i,v in enumerate(varList):
                if i != 0: msg += ', '
                msg += v
            return False,msg,"Recast as constraint."

        # zero multiplier
        varList = []
        valid = 0
        for m,v in terms:
            if m == 0:
                varList.append(str(v))
            else:
                valid += 1
        if varList and valid > 1:
            msg += "Parameter(s) with zero multipliers: "
            for i,v in enumerate(varList):
                if i != 0: msg += ', '
                msg += v
            msg += " will be ignored"
        elif varList:
            msg += "Parameter(s) with zero multipliers:"
            for i,v in enumerate(varList):
                if i != 0: msg += ', '
                msg += v
            return False,msg,"Equivalence Ignored."
        
        # hold parameters
        s = ''
        for m,v in terms:
            if v in holdParmList and (
                    "User supplied" in holdParmType.get(v,'') or
                    "symmetry" in holdParmType.get(v,'') or
                    "rigid body" in holdParmType.get(v,'')):
                if s: s += ', '
                s += str(v)
        if s:
            if msg: msg += '; '
            msg += "Parameter(s) set as Hold: "+s
            msg += "\nThis equivalence will be ignored; all parameters are held."
            return False,msg,'Has holds: Not varied.'
                
        # unrefined parameters
        gotVary = False
        gotNotVary = False
        s = ''
        for m,v in terms:
            if v in unvariedParmsList:
                gotNotVary = True
                if s: s += ', '
                s += str(v)
            else:
                gotVary = True
        if gotNotVary and gotVary:  # mix of varied and unvaried parameters
            if msg: msg += '. '
            msg += 'Unvaried parameter(s): '+s+"; remainder set Hold. All parameters fixed."
            return False,msg, 'Equivalence Ignored.'
        elif gotNotVary:
            if msg: msg += '. '
            msg += 'All parameters not varied. Equivalence Ignored.'
            return False,msg, 'Equivalence Ignored.'
            
        # undefined parameters
        undef = 0
        s = ''
        for m,v in terms[1:]:
            if v in undefinedVars:
                undef += 1
                if s: s += ', '
                s += str(v)
        if undef == len(terms[1:]):
            msg += 'Ignored: None of the dependent parameters are defined'
        elif terms[0][1] in undefinedVars:
            if s:
                s = terms[0][1] + ', ' + s
            else:
                s = terms[0][1]
            msg += 'Undefined parameter(s): '+s+'. Remainder will be fixed'
        elif undef:
            msg += 'Undefined parameter(s): '+s+' will be dropped'
    elif constrLst[-1] == 'h':
        v = terms[0][1]
        if v in undefinedVars: return False,"Parameter is undefined","Ignored"
        if v in unvariedParmsList: return False,"Parameter is not refined","Ignored"
    else:
        # check for post-grouping errors in new var & constr. eq. groups
        for m,v in terms:
            if v in groupErrors:
                return True,'Constraint singularity: see error listing','Singular'
        zeroList = []
        toBeUsed = []
        undef = []
        unvar = []
        hold = []
        for m,v in terms: # check for zero multiplier, undefined, unvaried or hold
            if constrLst[-1] == 'f':    # new var expressions; 
                if m == 0:
                    zeroList.append(str(v))
                elif v in undefinedVars:
                    undef.append(str(v))
                elif (v in holdParmList and constrLst[-2] and
                        "User supplied" in holdParmType.get(v,'') or
                        "symmetry" in holdParmType.get(v,'') or
                        "rigid body" in holdParmType.get(v,'')):                        
                    hold.append(str(v))
            else:                        # constraint equation
                if m == 0:
                    zeroList.append(str(v))
                elif v in undefinedVars:
                    undef.append(str(v))
                elif v in unvariedParmsList:
                    unvar.append(str(v))
                elif v in holdParmList and holdParmType.get(v,'') != 'dependent param':
                    hold.append(str(v))
                else:
                    toBeUsed.append(str(v))

        s = ''
        for v in zeroList:
            if s: s += ', '
            s += str(v)
        if s:
            if msg: msg += '; '
            msg += "Parameter(s) with zero multipliers: "+s
        s = ''
        for v in undef:
            if s: s += ', '
            s += str(v)
        if s:
            if msg: msg += '; '
            msg += "Undefined parameter(s): "+s
        s = ''
        for v in unvar:
            if s: s += ', '
            s += str(v)
        if s:
            if msg: msg += '; '
            msg += "Unrefined parameter(s) will be dropped: "+s
        s = ''
        for v in hold:
            if s: s += ', '
            s += str(v)
        if s:
            if msg: msg += '; '
            msg += '"Hold" parameter(s): '+s
        if hold and constrLst[-1] == 'f':
            if msg: msg += '; '
            msg += "\nNew var with holds cannot be processed; constraint ignored."
            note = 'Ignored'
        elif undef and toBeUsed:
            s = ''
            for v in toBeUsed:
                if s: s += ', '
                s += str(v)
            if msg: msg += '; '
            msg += "Adding Holds on "+s+"; Constraint Ignored."
            note = 'Ignored'
        elif undef or (len(toBeUsed) == 0 and (zeroList or unvar or hold)):
            if msg: msg += '; '
            msg += "No parameters remain. Constraint Ignored."
            note = 'Ignored'
        elif len(toBeUsed) == 1 and (zeroList or unvar or hold):
            if msg: msg += '; '
            msg += "One parameter is retained; converted to fixed value."
            note = 'Converted'
        elif zeroList or unvar or hold:
            if msg: msg += ': '
            msg += '\nConstraint adjusted for fixed values and retained.'
            note = 'Parameter(s) removed'
    return False,msg,note
        
def ComputeDepESD(covMatrix,varyList,allvars=True):
    '''Compute uncertainties for dependent parameters from independent ones
    returns a dictionary containing the esd values for dependent parameters
    
    :param np.array covMatrix: the full covariance matrix
    :param list varyList: the names of the variables matching the columns 
      and rows in covMatrix
    :param bool allvars: When True (default) s.u. values for all parameters 
      are placed in the returned dict. When False the number of s.u. values
      attempts to match the number of refined degrees of freedom. The s.u.'s 
      for dependent params from equivalences are not computed and 
      the number of dependent params from new var and generated var
      constraints matches the number of refined independent parameters. 
    '''
    sigmaDict = {}
    for varlist,mapvars,multarr,invmultarr in zip(dependentParmList,indParmList,arrayList,invarrayList):
        if multarr is None and not allvars: continue # probably not needed
        varied = 0
        # get the v-covar matrix for independent parameters 
        vcov = np.zeros((len(mapvars),len(mapvars)))
        for i1,name1 in enumerate(mapvars):
            if name1 not in varyList: continue
            varied += 1
            iv1 = varyList.index(name1)
            for i2,name2 in enumerate(mapvars):
                if name2 not in varyList: continue
                iv2 = varyList.index(name2)
                vcov[i1][i2] = covMatrix[iv1][iv2]
        # vec is the vector that multiplies each of the independent values
        for i,(v,vec) in enumerate(zip(varlist,invmultarr)):
            if i == varied: break
            sigmaDict[v] = np.sqrt(np.inner(vec.T,np.inner(vcov,vec)))
    return sigmaDict

def _FormatConstraint(RelDict,RelVal):
    '''Formats a Constraint or Function for use in a convenient way'''
    linelen = 65
    s = [""]
    for var,val in RelDict.items():
        if var.startswith('_'): continue
        if len(s[-1]) > linelen: s.append(' ')
        m = val
        if s[-1] != "" and m >= 0:
            s[-1] += ' + '
        elif s[-1] != "":
            s[-1] += ' - '
            m = abs(m)
        if m == 1:
            s[-1] += '%s '%var
        else:
            s[-1] += '%.3f*%s '%(m,var)
    if len(s[-1]) > linelen: s.append(' ')
    if RelVal is None:
        s[-1] += ' = New variable'
    else:
        s[-1] += ' = ' + str(RelVal)
    s1 = ''
    for s2 in s:
        if s1 != '': s1 += '\n  '
        s1 += s2
    return s1

def _showEquiv(varlist,mapvars,invmultarr,longmsg=False):
    '''Format an equivalence relationship
    note that 
    varlist,           mapvars,     invmultarr 
    are elements of
    dependentParmList, indParmList, invarrayList
    '''
    for i,mv in enumerate(mapvars):
        s1 = str(mv)
        if not longmsg:
            s1 += ' ==> '
        elif len(varlist) == 1:
            s1 += ' is equivalent to '
        else:
            s1 += ' is equivalent to parameters: '
        j = 0
        for v,m in zip(varlist,invmultarr):
            if debug: print ('v,m[0]: ',v,m[0])
            if len(s1.split('\n')[-1]) > 75: s1 += '\n        '
            if j > 0: s1 += ' & '
            j += 1
            s1 += str(v)
            if m != 1:
                s1 += " / " + str(m[0])
    return s1

def VarRemapShow(varyList=None,inputOnly=False,linelen=60):
    '''List out the saved relationships. This should be done after the constraints have been
    defined using :func:`StoreEquivalence`, :func:`GroupConstraints` and :func:`GenerateConstraints`.

    :returns: a string containing the details of the contraint relationships
    '''
    if varyList is None:
        varyList = saveVaryList
    s = ''
    if len(dependentVars) > 0:
        s += '\nDependent parameters, determined by constraints (also set as Hold):\n'
        for v in sorted(dependentVars):
            s += '    ' + v + '\n'
    first = True
    for v in sorted(holdParmList):
        if v in dependentVars: continue
        if v in holdParmType: v += '\t'+holdParmType[v]
        if first:
            s += '\nParameters set as Hold:\n'
        first = False
        s += '    ' + v + '\n'
        
    userOut = ''
    symOut = ''
    consOut = ''
    varOut = ''
    freeOut = ''
    global dependentParmList,arrayList,invarrayList,indParmList,symGenList

    for varlist,mapvars,multarr,invmultarr,symFlag in zip(
        dependentParmList,indParmList,arrayList,invarrayList,symGenList):
        for i,mv in enumerate(mapvars):
            if multarr is None:
#                s1 = '  ' + str(mv) + ' is equivalent to parameter(s): '
                if len(varlist) == 1:
                    s1 = '   ' + str(mv) + ' is equivalent to '
                else:
                    s1 = '   ' + str(mv) + ' is equivalent to parameters: '
                j = 0
                for v,m in zip(varlist,invmultarr):
                    if debug: print ('v,m[0]: ',v,m[0])
                    if len(s1.split('\n')[-1]) > 75: s1 += '\n        '
                    if j > 0: s1 += ' & '
                    j += 1
                    s1 += str(v)
                    if m != 1:
                        s1 += " / " + str(m[0])
                if symFlag:
                    symOut += s1 + '\n'
                else:
                    userOut += s1 + '\n'
                continue
            ln = ''
            # if mv in varyList: 
            #     lineOut = '  (V)* {} = '.format(mv)
            # else:
            #     lineOut = '  {} = '.format(mv)
            j = 0 
            lineOut = '  {} = '.format(mv)
            for (m,v) in zip(multarr[i,:],varlist):
                if m == 0: continue
                if m < 0:
                    lineOut += ' - '
                    m *= -1
                elif j != 0:
                    lineOut += ' + '
                j += 1
                if len(lineOut) > linelen:
                    ln += lineOut
                    lineOut = '\n  '
                if m == 1:
                    lineOut += '{}'.format(v)
                else:
                    lineOut += '({:.4g} * {})'.format(m,v)
            if mv in varyList: 
                lineOut += '\t *VARIED*'
            ln += lineOut
            
            if type(mv) is float:
                consOut += ln + '\n'
            elif '::nv-' in mv:
                varOut += ln + '\n'
            else:
                freeOut += ln + '\n'
    if userOut:
        s += '\nEquivalences:\n' + userOut
    if consOut:
        s += '\nConstraint Equations:\n' + consOut
    if varOut:
        s += '\nNew Variable assignments:\n' + varOut
    if freeOut:
        s += '\nGenerated Degrees of Freedom:\n' + freeOut
    if symOut:
        s += '\nSymmetry-generated equivalences:\n' + symOut
    if not (userOut or consOut or varOut or symOut):
        return s + '\nNo constraints or equivalences in use'
    elif inputOnly:
        return s
        
    s += '\nInverse parameter mapping relations:\n'
    lineDict = {} # store so we can sort them
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        varied = False
        for i,mv in enumerate(varlist):
            lineOut = '  {} = '.format(mv)
            j = 0 
            for m,v in zip(invmultarr[i,:],mapvars):
                if m == 0: continue
                if v == 0: continue
                if v in varyList: varied = True
                if m < 0:
                    lineOut += ' - '
                    m *= -1
                elif j != 0:
                    lineOut += ' + '
                j += 1
                if len(lineOut) > linelen:
                    s += lineOut
                    lineOut = '\n  '
                if m == 1:
                    lineOut += '{}'.format(v)
                else:
                    try:
                        lineOut += '{:.4g}'.format(m*v)
                    except:
                        lineOut += '({:.4g} * {})'.format(m,v)
            if j == 0: lineOut += '0'
            if varied: lineOut += '\t *VARIED*'
            lineDict[mv] = lineOut
    for key in sorted(lineDict):
        s += lineDict[key] + '\n'
    return s

def GetSymEquiv(seqmode,seqhistnum):
    '''Return the automatically generated (equivalence) relationships.

    :returns: a list of strings containing the details of the contraint relationships
    '''
    symout = []
    symerr = []
    symhelp = []
    global dependentParmList,arrayList,invarrayList,indParmList,symGenList

    for varlist,mapvars,multarr,invmultarr,symFlag in zip(
        dependentParmList,indParmList,arrayList,invarrayList,symGenList):
        if not symFlag: continue
        for i,mv in enumerate(mapvars):
            cnstr = [[1,G2obj.G2VarObj(mv)]]
            if multarr is None:
                s1 = ''
                s2 = ' = ' + str(mv)
                j = 0
                helptext = 'Variable {:} '.format(mv) + " ("+ G2obj.fmtVarDescr(mv) + ")"
                if len(varlist) == 1:
                    cnstr.append([invmultarr[0][0],G2obj.G2VarObj(varlist[0])])
                    # format the way Bob prefers
                    if invmultarr[0][0] == 1: 
                        s1 = str(varlist[0]) + ' = ' + str(mv)
                    else:
                        s1 = str(varlist[0]) + ' = ' + str(
                            invmultarr[0][0]) + ' * '+ str(mv)
                    s2 = ''
                    
                    m = 1./invmultarr[0][0]
                    var1 = str(varlist[0])
                    helptext += "\n\nis equivalent to "
                    if m == 1:
                        helptext += '\n  {:} '.format(var1) + " ("+ G2obj.fmtVarDescr(var1) + ")"
                    else:
                        helptext += '\n  {:3g} * {:} '.format(m,var1) + " ("+ G2obj.fmtVarDescr(var1) + ")"
                else:
                    helptext += "\n\nis equivalent to the following:"
                    for v,m in zip(varlist,invmultarr):
                        cnstr.append([m,G2obj.G2VarObj(v)])
                        #if debug: print ('v,m[0]: ',v,m[0])
                        if len(s1.split('\n')[-1]) > 75: s1 += '\n        '
                        if j > 0: s1 += ' =  '
                        j += 1
                        s1 += str(v)
                        if m != 1:
                            s1 += " / " + str(m[0])
                            helptext += '\n  {:3g} * {:} '.format(m,v) + " ("+ G2obj.fmtVarDescr(v) + ")"
                        else:
                            helptext += '\n  {:} '.format(v) + " ("+ G2obj.fmtVarDescr(v) + ")"
                err,msg,note = getConstrError(cnstr+[None,None,'e'],seqmode,seqhistnum)
                symerr.append([msg,note])
                symout.append(s1+s2)
                symhelp.append(helptext)
            else:
                s = '  %s = ' % mv
                j = 0
                for m,v in zip(multarr[i,:],varlist):
                    if m == 0: continue
                    if j > 0: s += ' + '
                    j += 1
                    s += '(%s * %s)' % (m,v)
                print ('unexpected sym op='+s)
    return symout,symerr,symhelp

def Dict2Deriv(varyList,derivDict,dMdv):
    '''Compute derivatives for Independent Parameters from the
    derivatives for the original parameters

    :param list varyList: a list of parameters names that will be varied

    :param dict derivDict: a dict containing derivatives for parameter values keyed by the
      parameter names.

    :param list dMdv: a Jacobian, as a list of np.array containing derivatives for dependent
      parameter computed from derivDict

    '''
    global dependentParmList,arrayList,invarrayList,indParmList,invarrayList
    for varlist,mapvars,multarr,invmultarr in zip(dependentParmList,indParmList,arrayList,invarrayList):
        for i,name in enumerate(mapvars):
            # grouped parameters: need to add in the derv. w/r
            # dependent variables to the independent ones
            if name not in varyList: continue # skip if independent var not varied
            if multarr is None:
                if debug: print ('start dMdv for',name,dMdv[varyList.index(name)])
                for v,m in zip(varlist,invmultarr):
                    if m[0] == 0: continue
                    dMdv[varyList.index(name)] += derivDict[v]/ m[0] 
            else:
                for v,m in zip(varlist,invmultarr[:,i]):
                    if m == 0: continue
                    dMdv[varyList.index(name)] += m * derivDict[v]

def Map2Dict(parmDict,varyList):
    '''Updates the parameter dictionary and the varyList using the 
    equivalence and constraint input. This should be called at least once, after 
    the constraints have been defined using :func:`StoreEquivalence`, 
    :func:`GroupConstraints` and :func:`GenerateConstraints` and before any 
    parameter refinement is done. 

    This completes the parameter dictionary by defining values for parameters 
    created by constraints based on the constraints that define them 
    using the values for the current parameters. It also removes all dependent 
    variables from the varyList

    :param dict parmDict: a dict containing parameter values keyed by the
      parameter names. For new variables created by constraints, entries
      will be added to the dictionary, if not alreay present, or the 
      values will be recomputed.

    :param list varyList: a list of parameters names. Will be modified.
    '''
    # remove fixed parameters from the varyList
    for item in holdParmList:
        if item in varyList: varyList.remove(item)
            
    # process the independent parameters:
    # * remove dependent ones from varylist
    # * for equivalences apply the independent parameters onto dependent variables
    global dependentParmList,arrayList,invarrayList,indParmList
    for varlist,mapvars,multarr,invmultarr in zip(dependentParmList,indParmList,arrayList,invarrayList):
        for item in varlist: # TODO: is this still needed?
            if item in varyList: varyList.remove(item)
        if multarr is None:
            #for v,val in zip(  # shows values to be set
            #    varlist,
            #    np.dot(invmultarr,np.array([parmDict[var] for var in mapvars]))
            #    ): print('parmDict set',v,':',val)
            parmDict.update(zip(
                varlist,
                np.dot(invmultarr,np.array([parmDict[var] for var in mapvars]))
                ))

    # * for the created parameters, compute them from their dependents
    for varlist,mapvars,multarr in zip(dependentParmList,indParmList,arrayList):
        if multarr is None: continue
        # evaluate constraints in the forward direction
        A = np.array([parmDict[var] for var in varlist])
        z = zip(mapvars,np.dot(multarr,A))
        # add/replace in parameter dict
        parmDict.update([i for i in z if type(i[0]) is not float])
    global saveVaryList
    saveVaryList = copy.copy(varyList)
            
def Dict2Map(parmDict):
    '''Applies the constraints defined using :func:`StoreEquivalence`,
    :func:`GroupConstraints` and :func:`GenerateConstraints` by changing
    values in a dict containing the parameters. This should be
    done after refinement and before the parameters are used for 
    any computations

    :param dict parmDict: a dict containing parameter values keyed by the
      parameter names. After this is called, all the dependent variables
      will be updated based on constraints and equivalences.
    '''
    global dependentParmList,arrayList,invarrayList,indParmList
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        if invmultarr is None:  # is this needed?
            if GSASIIpath.GetConfigValue('debug'): 
                print('Why does this constraint have None for invmultarr?',varlist,mapvars)
            continue
        valslist = np.array([parmDict.get(var,var) for var in mapvars])
        #for v,val in zip(varlist,np.dot(invmultarr,np.array(valslist))): print(v,val) # shows what is being set
        parmDict.update(zip(varlist,np.dot(invmultarr,valslist)))
        
#======================================================================
# internal routines follow (these routines are unlikely to be called
# from outside the module)

def GramSchmidtOrtho(a,nkeep=0):
    '''Use the Gram-Schmidt process (http://en.wikipedia.org/wiki/Gram-Schmidt) to
    find orthonormal unit vectors relative to first row.

    If nkeep is non-zero, the first nkeep rows in the array are not changed
    
    input: 
       arrayin: a 2-D non-singular square array
    returns:
       a orthonormal set of unit vectors as a square array
    '''
    def proj(a,b):
        'Projection operator'
        return a*(np.dot(a,b)/np.dot(a,a))
    for j in range(nkeep,len(a)):
        for i in range(j):
            a[j] -= proj(a[i],a[j])
        if np.allclose(np.linalg.norm(a[j]),0.0):
            raise ConstraintException("Singular input to GramSchmidtOrtho")
        a[j] /= np.linalg.norm(a[j])
    return a

def _FillArray(sel,d,collist):
    '''Construct a n by n matrix [n = len(collist)]
    with the initial m rows [m = len(sel)] using the 
    relationships defined in the expressions dict, d.
    Since m may be smaller than n, the remaining rows
    are filled with rows that are tested to not create
    a singular matrix. 

    :param list sel:     a list of indices in dict d
    :param list d:       a list of dict's where each dict describes an 
      expression from a constraint equation or a new var
    :param list collist: a list parameter names.
    :returns: an n by n numpy.array matrix
    '''
    n = len(collist)
    m = len(sel)
    arr = np.zeros(2*[n,])
    # fill the top rows
    for i,cnum in enumerate(sel):
        for j,var in enumerate(collist):
            arr[i,j] = d[cnum].get(var,0)
    try:
        _RowEchelon(m,copy.copy(arr),collist)
    except:
        raise Exception('Initial constraints singular')
    for i in range(m,n):
        arr[i][i] = 1   # add a diagonal element
        try:
            _RowEchelon(i+1,copy.copy(arr),collist)
            continue
        except:
            pass
        for j in range(n):
            if j == i: continue
            arr[i][j] = 1   # add another element
            try:
                _RowEchelon(i+1,copy.copy(arr),collist)
                break
            except:
                pass
            arr[i][j] = -1   # try a different valuefor this element
            try:
                _RowEchelon(i+1,copy.copy(arr),collist)
                break
            except:
                arr[i][j] = 0   # reset to add another element
        else:
            raise Exception('Unable to create non-singular matrix')
    return arr

def _SwapColumns(i,m,v):
    '''Swap columns in matrix m as well as the labels in v 
    so that element (i,i) is replaced by the first non-zero element in row i after that element

    Throws an exception if there are no non-zero elements in that row
    '''
    for j in range(i+1,len(v)):
        if not np.allclose(m[i,j],0):
            m[:,(i,j)] = m[:,(j,i)]
            v[i],v[j] = v[j],v[i]
            return
    else:
        raise ConstraintException('Singular input')

def _RowEchelon(m,arr,collist):
    '''Convert the first m rows in Matrix arr to row-echelon form
    exchanging columns in the matrix and collist as needed.

    throws an exception if the matrix is singular because
    the first m rows are not linearly independent
    '''
    for i in range(m):
        if np.allclose(arr[i,i],0):
            _SwapColumns(i,arr,collist)
        arr[i,:] /= arr[i,i] # normalize row
        # subtract current row from subsequent rows to set values to left of diagonal to 0
        for j in range(i+1,m):
            arr[j,:] -= arr[i,:] * arr[j,i]

if __name__ == "__main__":
    d = [{'a': 0, 'b': 1.5,'d':0}, {'d': -1}]
    lbls = ['a','c','b','d']
    sel = (0,1)
    try:
        arr2 = _FillArray(sel,d,lbls)
    except Exception as err:
        if 'Initial' in str(err): 
            print('initial error')
        else:
            print('unable to extend matrix error')
        import sys; sys.exit()
    print(arr2)

    d = [{'a': 1, 'b': 1,}]
    lbls = ['a','b']
    sel = (0,)
    arr1 = _FillArray(sel,d,lbls)
    print(arr1)
    print(GramSchmidtOrtho(arr1,1))
    
    import sys; sys.exit()
    
    parmdict = {}
    constrDict = [
        {'0:12:Scale': 2.0, '0:11:Scale': 1.0, '0:14:Scale': 4.0, '0:13:Scale': 3.0, '0:0:Scale': 0.5},
        {'0:0:eA': 0.0},
        {'2::C(10,6,1)': 1.0, '1::C(10,6,1)': 1.0},
        {'1::C(10,0,1)': 1.0, '2::C(10,0,1)': 1.0},
        {'1::AUiso:0': 1.0, '0::AUiso:0': 1.0},
        {'0::A0': 0.0}
        ]
    fixedList = ['5.0', '0', None, None, '1.0', '0']
    StoreEquivalence('2::atomx:3',('2::atomy:3', ('2::atomz:3',2,), ))
    #StoreEquivalence('1::atomx:3',('2::atomx:3', ('2::atomz:3',2,), )) # error: dependent & independent vars mixed
    #StoreEquivalence('1::atomx:3',('2::atomy:3', ('2::atomz:3',2,), )) # error: dependent vars repeated
    #StoreEquivalence('0:1:eA',('0:0:eA',)) # error: equiv & fixed
    #StoreEquivalence('0:99:Scale',('0:12:Scale',)) # error: equiv & constrained
    #StoreEquivalence('0:12:Scale',('0:99:Scale',)) # error: equiv & constrained
    varyList = ['2::atomx:3',
                '2::C(10,6,1)', '1::C(10,6,1)',
                '2::atomy:3', '2::atomz:3',
                '0:12:Scale', '0:11:Scale', '0:14:Scale', '0:13:Scale', '0:0:Scale']
#    varyList = ['0::A0', '0::AUiso:0', '0::Afrac:1', '0::Afrac:2', '0::Afrac:3', '0::Afrac:4',
#       '0::dAx:5', '0::dAy:5', '0::dAz:5', '0::AUiso:5', ':0:Back;0', ':0:Back;1', ':0:Back;2', ':0:Back;3', 
#       ':0:Back;4', ':0:Back;5', ':0:Back;6', ':0:Back;7', ':0:Back;8', ':0:Back;9', ':0:Back;10', ':0:Back;11'
#       :0:U', ':0:V', ':0:W', ':0:X', ':0:Y', ':0:Scale', ':0:DisplaceX', ':0:DisplaceY']
#    constrDict = [
#        {'0::Afrac:4': 24.0, '0::Afrac:1': 16.0, '0::Afrac:3': 24.0, '0::Afrac:2': 16.0},
#        {'0::Afrac:1': 1.0, '0::Afrac:2': 1.0},
#        {'0::Afrac:4': 1.0, '0::Afrac:3': 1.0}]
#    fixedList = ['40.0', '1.0', '1.0']

    msg,warning,groups,parmlist = GenerateConstraints(varyList,constrDict,fixedList,parmdict)
    print (VarRemapShow(varyList))
    parmdict.update( {
        '0:12:Scale': 1.0, '0:11:Scale': 1.0, '0:14:Scale': 1.0, '0:13:Scale': 1.0, '0:0:Scale': 2.0,
        '0:0:eA': 0.0,
        '2::C(10,6,1)': 0.2, '1::C(10,6,1)': 0.3,
        '1::C(10,0,1)': 0.2, '2::C(10,0,1)': 0.3,
        '1::AUiso:0': 0.02, '0::AUiso:0': 0.03,
        '0::A0': 0.0,
        '2::atomx:3':0.23,'2::atomy:3':-.23, '2::atomz:3':-0.11,
        })
    print ('parmdict start',parmdict)
    print ('varylist start',varyList)
    before = parmdict.copy()
    Map2Dict(parmdict,varyList)
    print ('parmdict before and after Map2Dict')
    print ('  key / before / after')
    for key in sorted(list(parmdict.keys())):
        print ('  '+key,'\t',before.get(key),'\t',parmdict[key])
    print ('varylist after',varyList)
    before = parmdict.copy()
    Dict2Map(parmdict)
    print ('after Dict2Map')
    print ('  key / before / after')
    for key in sorted(list(parmdict.keys())):
        print ('  '+key,'\t',before.get(key),'\t',parmdict[key])
#    dMdv = len(varylist)*[0]
#    deriv = {}
#    for i,v in enumerate(parmdict.keys()): deriv[v]=i
#    Dict2Deriv(varylist,deriv,dMdv)
