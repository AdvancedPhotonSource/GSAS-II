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

Types of constraints
--------------------

There are four ways to specify constraints, as listed below. 
Note that parameters are initially stored in the 
main section of the GSAS-II data tree under heading ``Constraints``. 
This dict has four keys, 'Hist', 'HAP', 'Global', and 'Phase', 
each containing a list of constraints. An additional set of constraints 
are generated for each phase based on symmetry considerations by calling 
:func:`GSASIIstrIO.GetPhaseData`. 

Note that in the constraints, as stored in the GSAS-II data tree, parameters 
are stored as :class:`GSASIIobj.G2VarObj` objects, as these objects allow for 
changes in numbering of phases, histograms and atoms. 
When they are interpreted (in :func:`GSASIIstrIO.ProcessConstraints`), 
references 
to numbered objects are resolved using the appropriate random ids and the
parameter object is converted to a string of form ``ph:hst:VARNAM:at``.

Alternate parameters (New Var)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameter redefinition ("New Var" constraints) 
is done by creating an expression that relates several 
parameters:

::

   Mx1 * Px + My1 * Py +...
   Mx2 * Px + Mz2 * Pz + ...

where Pj is a GSAS-II parameter name and Mjk is a constant (float) multiplier. 
Alternately, multipliers Mjk can contain a formula (str) that will be evaluated prior 
to the start of the refinement. In a formula, GSAS-II parameters will be replaced by the 
value of the parameter before the formula is evaluated, so ``'np.cos(0::Ax:2)'`` is a valid 
multiplier. At present, only phase (atom/cell) parameters are available for use in 
a formula, but this can be expanded if needed. 

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

These "New Var" constraints are stored as described for type "f" in the 
:ref:`constraint definitions table <Constraint_definitions_table>`.

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

A "Const" constraint is stored as described for type "c" in the 
:ref:`constraint definitions table <Constraint_definitions_table>`.

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
parameter before the formula is evaluate, so ``'np.cos(0::Ax:2)'`` is a valid multiplier. 
At present, only phase (atom/cell) parameters are available for use in 
a formula, but this can be expanded if needed. 
Note that 
the latter constraint expression is conceptually identical to 
defining constraint equations. In practice, however, 
equivalenced parameters are processed in a 
different and more direct manner than constraint equations. The previous 
set of equalities could also be written in this way as a set of constraint 
equations: ::

  C1 * P1 - C2 * P2 = 0
  C1 * P1 - C3 * P3 = 0
  ...


The first parameter (P1 above) 
is considered the independent variable 
and the remaining parameters are dependent variables. The dependent variables 
are set from the independent variable. 
An example of how this may be used woul be if, for example, 
a material has a number of O atoms, all in fairly similar bonding environments
and the diffraction data are sparse, one my reduce the complexity of the model
by defining Uiso for the first O atoms to be identical to the remaining atoms. 
The results of this refinement will be simpler to understand than if a set of
constraint equations is used because the refined parameter will be the independent variable, which will be as ph::Uiso:n, corresponding to the first O atom. 

A parameter can be used in multiple 
equivalences as independent variable,
but if parameter is used as both a dependent and independent variable or
a parameter is used in equivalences and in "New Var" or "Const" constraints, 
this create conflicts that cannot be resolved within the equivalences implementation
but can be handled as constraint equations. 
The equivalences that violate this are discovered in :func:`CheckEquivalences`
and then :func:`MoveConfEquiv` is used to change these equivalences to "Const" equations.

Equivalenced parameters ("EQUIV" constraints), when defined by users, 
are stored as described for type "e" in the 
:ref:`constraint definitions table <Constraint_definitions_table>`.
Other equvalences are generated by symmetry prior 
to display or refinement in :func:`GSASIIstrIO.GetPhaseData`.
These are not stored. 

Fixed parameters (Hold)
^^^^^^^^^^^^^^^^^^^^^^^^

When parameters are refined where a single refinement flag determines that several variables 
are refined at the same time (examples are: cell parameters, atom positions, anisotropic
displacement parameters, magnetic moments,...) it can be useful to specify that a 
specific parameter should not be varied. These will most commonly be generated due to symmetry, 
but under specific conditions, there may be other good reasons to constrain a parameter. 

A "Hold" constraint is stored as described for type "h" in the 
:ref:`constraint definitions table <Constraint_definitions_table>`.

.. _Constraints_processing:

Constraint Processing
---------------------

When constraints will be used or edited, they are processed using a series of
calls:

* First all of the stored constraints are appended into a single list. They are
  initially stored in separate lists only to improve their creation and display
  in the GUI.

* Then :func:`InitVars` is used to initialize the global variables in 
  this module (:mod:`GSASIImapvars`).

* Then :func:`GSASIIstrIO.ProcessConstraints` is used to initially process the 
  constraints, as described below. 

* Symmetry-generated equivalences are then created in 
  :func:`GSASIIstrIO.GetPhaseData`, which also calls 
  :func:`GSASIIstrIO.cellVary` and for Pawley refinements 
  :func:`GSASIIstrIO.GetPawleyConstr`. These are entered directly into this
  module's globals using :func:`StoreEquivalence`.
* Constraints/equivalences are then checked for possible conflicts with
  :func:`CheckConstraints`, this requires grouping the constraints, 
  as described below.
* In refinements, :func:`GenerateConstraints` is then called to 
  create the constraints that will be used, see below for 
* For debugging constraints, :func:`VarRemapShow` can be called after 
  :func:`GenerateConstraints` to display the generated constraints. 

Constraint Reorganization (:func:`~GSASIIstrIO.ProcessConstraints`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`GSASIIstrIO.ProcessConstraints` is used to initially process the 
constraints. This does these things: 

1. The "Hold", "Const" and "New Var" expressions are split between two paired lists, 
   :data:`constDictList` and :data:`fixedList` which are set:
 
   * For "Hold" entries a dict with a single entry is placed in constDictList where 
     the key is the parameter name (associated value is 0.0) and fixedList gets a 
     value of 0.0.
   * For "Const" entries, a dict with multiple entries is placed in constDictList where 
     the key is the parameter name and the value is the multiplier for the parameter, 
     while fixedList gets a string value corresponding to the constant value for 
     the expression. 
   * For "New Var" entries, a dict with multiple entries is placed in constDictList 
     where the key is the parameter name and the value is the multiplier 
     for the parameter; an additional key "_vary" is given the value of True or False
     depending on the refinement flag setting. The corresponding entry in 
     fixedList is None

   The output from this will have this form where the first entry is a "Const", the 
   second is a "New Var" and the final is a "Hold". 

  .. code-block:: python

    constrDict = [
         {'0:12:Scale': 2.0, '0:14:Scale': 4.0, '0:13:Scale': 3.0, '0:0:Scale': 0.5},
         {'2::C(10,6,1)': 1.0, '1::C(10,6,1)': 1.0, '_vary':True},
         {'0::A0': 0.0}]
    fixedList = ['5.0', None, '0']


2. Equivalences are stored using :func:`StoreEquivalence`
   into this module's globals (:data:`~GSASIImapvars.arrayList`, 
   :data:`~GSASIImapvars.invarrayList`,    :data:`~GSASIImapvars.indParmList`, 
   :data:`~GSASIImapvars.dependentParmList` and  :data:`~GSASIImapvars.symGenList`)


Parameter Grouping (:func:`GenerateConstraints`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions :func:`CheckConstraints` and :func:`GenerateConstraints` are used to
process the parameter equivalences and constraint lists created in 
:func:`~GSASIIstrIO.ProcessConstraints`. The former is used to generate 
error messages and the latter to generate the internal information used to 
apply the constraints. 

Initially, in both a list of parameters that are fixed and those used in 
constraint relations are tabulated in :func:`CheckEquivalences`. The equivalence
relations are the scanned for the following potential problems:

1. a parameter is used as a dependent variable in more than one 
   equivalence relation
2. a parameter is fixed and used in in an equivalence relation either 
   as a dependent or independent variable
3. a parameter is used as a dependent variable in one equivalence relation and 
   as a independent variable in another
4. a parameter is used in in an equivalence relation (either 
   as a dependent or independent variable) and is used in a constraint 
   expression
5. a parameter is not defined in a particular refinement, but is used in an
   equivalence relation
6. a parameter uses a wildcard for the histogram number (sequential refinements)

Cases 1 & 2 above cannot be corrected, and result in errors. Cases 3 & 4 are 
potentially corrected with :func:`MoveConfEquiv`, as described below. 
Case 5 causes the equivalence to 
be dropped. Case 6 causes the current histogram number to be substituted for
the wildcard. 

For cases 3 & 4, :func:`MoveConfEquiv` is used to change these equivalences 
into "Const" equations. This can potentially mean that additional 
equivalences will be problematic, so if there are changes made by 
:func:`MoveConfEquiv`, :func:`CheckEquivalences` is repeated. If any problem 
cases are noted, the refinement cannot be performed. 

Constraint expressions ("Const" and "New Var") are sorted into 
groups so that each group contains the minimum number of entries that 
ensures each parameter is referenced in only one group in 
:func:`GroupConstraints`. This is done by scanning the 
list of dicts in :data:`constDictList` one by one and making a list 
of parameters used in that constraint expression. Any expression that contains
a parameter in is in that list is added to the current group and those 
parameters are added to this list of parameters. The list of ungrouped 
expressions is then scanned again until no more expressions are added to the 
current group. This process is repeated until every expression has been 
placed in a group. Function :func:`GroupConstraints` returns two lists of lists.
The first has, for each group, a list of the indices in :data:`constDictList` 
that comprise the group (there can be only one). The second list contains, 
for each group, the unique parameter names in that group. 

Each constraint group is then processed. First, wildcard parameters are
renamed (in a sequential refinement). Any fixed parameters that are used 
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
:func:`_FillArray` with parameter ``FillDiagonals=False`` (the default). Any 
unspecified rows are left as all zero. The first Nc rows in the
array are then coverted to row-echelon form using :func:`_RowEchelon`. This 
will create an Exception if any two rows are linearly dependent (which means 
that no matter what values are used for the remaining rows, that the matrix 
will be singular). :func:`_FillArray` is then called with parameter 
``FillDiagonals=True``, which again creates a square matrix but where 
unspecified rows are zero except for the diagonal elements. The  
`Gram-Schmidt process <http://en.wikipedia.org/wiki/Gram-Schmidt>`_, 
implemented  in :func:`GramSchmidtOrtho`, is used to find orthonormal unit 
vectors for the remaining Np-Nc rows of the matrix. This will fail with 
a ``ConstraintException`` if this is not possible (singular matrix) or 
the result is placed in :data:`constrArr` as a numpy array. 

Rows in the matrix corresponding to "New Var" constraints and those that 
were generated by the Gram-Schmidt process are provided with parameter names
(this can be specified if a "New Var" entry by using a ``"_name"`` element
in the constraint dict, but at present this is not implemented.) Names are 
generated using :data:`paramPrefix` which is set to ``"::constr"``, plus a 
number to make the new parameter name unique. Global dict :data:`genVarLookup`
provides a lookup table, where the names of the parameters related to this new
parameter can be looked up easily. 

Finally the parameters used as input to the constraint are placed in 
this module's globals
:data:`~GSASIImapvars.dependentParmList` and the constraint matrix is 
placed in into  :data:`~GSASIImapvars.arrayList`. This can be used to compute 
the initial values for "New Var" parameters. The inverse of the 
constraint matrix is placed in :data:`~GSASIImapvars.invarrayList` and a list 
of the "New Var" parameters and a list of the fixed values (as str's) 
is placed in :data:`~GSASIImapvars.indParmList`. A lookup table for
fixed values as floats is placed in :data:`~GSASIImapvars.fixedDict`. 
Finally the appropriate entry in :data:`~GSASIImapvars.symGenList` is set to 
False to indicate that this is not a symmetry generated constraint. 


*Externally-Accessible Routines*
---------------------------------

To define a set of constrained and unconstrained relations, one
defines a list of dictionary defining constraint parameters and their
values, a list of fixed values for each constraint and a list of
parameters to be varied. In addition, one uses
:func:`StoreEquivalence` to define parameters that are equivalent. 
Use :func:`EvaluateMultipliers` to convert formula-based constraint/equivalence
multipliers to numbers and then 
use :func:`CheckConstraints` to check that the input is
internally consistent and finally :func:`GroupConstraints` and
:func:`GenerateConstraints` to generate the internally used
tables. Routine :func:`Map2Dict` is used to initialize the parameter
dictionary and routine :func:`Dict2Map`, :func:`Dict2Deriv`, and
:func:`ComputeDepESD` are used to apply constraints. Routine
:func:`VarRemapShow` is used to print out the constraint information,
as stored by :func:`GenerateConstraints`. Further information on each routine
is below: 

:func:`InitVars`
  This is optionally used to clear out all defined previously defined constraint information
  
:func:`StoreEquivalence`
  To implement parameter redefinition, one calls StoreEquivalence. This should be called for every set of
  equivalence relationships. There is no harm in using StoreEquivalence with the same independent variable:

  .. code-block:: python

       StoreEquivalence('x',('y',))
       StoreEquivalence('x',('z',))

  or equivalently 

  .. code-block:: python

       StoreEquivalence('x',('y','z'))

  The latter will run more efficiently. Note that mixing independent 
  and dependent variables would require assignments, such as 

  .. code-block:: python

        StoreEquivalence('x',('y',))
        StoreEquivalence('y',('z',))

  would require that equivalences be applied in a particular order and 
  thus is implemented as a constraint equation rather than an equivalence. 
        
  Use StoreEquivalence before calling GenerateConstraints or CheckConstraints

:func:`CheckConstraints`
   check that input in internally consistent

:func:`GenerateConstraints` 
   generate the internally used tables from constraints and equivalences

:func:`EvaluateMultipliers`
   Convert any string-specified (formula-based) multipliers to numbers. Call this before 
   using :func:`CheckConstraints` or :func:`GenerateConstraints`. 
   At present, the code may pass only the dict for phase (atom/cell) parameters, 
   but this could be expanded if needed. 

:func:`Map2Dict`
   To determine values for the parameters created in this module, one
   calls Map2Dict. This will not apply contraints.

:func:`Dict2Map`
   To take values from the new independent parameters and constraints,
   one calls Dict2Map and set the parameter array, thus appling contraints.

:func:`Dict2Deriv`
   Use Dict2Deriv to determine derivatives on independent parameters
   from those on dependent ones.

:func:`ComputeDepESD`      
   Use ComputeDepESD to compute uncertainties on dependent variables.

:func:`VarRemapShow`
   To show a summary of the parameter remapping, one calls VarRemapShow.

*Global Variables*
------------------

dependentParmList:
   a list containing group of lists of
   parameters used in the group. Note that parameters listed in
   dependentParmList should not be refined as they will not affect
   the model

indParmList:
     a list containing groups of Independent parameters defined in
     each group. This contains both parameters used in parameter
     redefinitions as well as names of generated new parameters.

arrayList:
     a list containing group of relationship matrices to relate
     parameters in dependentParmList to those in indParmList. Unlikely
     to be used externally.

invarrayList:
     a list containing group of relationship matrices to relate
     parameters in indParmList to those in dependentParmList. Unlikely
     to be used externally.

fixedVarList:
     a list of parameters that have been 'fixed'
     by defining them as equal to a constant (::var: = 0). Note that
     the constant value is ignored at present. These parameters are
     later removed from varyList which prevents them from being refined. 
     Unlikely to be used externally.

fixedDict:
     a dictionary containing the fixed values corresponding
     to parameter equations.  The dict key is an ascii string, but the
     dict value is a float.  Unlikely to be used externally.

symGenList:
     a list of boolean values that will be True to indicate that a constraint
     (only equivalences) is generated by symmetry (or Pawley overlap)
     
problemVars:
     a list containing parameters that show up in constraints producing errors



*Routines/variables*
---------------------

Note that parameter names in GSAS-II are strings of form ``<ph#>:<hst#>:<nam>`` or ``<ph#>::<nam>:<at#>``. 
"""

from __future__ import division, print_function
import numpy as np
import sys
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
# data used for constraints; 
debug = False # turns on printing as constraint input is processed

# note that constraints are stored listed by contraint groups,
# where each constraint
# group contains those parameters that must be handled together
dependentParmList = []
'''a list of lists where each item contains a list of parameters in each constraint group. 
note that parameters listed in dependentParmList should not be refined directly.'''
indParmList = [] # a list of names for the new parameters
'''a list of lists where each item contains a list for each constraint group with 
fixed values for constraint equations and names of generated (New Var) parameters.
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
fixedDict = {}
'''A dict lookup-table containing the fixed values corresponding 
to defined parameter equations. Note the key is the original ascii string
and the value in the dict is a float.
'''
fixedVarList = []
'''List of parameters that should not be refined.
'''
symGenList = []
'''A list of flags that if True indicates a constraint was generated by symmetry
'''
problemVars = []
'''a list of parameters causing errors
'''
dependentVars = []
'A list of dependent variables, taken from (:data:`dependentParmList`).'
independentVars = []
'A list of dependent variables, taken from (:data:`indParmList`).'
genVarLookup = {}
'provides a list of parameters that are related to each generated parameter'
paramPrefix = "::constr"
'A prefix for generated parameter names'
consNum = 0
'The number to be assigned to the next constraint to be created'

class ConstraintException(Exception):
    '''Defines an Exception that is used when an exception is raised processing constraints
    '''
    pass

def InitVars():
    '''Initializes all constraint information'''
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict,consNum,symGenList
    dependentParmList = [] # contains a list of parameters in each group
    arrayList = [] # a list of of relationship matrices 
    invarrayList = [] # a list of inverse relationship matrices 
    indParmList = [] # a list of names for the new parameters
    fixedDict = {} # a dictionary containing the fixed values corresponding to defined parameter equations
    symGenList = [] # Flag if constraint is generated by symmetry
    consNum = 0 # number of the next constraint to be created
    global genVarLookup
    genVarLookup = {}

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
    """divide the constraints into groups that share no parameters.

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
    for i,consi in enumerate(constrDict):
        if i in assignedlist: continue # already in a group, skip
        # starting a new group
        grouplist = [i,]
        assignedlist.append(i)
        groupset = set(VarKeys(consi))
        changes = True # always loop at least once
        while(changes): # loop until we can't find anything to add to the current group
            changes = False # but don't loop again unless we find something
            for j,consj in enumerate(constrDict):
                if j in assignedlist: continue # already in a group, skip
                if len(set(VarKeys(consj)) & groupset) > 0: # true if this needs to be added
                    changes = True
                    grouplist.append(j)
                    assignedlist.append(j)
                    groupset = groupset | set(VarKeys(consj))
        group = sorted(grouplist)
        varlist = sorted(list(groupset))
        groups.append(group)
        ParmList.append(varlist)
    return groups,ParmList

def CheckConstraints(varyList,constrDict,fixedList):
    '''Takes a list of relationship entries comprising a group of
    constraints and checks for inconsistencies such as conflicts in
    parameter/variable definitions and or inconsistently varied parameters.

    :param list varyList: a list of parameters names that will be varied

    :param dict constrDict: a list of dicts defining relationships/constraints
      (as created in :func:`GSASIIstrIO.ProcessConstraints` and
      documented in :func:`GroupConstraints`)

    :param list fixedList: a list of values specifying a fixed value for each
      dict in constrDict. Values are either strings that can be converted to
      floats or ``None`` if the constraint defines a new parameter rather
      than a constant.

    :returns: two strings: 

      * the first lists conflicts internal to the specified constraints
      * the second lists conflicts where the varyList specifies some
        parameters in a constraint, but not all
        
      If there are no errors, both strings will be empty
    '''
    import re
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    global problemVars
    # Process the equivalences
    #    If there are conflicting parameters, move them into constraints. This
    #    may create new conflicts, requiring movement of other parameters, so
    #    loop until there are no more changes to make.
    parmsChanged = True
    while parmsChanged:
        parmsChanged = 0
        errmsg,warnmsg,fixVlist,dropVarList,translateTable = CheckEquivalences(
            constrDict,varyList)
        #print('debug: using MoveConfEquiv to address =',errmsg)
        if problemVars: parmsChanged = MoveConfEquiv(constrDict,fixedList)
#    GSASIIpath.IPyBreak()

    groups,parmlist = GroupConstraints(constrDict)
    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create a warning message.
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1:   # process fixed (held) variables
            var = varlist[0]
            if var not in fixedVarList:
                fixedVarList.append(var)
            continue
        for rel in group:
            varied = 0
            notvaried = ''
            for var in constrDict[rel]:
                if var.startswith('_'): continue
                if not re.match('[0-9]*:[0-9\*]*:',var):
                    warnmsg += "\nParameter "+str(var)+" does not begin with a ':'"
                if var in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += var
                if var in fixVlist:
                    errmsg += '\nParameter '+var+" is Fixed and used in a constraint:\n\t"
                    if var not in problemVars: problemVars.append(var)
                    errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
            if varied > 0 and varied != len(VarKeys(constrDict[rel])):
                warnmsg += "\nNot all parameters refined in constraint:\n\t"
                warnmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                warnmsg += '\nNot refined: ' + notvaried + '\n'
    if errmsg or warnmsg:
        return errmsg,warnmsg

    # now look for process each group and create the relations that are needed to form
    # non-singular square matrix
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue # a constraint group with a single parameter can be ignored
        if len(varlist) < len(group): # too many relationships -- no can do
            errmsg += "\nOver-constrained input. "
            errmsg += "There are more constraints " + str(len(group))
            errmsg += "\n\tthan parameters " + str(len(varlist)) + "\n"
            for rel in group:
                errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
                continue
        try:
            multarr = _FillArray(group,constrDict,varlist)
            _RowEchelon(len(group),multarr,varlist)
        except:
            errmsg += "\nSingular input. "
            errmsg += "There are internal inconsistencies in these constraints\n"
            for rel in group:
                errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
            continue
        try:
            multarr = _FillArray(group,constrDict,varlist,FillDiagonals=True)
            GramSchmidtOrtho(multarr,len(group))
        except:
            errmsg += "\nUnexpected singularity with constraints (in Gram-Schmidt)\n"
            for rel in group:
                errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
            continue
        mapvar = []
        group = group[:]
        # scan through all generated and input parameters
        # Check again for inconsistent parameter use
        # for new parameters -- where varied and unvaried parameters get grouped
        # together. I don't think this can happen when not flagged before, but
        # it does not hurt to check again. 
        for i in range(len(varlist)):
            varied = 0
            notvaried = ''
            if len(group) > 0:
                rel = group.pop(0)
                fixedval = fixedList[rel]
                for var in VarKeys(constrDict[rel]):
                    if var in varyList:
                        varied += 1
                    else:
                        if notvaried: notvaried += ', '
                        notvaried += var
            else:
                fixedval = None
            if fixedval is None:
                varname = paramPrefix + str(consNum) # assign a name to a parameter
                mapvar.append(varname)
                consNum += 1
            else:
                mapvar.append(fixedval)
            if varied > 0 and notvaried != '':
                warnmsg += "\nNot all parameters refined in generated constraint"
                warnmsg += '\nPlease report this unexpected error\n'
                for rel in group:
                    warnmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                    warnmsg += "\n"
                warnmsg += '\n\tNot refined: ' + notvaried + '\n'
        try:
            np.linalg.inv(multarr)            
        except:
            errmsg += "\nSingular input. "
            errmsg += "The following constraints are not "
            errmsg += "linearly independent\n\tor do not "
            errmsg += "allow for generation of a non-singular set\n"
            errmsg += 'Please report this unexpected error\n'
            for rel in group:
                errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
    _setVarLists([])
    return errmsg,warnmsg

def GenerateConstraints(varyList,constrDict,fixedList,parmDict=None,SeqHist=None):
    '''Takes a list of relationship entries comprising a group of
    constraints and builds the relationship lists and their inverse
    and stores them in global parameters Also checks for internal
    conflicts or inconsistencies in parameter/variable definitions.

    :param list varyList: a list of parameters names (strings of form
      ``<ph>:<hst>:<nam>``) that will be varied. Note that this is changed here. 
    
    :param dict constrDict: a list of dicts defining relationships/constraints
      (as defined in :func:`GroupConstraints`)

    :param list fixedList: a list of values specifying a fixed value for each
      dict in constrDict. Values are either strings that can be converted to
      floats, float values or None if the constraint defines a new parameter.
      
    :param dict parmDict: a dict containing all parameters defined in current
      refinement.

    :param int SeqHist: number of current histogram, when used in a sequential
      refinement. None (default) otherwise. Wildcard parameter names are
      set to the current histogram, when found if not None.
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    global genVarLookup
    msg = ''
    shortmsg = ''
    changed = False
    
    # Process the equivalences
    #    If there are conflicting parameters, move them into constraints. This
    #    may create new conflicts, requiring movement of other parameters, so
    #    loop until there are no more changes to make.
    parmsChanged = True
    while parmsChanged:
        parmsChanged = 0
        errmsg,warnmsg,fixVlist,dropVarList,translateTable = CheckEquivalences(
            constrDict,varyList,parmDict,SeqHist)
        if problemVars:
            parmsChanged = MoveConfEquiv(constrDict,fixedList)
            changed = True
    if errmsg:
        msg = errmsg
    if warnmsg:
        if msg: msg += '\n'
        msg += warnmsg

    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create an error message. 
    groups,parmlist = GroupConstraints(constrDict)
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1:   # process fixed (held) variables
            var = varlist[0]
            if var not in fixedVarList:
                fixedVarList.append(var)
            continue
        for rel in group:
            varied = 0
            notvaried = ''
            unused = 0
            notused = ''
            for var in constrDict[rel]:
                if var.startswith('_'): continue
                if var.split(':')[1] == '*' and SeqHist is not None:
                    # convert wildcard var to reference current histogram; save translation in table
                    sv = var.split(':')
                    sv[1] = str(SeqHist)
                    translateTable[var] = ':'.join(sv)
                    var = translateTable[var]
                if parmDict is not None and var not in parmDict:
                    unused += 1
                    if notvaried: notused += ', '
                    notused += var
                if var in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += var
                if var in fixedVarList:
                    msg += '\nError: parameter '+var+" is Fixed and used in a constraint:\n\t"
                    msg += _FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
            #if unused > 0:# and unused != len(VarKeys(constrDict[rel])):
            if unused > 0 and unused != len(VarKeys(constrDict[rel])):
                #msg += "\nSome (but not all) parameters in constraint are not defined:\n\t"
                #msg += _FormatConstraint(constrDict[rel],fixedList[rel])
                #msg += '\nNot used: ' + notused + '\n'
                shortmsg += notused+" not used in constraint "+_FormatConstraint(constrDict[rel],fixedList[rel])
            elif varied > 0 and varied != len(VarKeys(constrDict[rel])):
                #msg += "\nNot all parameters refined in constraint:\n\t"
                #msg += _FormatConstraint(constrDict[rel],fixedList[rel])
                #msg += '\nNot refined: ' + notvaried + '\n'
                shortmsg += notvaried+" not varied in constraint "+_FormatConstraint(constrDict[rel],fixedList[rel])
    # if there were errors found, go no farther
    if shortmsg and SeqHist is not None:
        if msg:
            print (' *** ERROR in constraint definitions! ***')
            print (msg)
            raise ConstraintException
        print ('*** Sequential refinement: ignoring constraint definition(s): ***')
        print (shortmsg)
        msg = ''
    elif shortmsg:
        msg += shortmsg
    if msg:
        print (' *** ERROR in constraint definitions! ***')
        print (msg)
        raise ConstraintException
                
    # now process each group and create the relations that are needed to form
    # a non-singular square matrix
    # If all are varied and this is a constraint equation, then set VaryFree flag
    # so that the newly created relationships will be varied
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        # for constraints, if all included parameters are refined,
        # set the VaryFree flag, and remaining degrees of freedom will be
        # varied (since consistency was checked, if any one parameter is
        # refined, then assume that all are)
        varsList = [] # make a list of all the referenced parameters as well
        VaryFree = False
        for rel in group:
            varied = 0
            unused = 0
            for var in VarKeys(constrDict[rel]):
                var = translateTable.get(var,var) # replace wildcards
                if parmDict is not None and var not in parmDict:
                    unused += 1                    
                if var not in varsList: varsList.append(var)
                if var in varyList: varied += 1
            if fixedList[rel] is not None and varied > 0:
                VaryFree = True
        if len(varlist) < len(group): # too many relationships -- no can do
            msg = 'too many relationships'
            break
        # Since we checked before, if any parameters are unused, then all must be. 
        # If so, this set of relationships can be ignored
        if unused:
            if debug: print('Constraint ignored (all parameters undefined)')
            if debug: print ('    '+_FormatConstraint(constrDict[rel],fixedList[rel]))
            continue
        # fill in additional degrees of freedom
        try:
            arr = _FillArray(group,constrDict,varlist)
            _RowEchelon(len(group),arr,varlist)
            constrArr = _FillArray(group,constrDict,varlist,FillDiagonals=True)
            GramSchmidtOrtho(constrArr,len(group))
        except:
            msg = 'Singular relationships'
            break
        mapvar = []
        group = group[:]
        # scan through all generated and input relationships, we need to add to the varied list
        # all the new parameters where VaryFree has been set or where a New Var is varied.
        #
        # If a group does not contain any fixed values (constraint equations)
        # and nothing in the group is varied, drop this group, so that the 
        # dependent parameters can be refined individually.
        unused = True
        for i in range(len(varlist)):
            if len(group) > 0: # get the original equation reference
                rel = group.pop(0)
                fixedval = fixedList[rel]
                varyflag = constrDict[rel].get('_vary',False)
                varname = constrDict[rel].get('_name','')
            else: # this relationship has been generated
                varyflag = False
                varname = ''
                fixedval = None
            if fixedval is None: # this is a new parameter, not a constraint
                if not varname:
                    varname = paramPrefix + str(consNum) # no assigned name, create one
                    consNum += 1
                mapvar.append(varname)
                genVarLookup[varname] = varlist # save list of parameters related to this new var
                # vary the new relationship if it is a degree of freedom in
                # a set of contraint equations or if a New Var is flagged to be varied.
                if VaryFree or varyflag: 
                    unused = False
                    varyList.append(varname)
                    # fix (prevent varying) of all the parameters inside the constraint group
                    # (dependent vars)
                    for var in varsList:
                        if var in varyList: varyList.remove(var)
            else:
                unused = False
                mapvar.append(fixedval)
        if unused: # skip over constraints that don't matter (w/o fixed value or any refined parameters)
            if debug: print('Constraint ignored (all parameters unrefined)')
            if debug: print ('   '+_FormatConstraint(constrDict[rel],fixedList[rel]))
            continue 
        dependentParmList.append([translateTable.get(var,var) for var in varlist])
        arrayList.append(constrArr)
        invarrayList.append(np.linalg.inv(constrArr))
        indParmList.append(mapvar)
        symGenList.append(False)
    if msg:
        print (' *** ERROR in constraint definitions! ***')
        print (msg)
        print (VarRemapShow(varyList))
        raise ConstraintException
    # setup dictionary containing the fixed values
    global fixedDict 
    # key is original ascii string, value is float
    for fixedval in fixedList:
        if fixedval:
            fixedDict[fixedval] = float(fixedval)
    _setVarLists(dropVarList)
    if changed:
        print(60*'=')
        print('Constraints were reclassified to avoid conflicts, as below:')
        print (VarRemapShow(varyList,True))
        print(60*'=')
    return groups,parmlist # saved for sequential fits
    
def _setVarLists(dropVarList):
    '''Make list of dependent and independent variables (after dropping unused vars in dropVarList)
    '''
    global dependentParmList,indParmList
    global dependentVars
    global independentVars
    dependentVars = []
    independentVars = []
    for varlist,mapvars in zip(dependentParmList,indParmList):  # process all constraints
        for mv in mapvars:
            if mv in dropVarList: continue
            if mv not in independentVars: independentVars.append(mv)
        for mv in varlist:
            if mv in dropVarList: continue
            if mv not in dependentVars: dependentVars.append(mv)
    if debug: # on debug, show what is parsed & generated, semi-readable
        print (50*'-')
        #print (VarRemapShow(varyList))
        #print ('Varied: ',varyList)
        print ('Not Varied: ',fixedVarList)

# def CheckEquivalences(constrDict,varyList):
#     global dependentParmList,arrayList,invarrayList,indParmList,consNum
#     global problemVars
#     warnmsg = ''
#     errmsg = ''
#     problemVars = []
#     # process fixed variables (holds)
#     fixVlist = [] # list of Fixed vars 
#     constrVars = [] # list of vars in constraint expressions
#     for cdict in constrDict:
#         # N.B. No "_" names in holds
#         if len(cdict) == 1:
#             fixVlist.append(list(cdict.keys())[0])
#         else:
#             constrVars += cdict.keys() # this will include _vary (not a problem)
#     # process equivalences: make a list of dependent and independent vars
#     #    and check for repeated uses (repetition of a parameter as an
#     #    independent var is OK)
#     indepVarList = []
#     depVarList = []
#     multdepVarList = []
#     for varlist,mapvars,multarr,invmultarr in zip(
#         dependentParmList,indParmList,arrayList,invarrayList):
#         if multarr is None: # an equivalence
#             zeromult = False
#             for mv in mapvars:
#                 varied = 0
#                 notvaried = ''
#                 if mv in varyList:
#                     varied += 1
#                 else:
#                     if notvaried: notvaried += ', '
#                     notvaried += mv
#                 if mv not in indepVarList: indepVarList.append(mv)
#                 for v,m in zip(varlist,invmultarr):
#                     if v in indepVarList:
#                         errmsg += '\nVariable '+v+' is used to set values in a constraint before its value is set in another constraint\n'
#                         if v not in problemVars: problemVars.append(v)
#                     if m == 0: zeromult = True
#                     if v in varyList:
#                         varied += 1
#                     else:
#                         if notvaried: notvaried += ', '
#                         notvaried += v
#                     if v in depVarList:
#                         multdepVarList.append(v)
#                     else:
#                         depVarList.append(v)
#             if varied > 0 and varied != len(varlist)+1:
#                 warnmsg += "\nNot all variables refined in equivalence:\n\t"
#                 s = ""
#                 for v in varlist:
#                     if s != "": s+= " & "
#                     s += str(v)            
#                 warnmsg += str(mv) + " => " + s
#                 warnmsg += '\nNot refined: ' + notvaried + '\n'
#             if zeromult:
#                 errmsg += "\nZero multiplier is invalid in equivalence:\n\t"
#                 s = ""
#                 for v in varlist:
#                     if s != "": s+= " & "
#                     s += str(v)            
#                 errmsg += str(mv) + " => " + s + '\n'
#     # check for errors:
#     if len(multdepVarList) > 0:
#         errmsg += "\nThe following parameters(s) are used in conflicting Equivalence relations as dependent variables:\n"
#         s = ''
#         for var in sorted(set(multdepVarList)):
#             if v not in problemVars: problemVars.append(v)
#             if s != "": s+= ", "
#             s += str(var)            
#         errmsg += '\t'+ s + '\n'
#     equivVarList = list(set(indepVarList).union(set(depVarList)))
#     if debug: print ('equivVarList',equivVarList)
#     # check for parameters that are both fixed and in an equivalence (not likely)
#     inboth = set(fixVlist).intersection(set(equivVarList))
#     if len(inboth) > 0:
#         errmsg += "\nThe following parameter(s) are used in both Equivalence and Fixed constraints:\n"
#         s = ''
#         for var in sorted(inboth):
#             if var not in problemVars: problemVars.append(var)
#             if s != "": s+= ", "
#             s += str(var)
#         errmsg += '\t'+ s + '\n'
#     # check for parameters that in both an equivalence and a constraint expression 
#     inboth = set(constrVars).intersection(set(equivVarList))
#     if len(inboth) > 0:
#         errmsg += "\nThe following parameter(s) are used in both Equivalence and Equiv or new var constraints:\n"
#         s = ''
#         for var in sorted(inboth):
#             if var not in problemVars: problemVars.append(var)
#             if s != "": s+= ", "
#             s += str(var)
#         errmsg += '\t'+ s + '\n'
#     return errmsg,warnmsg,fixVlist

def CheckEquivalences(constrDict,varyList,parmDict=None,SeqHist=None):
    '''Process equivalence constraints, looking for conflicts such as 
    where a parameter is used in both an equivalence and a constraint expression
    or where chaining is done (A->B and B->C). 
    When called during refinements, parmDict is defined, and for sequential refinement 
    SeqHist ia also defined.

      * parmDict is used to remove equivalences where a parameter is not present 
        in a refinement
      * SeqHist is used to rename wild-card parameter names in sequential 
        refinements to use the current histogram.
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    global problemVars
    warnmsg = ''
    errmsg = ''
    problemVars = []
    # process fixed parameters (holds)
    fixVlist = [] # list of Fixed vars 
    constrVars = [] # list of vars in constraint expressions
    for cdict in constrDict:
        # N.B. No "_" names in holds
        if len(cdict) == 1:
            fixVlist.append(list(cdict.keys())[0])
        else:
            constrVars += cdict.keys() # this will include _vary (not a problem)
    # process equivalences: make a list of dependent and independent vars
    #    and check for repeated uses (repetition of a parameter as an
    #    independent var is OK)
    indepVarList = []
    depVarList = []
    multdepVarList = []
    dropVarList = []
    translateTable = {} # lookup table for wildcard referenced parameters
    for varlist,mapvars,multarr,invmultarr in zip(
        dependentParmList,indParmList,arrayList,invarrayList):
        if multarr is None: # an equivalence
            zeromult = False
            for i,mv in enumerate(mapvars):
                if mv.split(':')[1] == '*' and SeqHist is not None:
                    # convert wildcard var to reference current histogram; save translation in table
                    sv = mv.split(':')
                    sv[1] = str(SeqHist)
                    mv = translateTable[mv] = ':'.join(sv)
                    mapvars[i] = mv
                varied = 0
                notvaried = ''
                if mv in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += mv
                if parmDict is not None and mv not in parmDict:
                    print ("Dropping equivalence for parameter "+str(mv)+". Not defined in this refinement")
                    if mv not in dropVarList: dropVarList.append(mv)
                if mv not in indepVarList: indepVarList.append(mv)
            for i,(v,m) in enumerate(zip(varlist,invmultarr)):
                if v.split(':')[1] == '*' and SeqHist is not None:
                    # convert wildcard var to reference current histogram; save translation in table
                    sv = v.split(':')
                    sv[1] = str(SeqHist)
                    varlist[i] = v = translateTable[v] = ':'.join(sv)
                if parmDict is not None and v not in parmDict:
                    print ("Dropping equivalence for dep. variable "+str(v)+". Not defined in this refinement")
                    if v not in dropVarList: dropVarList.append(v)
                    continue
                if m == 0: zeromult = True
                if v in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += v
                if v in indepVarList:
                    errmsg += '\nParameter '+v+' is used to set values in a constraint before its value is set in another constraint\n'
                    if v not in problemVars: problemVars.append(v)
                if v in depVarList:
                    multdepVarList.append(v)
                else:
                    depVarList.append(v)
            if varied > 0 and varied != len(varlist)+1:
                warnmsg += "\nNot all parameters refined in equivalence:\n\t"
                s = ""
                for v in varlist:
                    if s != "": s+= " & "
                    s += str(v)            
                warnmsg += str(mv) + " => " + s
                warnmsg += '\nNot refined: ' + notvaried + '\n'
            if zeromult:
                errmsg += "\nZero multiplier is invalid in equivalence:\n\t"
                s = ""
                for v in varlist:
                    if s != "": s+= " & "
                    s += str(v)            
                errmsg += str(mv) + " => " + s + '\n'
    # check for errors:
    if len(multdepVarList) > 0:
        errmsg += "\nThe following parameters(s) are used in conflicting Equivalence relations as dependent variables:\n"
        s = ''
        for var in sorted(set(multdepVarList)):
            if v not in problemVars: problemVars.append(v)
            if s != "": s+= ", "
            s += str(var)            
        errmsg += '\t'+ s + '\n'
    equivVarList = list(set(indepVarList).union(set(depVarList)))
    if debug: print ('equivVarList',equivVarList)
    # check for parameters that are both fixed and in an equivalence (not likely)
    inboth = set(fixVlist).intersection(set(equivVarList))
    if len(inboth) > 0:
        errmsg += "\nThe following parameter(s) are used in both Equivalence and Fixed constraints:\n"
        s = ''
        for var in sorted(inboth):
            if var not in problemVars: problemVars.append(var)
            if s != "": s+= ", "
            s += str(var)
        errmsg += '\t'+ s + '\n'
    # check for parameters that in both an equivalence and a constraint expression 
    inboth = set(constrVars).intersection(set(equivVarList))
    if len(inboth) > 0:
        errmsg += "\nThe following parameter(s) are used in both Equivalence and Equiv or new var constraints:\n"
        s = ''
        for var in sorted(inboth):
            if var not in problemVars: problemVars.append(var)
            if s != "": s+= ", "
            s += str(var)
        errmsg += '\t'+ s + '\n'
    return errmsg,warnmsg,fixVlist,dropVarList,translateTable

def MoveConfEquiv(constrDict,fixedList):
    '''Address conflicts in Equivalence constraints by creating an constraint equation 
    that has the same action as the equivalence and removing the Equivalence
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    global problemVars
    parmsChanged = 0
    for i,(varlist,mapvars) in enumerate(zip(dependentParmList,indParmList)):
        conf = False
        for mv in mapvars:
            if mv in problemVars:
                conf = True
                break
        for v in varlist:
            if v in problemVars:
                conf = True
                break
        if conf:
            parmsChanged += 1
            indvar = indParmList[i][0]
            for dep,mult in zip(dependentParmList[i],invarrayList[i]):
                #print('debug replacing equiv with constraint equation 0=',{indvar:-1.,dep:mult[0]})
                constrDict += [{indvar:-1.,dep:mult[0]}]
                fixedList += ['0.0']
            dependentParmList[i] = None
    if parmsChanged:
        for i in range(len(dependentParmList)-1,-1,-1):
            if dependentParmList[i] is None:
                del dependentParmList[i],indParmList[i],arrayList[i],invarrayList[i]
    return parmsChanged

def StoreEquivalence(independentVar,dependentList,symGen=True):
    '''Takes a list of dependent parameter(s) and stores their
    relationship to a single independent parameter (independentVar).

    Called with user-supplied constraints by :func:`GSASIIstrIO.ProcessConstraints,
    with Pawley constraints from :func:`GSASIIstrIO.GetPawleyConstr`, 
    with Unit Cell constraints from :func:`GSASIIstrIO.cellVary`
    with symmetry-generated atom constraints from :func:`GSASIIstrIO.GetPhaseData`

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
    def SubfromParmDict(s,prmDict):
        for key in prmDict:
            if key in s:
                s = s.replace(key,str(prmDict[key]))
        return eval(s)
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
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    printlist = []
    mapvars = GetIndependentVars()
    for i,name in enumerate(mapvars):
        if name in fixedDict: continue
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

def ComputeDepESD(covMatrix,varyList,parmDict):
    '''Compute uncertainties for dependent parameters from independent ones
    returns a dictionary containing the esd values for dependent parameters
    '''
    sigmaDict = {}
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        #if invmultarr is None: continue # probably not needed
#        try: 
#            valuelist = [parmDict[var] for var in mapvars]
#        except KeyError:
#            continue
        # get the v-covar matrix for independent parameters 
        vcov = np.zeros((len(mapvars),len(mapvars)))
        for i1,name1 in enumerate(mapvars):
            if name1 not in varyList: continue
            iv1 = varyList.index(name1)
            for i2,name2 in enumerate(mapvars):
                if name2 not in varyList: continue
                iv2 = varyList.index(name2)
                vcov[i1][i2] = covMatrix[iv1][iv2]
        # vec is the vector that multiplies each of the independent values
        for v,vec in zip(varlist,invmultarr):
            sigmaDict[v] = np.sqrt(np.inner(vec.T,np.inner(vcov,vec)))
    return sigmaDict

def _FormatConstraint(RelDict,RelVal):
    '''Formats a Constraint or Function for use in a convenient way'''
    linelen = 45
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
        s[-1] += '%.3f*%s '%(m,var)
    if len(s[-1]) > linelen: s.append(' ')
    if RelVal is None:
        s[-1] += ' = New variable'
    else:
        s[-1] += ' = ' + RelVal
    s1 = ''
    for s2 in s:
        if s1 != '': s1 += '\n\t'
        s1 += s2
    return s1

def VarRemapShow(varyList,inputOnly=False):
    '''List out the saved relationships. This should be done after the constraints have been
    defined using :func:`StoreEquivalence`, :func:`GroupConstraints` and :func:`GenerateConstraints`.

    :returns: a string containing the details of the contraint relationships
    '''
    s = ''
    if len(fixedVarList) > 0:
        s += 'Fixed Parameters:\n'
        for v in fixedVarList:
            s += '    ' + v + '\n'
    if not inputOnly:
        s += 'User-supplied parameter mapping relations:\n'
    symout = ''
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict,symGenList

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
                    symout += s1 + '\n'
                else:
                    s += s1 + '\n'
                continue
            s += '  %s = ' % mv
            j = 0
            for m,v in zip(multarr[i,:],varlist):
                if m == 0: continue
                if j > 0: s += ' + '
                j += 1
                s += '(%s * %s)' % (m,v)
            if mv in varyList: s += ' VARY'
            s += '\n'
    if symout:
        s += 'Symmetry-generated relations:\n' + symout
    if inputOnly: return s
    s += 'Inverse parameter mapping relations:\n'
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        for i,mv in enumerate(varlist):
            s += '  %s = ' % mv
            j = 0
            for m,v in zip(invmultarr[i,:],mapvars):
                if m == 0: continue
                if j > 0: s += ' + '
                j += 1
                s += '(%s * %s)' % (m,v)
            s += '\n'
    return s

def GetSymEquiv():
    '''Return the automatically generated (equivalence) relationships.

    :returns: a list of strings containing the details of the contraint relationships
    '''
    symout = []
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict,symGenList

    for varlist,mapvars,multarr,invmultarr,symFlag in zip(
        dependentParmList,indParmList,arrayList,invarrayList,symGenList):
        for i,mv in enumerate(mapvars):
            if not symFlag: continue
            if multarr is None:
                #s1 = str(mv) + ' = '
                s1 = ''
                s2 = ' = ' + str(mv)
                j = 0
                for v,m in zip(varlist,invmultarr):
                    if debug: print ('v,m[0]: ',v,m[0])
                    if len(s1.split('\n')[-1]) > 75: s1 += '\n        '
                    if j > 0: s1 += ' =  '
                    j += 1
                    s1 += str(v)
                    if m != 1:
                        s1 += " / " + str(m[0])
                symout.append(s1+s2)
                continue
            else:
                s = '  %s = ' % mv
                j = 0
                for m,v in zip(multarr[i,:],varlist):
                    if m == 0: continue
                    if j > 0: s += ' + '
                    j += 1
                    s += '(%s * %s)' % (m,v)
                print ('unexpected sym op='+s)
    return symout

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
                for v,m in zip(varlist,invmultarr):
                    if debug: print ('start dMdv',dMdv[varyList.index(name)])
                    if debug: print ('add derv',v,'/',m[0],'to derv',name,'add=',derivDict[v] / m[0])
                    if m == 0: continue
                    dMdv[varyList.index(name)] += derivDict[v] / m[0]
            else:
                for v,m in zip(varlist,multarr[i,:]):
                    if debug: print ('start dMdv',dMdv[varyList.index(name)])
                    if debug: print ('add derv',v,'*',m,'to derv',name,'add=',m * derivDict[v])
                    if m == 0: continue
                    dMdv[varyList.index(name)] += m * derivDict[v]

def Map2Dict(parmDict,varyList):
    '''Create (or update) the Independent Parameters from the original
    set of Parameters

    Removes dependent variables from the varyList

    This should be done once, after the constraints have been
    defined using :func:`StoreEquivalence`,
    :func:`GroupConstraints` and :func:`GenerateConstraints` and
    before any parameter refinement is done. This completes the parameter
    dictionary by defining independent parameters and it satisfies the
    constraint equations in the initial parameters

    :param dict parmDict: a dict containing parameter values keyed by the
      parameter names.
      This will contain updated values for both dependent and independent
      parameters after Dict2Map is called. It will also contain some
      unexpected entries of every constant value {'0':0.0} & {'1.0':1.0},
      which do not cause any problems. 

    :param list varyList: a list of parameters names that will be varied
    

    '''
    # process the Independent vars: remove dependent ones from varylist
    # and then compute values for the Independent ones from their dependents
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    for varlist,mapvars,multarr in zip(dependentParmList,indParmList,arrayList):
        for item in varlist:
            try:
                varyList.remove(item)
            except ValueError:
                pass
        if multarr is None: continue
        valuelist = [parmDict[var] for var in varlist]
        parmDict.update(zip(mapvars,
                            np.dot(multarr,np.array(valuelist)))
                        )
    # now remove fixed parameters from the varyList
    global fixedVarList
    for item in fixedVarList:
        try:
            varyList.remove(item)
        except ValueError:
            pass
    # Set constrained parameters to their fixed values
    parmDict.update(fixedDict)

def Dict2Map(parmDict,varyList):
    '''Applies the constraints defined using :func:`StoreEquivalence`,
    :func:`GroupConstraints` and :func:`GenerateConstraints` by changing
    values in a dict containing the parameters. This should be
    done before the parameters are used for any computations

    :param dict parmDict: a dict containing parameter values keyed by the
      parameter names.
      This will contain updated values for both dependent and independent
      parameters after Dict2Map is called. It will also contain some
      unexpected entries of every constant value {'0':0.0} & {'1.0':1.0},
      which do not cause any problems. 

    :param list varyList: a list of parameters names that will be varied
    
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    # reset fixed values (should not be needed, but very quick) 
    # - this seems to update parmDict with {'0':0.0} & {'1.0':1.0} - probably not what was intended
    # not needed, but also not a problem - BHT
    parmDict.update(fixedDict)
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        #if invmultarr is None: continue
        try: 
            valuelist = [parmDict[var] for var in mapvars]
        except KeyError:
            continue
        parmDict.update(zip(varlist,np.dot(invmultarr,np.array(valuelist))))

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

def _FillArray(sel,dict,collist,FillDiagonals=False):
    '''Construct a n by n matrix (n = len(collist)
    filling in the rows using the relationships defined
    in the dictionaries selected by sel

    If FillDiagonals is True, diagonal elements in the
    array are set to 1.0
    '''
    n = len(collist)
    if FillDiagonals:
        arr = np.eye(n)
    else:
        arr = np.zeros(2*[n,])
    # fill the top rows
    for i,cnum in enumerate(sel):
        for j,var in enumerate(collist):
            arr[i,j] = dict[cnum].get(var,0)
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
    varylist = ['2::atomx:3',
                '2::C(10,6,1)', '1::C(10,6,1)',
                '2::atomy:3', '2::atomz:3',
                '0:12:Scale', '0:11:Scale', '0:14:Scale', '0:13:Scale', '0:0:Scale']
#    e,w = CheckConstraints([,
#                     [{'2:0:Scale': 1.0, '5:0:Scale': 1.0, '10:0:Scale': 1.0, '6:0:Scale': 1.0, '9:0:Scale': 1.0, '8:0:Scale': 1.0,# '3:0:Scale': 1.0, '4:0:Scale': 1.0, '7:0:Scale': 1.0, '1:0:Scale': 1.0, '0:0:Scale': 1.0}],
#                     ['1.0'])
#    if e: print 'error=',e
#    if w: print 'error=',w
#    varyList = ['0::A0', '0::AUiso:0', '0::Afrac:1', '0::Afrac:2', '0::Afrac:3', '0::Afrac:4',
#       '0::dAx:5', '0::dAy:5', '0::dAz:5', '0::AUiso:5', ':0:Back;0', ':0:Back;1', ':0:Back;2', ':0:Back;3', 
#       ':0:Back;4', ':0:Back;5', ':0:Back;6', ':0:Back;7', ':0:Back;8', ':0:Back;9', ':0:Back;10', ':0:Back;11'
#       :0:U', ':0:V', ':0:W', ':0:X', ':0:Y', ':0:Scale', ':0:DisplaceX', ':0:DisplaceY']
#    constrDict = [
#        {'0::Afrac:4': 24.0, '0::Afrac:1': 16.0, '0::Afrac:3': 24.0, '0::Afrac:2': 16.0},
#        {'0::Afrac:1': 1.0, '0::Afrac:2': 1.0},
#        {'0::Afrac:4': 1.0, '0::Afrac:3': 1.0}]
#    fixedList = ['40.0', '1.0', '1.0']

    errmsg, warnmsg = CheckConstraints(varylist,constrDict,fixedList)
    if errmsg:
        print ("*** Error ********************")
        print (errmsg)
    if warnmsg:
        print ("*** Warning ********************")
        print (warnmsg)
    if errmsg or warnmsg:
        sys.exit()
    groups,parmlist = GroupConstraints(constrDict)
    GenerateConstraints(groups,parmlist,varylist,constrDict,fixedList)
    print (VarRemapShow(varylist))
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
    print ('varylist start',varylist)
    before = parmdict.copy()
    Map2Dict(parmdict,varylist)
    print ('parmdict before and after Map2Dict')
    print ('  key / before / after')
    for key in sorted(list(parmdict.keys())):
        print ('  '+key,'\t',before.get(key),'\t',parmdict[key])
    print ('varylist after',varylist)
    before = parmdict.copy()
    Dict2Map(parmdict,varylist)
    print ('after Dict2Map')
    print ('  key / before / after')
    for key in sorted(list(parmdict.keys())):
        print ('  '+key,'\t',before.get(key),'\t',parmdict[key])
#    dMdv = len(varylist)*[0]
#    deriv = {}
#    for i,v in enumerate(parmdict.keys()): deriv[v]=i
#    Dict2Deriv(varylist,deriv,dMdv)
