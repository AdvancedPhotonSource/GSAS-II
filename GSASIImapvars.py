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

Parameter redefinition (new vars) is done by creating one or more relationships
between a set of parameters

::

   Mx1 * Px + My1 * Py +...
   Mx2 * Px + Mz2 * Pz + ...

where Pj is a parameter name and Mjk is a constant.

Constant constraint Relations can also be supplied in the form of an equation:

::

  nx1 * Px + ny1 * Py +... = C1

where Cn is a constant. These equations define an algebraic
constrant.

Parameters can also be "fixed" (held), which prevents them from being refined.

All of the above three cases are input using routines
GroupConstraints and GenerateConstraints. The input consists of a list of
relationship dictionaries:

.. code-block:: python

    constrDict = [
         {'0:12:Scale': 2.0, '0:14:Scale': 4.0, '0:13:Scale': 3.0, '0:0:Scale': 0.5},
         {'2::C(10,6,1)': 1.0, '1::C(10,6,1)': 1.0},
         {'0::A0': 0.0}]
    fixedList = ['5.0', None, '0']

Where the dictionary defines the first part of an expression and the corresponding fixedList
item is either None (for parameter redefinition) or the constant value, for a constant
constraint equation. A dictionary that contains a single term defines a variable that
will be fixed (held). The multiplier and the fixedList value in this case are ignored.

Parameters can also be equivalenced or "slaved" to another parameter, such that one
(independent) parameter is equated to several (now dependent) parameters. In
algebraic form this is:

::

   P0 = M1 * P1 = M2 * P2 = ...

Thus parameters P0, P1 and P2,... are linearly equivalent. Routine StoreEquivalence is
used to specify these equivalences.

Parameter redefinition (new vars) describes a new, independent, parameter, which
is defined in terms of dependent parameters that are defined in the
Model, while fixed constrained relations effectively reduce the complexity
of the Model by removing a degree of freedom. It is possible for a parameter to appear
in both a parameter redefinition expression and a fixed constraint equation, but a
parameter cannot be used a parameter equivalance cannot be used elsewhere (not fixed,
constrained or redefined). Likewise a fixed parameter cannot be used elsewhere (not
equivalanced, constrained or redefined).

Relationships are grouped so that a set of dependent parameters appear
in only one group (done in routine GroupConstraints.) Note that if a
group contains relationships/equations that involve N dependent
parameters, there must exist N-C new parameters, where C is the number
of contraint equations in the group. Routine GenerateConstraints takes
the output from GroupConstraints and generates the
"missing" relationships and saves that information in the module's
global variables. Each generated parameter is named sequentially using paramPrefix.

A list of parameters that will be varied is specified as input to GenerateConstraints
(varyList). A fixed parameter will simply be removed from this list preventing that
parameter from being varied. Note that all parameters in a constraint relationship
must specified as varied (appear in varyList) or none can be varied. This is
checked in GenerateConstraints. Likewise, if all parameters in a constraint are
not referenced in a refinement, the constraint is ignored, but if some parameters
in a constraint group are not referenced in a refinement, but others are this
constitutes and error. 

* When a new variable is created, the variable is assigned the name associated
  in the constraint definition or it is assigned a default name of form
  ``::constr<n>`` (see paramPrefix). The vary setting for variables used in the
  constraint are ignored.
  Note that any generated "missing" relations are not varied. Only
  the input relations can be are varied.
  
* If all parameters in a fixed constraint equation are varied, the generated "missing"
  relations in the group are all varied. This provides the N-C degrees of freedom. 

*External Routines*
-------------------

To define a set of constrained and unconstrained relations, one
defines a list of dictionary defining constraint parameters and their
values, a list of fixed values for each constraint and a list of
parameters to be varied. In addition, one uses
:func:`StoreEquivalence` to define parameters that are equivalent. One
can then use :func:`CheckConstraints` to check that the input is
internally consistent and finally :func:`GroupConstraints` and
:func:`GenerateConstraints` to generate the internally used
tables. Routines :func:`Map2Dict` is used to initialize the parameter
dictionary and :func:`Dict2Map`, :func:`Dict2Deriv`, and
:func:`ComputeDepESD` are used to apply constraints. Routine
:func:`VarRemapShow` is used to print out the constraint information,
as stored by :func:`GenerateConstraints`.

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

  The latter will run more efficiently. Note that mixing independent and dependent variables is
  problematic. This is not allowed:

  .. code-block:: python

        StoreEquivalence('x',('y',))
        StoreEquivalence('y',('z',))
        
  Use StoreEquivalence before calling GenerateConstraints or CheckConstraints

:func:`CheckConstraints`
   To check that input in internally consistent, use CheckConstraints

:func:`Map2Dict`
   To determine values for the parameters created in this module, one
   calls Map2Dict. This will not apply contraints.

:func:`Dict2Map`
   To take values from the new independent parameters and constraints,
   one calls Dict2Map. This will apply contraints.

:func:`Dict2Deriv`
   Use Dict2Deriv to determine derivatives on independent parameters
   from those on dependent ones

:func:`ComputeDepESD`      
   Use ComputeDepESD to compute uncertainties on dependent variables

:func:`VarRemapShow`
   To show a summary of the parameter remapping, one calls VarRemapShow

*Global Variables*
------------------

dependentParmList:
   contains a list by group of lists of
   parameters used in the group. Note that parameters listed in
   dependentParmList should not be refined as they will not affect
   the model

indParmList:
     a list of groups of Independent parameters defined in
     each group. This contains both parameters used in parameter
     redefinitions as well as names of generated new parameters.

fixedVarList:
     a list of variables that have been 'fixed'
     by defining them as equal to a constant (::var: = 0). Note that
     the constant value is ignored at present. These variables are
     later removed from varyList which prevents them from being refined. 
     Unlikely to be used externally.

arrayList:
     a list by group of relationship matrices to relate
     parameters in dependentParmList to those in indParmList. Unlikely
     to be used externally.

invarrayList:
     a list by group of relationship matrices to relate
     parameters in indParmList to those in dependentParmList. Unlikely
     to be used externally.

fixedDict:
     a dictionary containing the fixed values corresponding
     to parameter equations.  The dict key is an ascii string, but the
     dict value is a float.  Unlikely to be used externally.

*Routines*
----------

Note that parameter names in GSAS-II are strings of form ``<ph>:<hst>:<nam>``

"""

import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
# data used for constraints; 
debug = False # turns on printing as constraint input is processed
# note that constraints are stored listed by contraint groups, where each constraint
# group contains those parameters that must be handled together
dependentParmList = [] # contains a list of parameters in each group
# note that parameters listed in dependentParmList should not be refined 
arrayList = [] # a list of of relationship matrices 
invarrayList = [] # a list of inverse relationship matrices 
indParmList = [] # a list of names for the new parameters
fixedDict = {} # a dictionary containing the fixed values corresponding to defined parameter equations
               # key is original ascii string, value is float
fixedVarList = [] # List of variables that should not be refined

# prefix for parameter names
paramPrefix = "::constr"
consNum = 0 # number of the next constraint to be created

def InitVars():
    '''Initializes all constraint information'''
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict,consNum
    dependentParmList = [] # contains a list of parameters in each group
    arrayList = [] # a list of of relationship matrices 
    invarrayList = [] # a list of inverse relationship matrices 
    indParmList = [] # a list of names for the new parameters
    fixedDict = {} # a dictionary containing the fixed values corresponding to defined parameter equations
    consNum = 0 # number of the next constraint to be created
    fixedVarList = []

def VarKeys(constr):
    """Finds the keys in a constraint that represent variables
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
    errmsg = ''
    warnmsg = ''
    fixVlist = []
    # process fixed variables (holds)
    for cdict in constrDict:
        # N.B. No "_" names in holds
        if len(cdict) == 1:
            fixVlist.append(cdict.keys()[0])
    
    # process equivalences: make a list of dependent and independent vars
    #    and check for repeated uses (repetition of a parameter as an
    #    independent var is OK)
    indepVarList = []
    depVarList = []
    multdepVarList = []
    for varlist,mapvars,multarr,invmultarr in zip(
        dependentParmList,indParmList,arrayList,invarrayList):
        if multarr is None: # an equivalence
            zeromult = False
            for mv in mapvars:
                varied = 0
                notvaried = ''
                if mv in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += mv
                if mv not in indepVarList: indepVarList.append(mv)
                for v,m in zip(varlist,invmultarr):
                    if v in indepVarList:
                        errmsg += '\nVariable '+v+' is used to set values in a constraint before its value is set in another constraint\n'
                    if m == 0: zeromult = True
                    if v in varyList:
                        varied += 1
                    else:
                        if notvaried: notvaried += ', '
                        notvaried += v
                    if v in depVarList:
                        multdepVarList.append(v)
                    else:
                        depVarList.append(v)
            if varied > 0 and varied != len(varlist)+1:
                warnmsg += "\nNot all variables refined in equivalence:\n\t"
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
            if s != "": s+= ", "
            s += str(var)            
        errmsg += '\t'+ s + '\n'
    equivVarList = list(set(indepVarList).union(set(depVarList)))
    if debug: print 'equivVarList',equivVarList
    inboth = set(fixVlist).intersection(set(equivVarList))
    if len(inboth) > 0:
        errmsg += "\nThe following parameter(s) are used in both Equivalence and Fixed constraints:\n"
        s = ''
        for var in sorted(inboth):
            if s != "": s+= ", "
            s += str(var)
        errmsg += '\t'+ s + '\n'

    groups,parmlist = GroupConstraints(constrDict)
    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create a warning message.
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        for rel in group:
            varied = 0
            notvaried = ''
            for var in constrDict[rel]:
                if var.startswith('_'): continue
                if not re.match('[0-9]*:[0-9\*]*:',var):
                    warnmsg += "\nVariable "+str(var)+" does not begin with a ':'"
                if var in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += var
                if var in fixVlist:
                    errmsg += '\nParameter '+var+" is Fixed and used in a constraint:\n\t"
                    errmsg += _FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
            if varied > 0 and varied != len(VarKeys(constrDict[rel])):
                warnmsg += "\nNot all variables refined in constraint:\n\t"
                warnmsg += _FormatConstraint(constrDict[rel],fixedList[rel])
                warnmsg += '\nNot refined: ' + notvaried + '\n'
    if errmsg or warnmsg:
        return errmsg,warnmsg

    # now look for process each group and create the relations that are needed to form
    # non-singular square matrix
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue # a constraint group with a single variable can be ignored
        if len(varlist) < len(group): # too many relationships -- no can do
            errmsg += "\nOver-constrained input. "
            errmsg += "There are more constraints " + str(len(group))
            errmsg += "\n\tthan variables " + str(len(varlist)) + "\n"
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
        # scan through all generated and input variables
        # Check again for inconsistent variable use
        # for new variables -- where varied and unvaried parameters get grouped
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
                varname = paramPrefix + str(consNum) # assign a name to a variable
                mapvar.append(varname)
                consNum += 1
            else:
                mapvar.append(fixedval)
            if varied > 0 and notvaried != '':
                warnmsg += "\nNot all variables refined in generated constraint"
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
    return errmsg,warnmsg

def GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict=None,SeqHist=None):
    '''Takes a list of relationship entries comprising a group of
    constraints and builds the relationship lists and their inverse
    and stores them in global variables Also checks for internal
    conflicts or inconsistencies in parameter/variable definitions.

    :param list groups: a list of grouped contraints where each constraint
      grouped containts a list of indices for constraint constrDict entries,
      created in :func:`GroupConstraints` (returned as 1st value)

    :param list parmlist: a list containing lists of parameter names
      contained in each group, created in :func:`GroupConstraints`
      (returned as 2nd value)

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
      refinement. None (default) otherwise. Wildcard variable names are
      set to the current histogram, when found if not None.
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    msg = ''

    # process fixed (held) variables
    for cdict in constrDict:
        if len(cdict) == 1:
            fixedVarList.append(cdict.keys()[0])
    
    # process equivalences: make a list of dependent and independent vars
    #    and check for repeated uses (repetition of a parameter as an
    #    independent var is OK [A=B; A=C], but chaining: [A=B; B=C] is not good)
    dropVarList = []
    translateTable = {} # lookup table for wildcard referenced variables
    for varlist,mapvars,multarr,invmultarr in zip(       # process equivalences
        dependentParmList,indParmList,arrayList,invarrayList):
        if multarr is None: # true only if an equivalence
            zeromult = False
            for mv in mapvars:
                #s = ''
                varied = 0
                notvaried = ''
                if mv in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += mv
                if parmDict is not None and mv not in parmDict:
                    print "Dropping equivalence for variable "+str(mv)+". Not defined in this refinement"
                    if mv not in dropVarList: dropVarList.append(mv)
                    #msg += "\nCannot equivalence to variable "+str(mv)+". Not defined in this refinement"
                    #continue
            for v,m in zip(varlist,invmultarr):
                if parmDict is not None and v not in parmDict:
                    print "Dropping equivalence for dep. variable "+str(v)+". Not defined in this refinement"
                    if v not in dropVarList: dropVarList.append(v)
                    continue
                if m == 0: zeromult = True
                if v in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += v
            if varied > 0 and varied != len(varlist)+1:
                msg += "\nNot all variables refined in equivalence:\n\t"
                s = ""
                for v in varlist:
                    if s != "": s+= " & "
                    s += str(v)            
                msg += str(mv) + " => " + s
                msg += '\nNot refined: ' + notvaried + '\n'
            if zeromult:
                msg += "\nZero multiplier is invalid in equivalence:\n\t"
                s = ""
                for v in varlist:
                    if s != "": s+= " & "
                    s += str(v)            
                msg += str(mv) + " => " + s + '\n'

    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create an error message. 
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
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
                msg += "\nSome (but not all) variables in constraint are not defined:\n\t"
                msg += _FormatConstraint(constrDict[rel],fixedList[rel])
                msg += '\nNot used: ' + notused + '\n'
            if varied > 0 and varied != len(VarKeys(constrDict[rel])):
                msg += "\nNot all variables refined in constraint:\n\t"
                msg += _FormatConstraint(constrDict[rel],fixedList[rel])
                msg += '\nNot refined: ' + notvaried + '\n'
    # if there were errors found, go no farther
    if msg:
        print ' *** ERROR in constraint definitions! ***'
        print msg
        raise Exception
                
    # now process each group and create the relations that are needed to form
    # a non-singular square matrix
    # If all are varied and this is a constraint equation, then set VaryFree flag
    # so that the newly created relationships will be varied
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        # for constraints, if all included variables are refined,
        # set the VaryFree flag, and remaining degrees of freedom will be
        # varied (since consistency was checked, if any one variable is
        # refined, then assume that all are)
        varsList = [] # make a list of all the referenced variables as well
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
        # Since we checked before, if any variables are unused, then all must be. 
        # If so, this set of relationships can be ignored
        if unused:
            if debug: print('Constraint ignored (all variables undefined)')
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
            if fixedval is None: # this is a new variable, not a constraint
                if not varname:
                    varname = paramPrefix + str(consNum) # no assigned name, create one
                    consNum += 1
                mapvar.append(varname)
                # vary the new relationship if it is a degree of freedom in
                # a set of contraint equations or if a New Var is flagged to be varied.
                if VaryFree or varyflag: 
                    unused = False
                    varyList.append(varname)
                    # fix (prevent varying) of all the variables inside the constraint group
                    # (dependent vars)
                    for var in varsList:
                        if var in varyList: varyList.remove(var)
            else:
                unused = False
                mapvar.append(fixedval)
        if unused: # skip over constraints that don't matter (w/o fixed value or any refined variables)
            if debug: print('Constraint ignored (all variables unrefined)')
            if debug: print ('   '+_FormatConstraint(constrDict[rel],fixedList[rel]))
            continue 
        dependentParmList.append([translateTable.get(var,var) for var in varlist])
        arrayList.append(constrArr)
        invarrayList.append(np.linalg.inv(constrArr))
        indParmList.append(mapvar)
    if msg:
        print ' *** ERROR in constraint definitions! ***'
        print msg
        print VarRemapShow(varyList)
        raise Exception
    # setup dictionary containing the fixed values
    global fixedDict 
    # key is original ascii string, value is float
    for fixedval in fixedList:
        if fixedval:
            fixedDict[fixedval] = float(fixedval)

    # make list of dependent and independent variables (after dropping unused)
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
        print 50*'-'
        print VarRemapShow(varyList)
        print 'Varied: ',varyList
        print 'Not Varied: ',fixedVarList

def StoreEquivalence(independentVar,dependentList):
    '''Takes a list of dependent parameter(s) and stores their
    relationship to a single independent parameter (independentVar)

    :param str independentVar: name of master parameter that will be used to determine the value
      to set the dependent variables

    :param list dependentList: a list of parameters that will set from
         independentVar. Each item in the list can be a string with the parameter
         name or a tuple containing a name and multiplier:
         ``['parm1',('parm2',.5),]``

    '''
    
    global dependentParmList,arrayList,invarrayList,indParmList
    mapList = []
    multlist = []
    for var in dependentList:
        if isinstance(var, basestring):
            mult = 1.0
        elif len(var) == 2:
            var,mult = var
        else:
            raise Exception("Cannot parse "+repr(var) + " as var or (var,multiplier)")
        mapList.append(var)
        multlist.append(tuple((mult,)))
    # added relationships to stored values
    arrayList.append(None)
    invarrayList.append(np.array(multlist))
    indParmList.append(tuple((independentVar,)))
    dependentParmList.append(mapList)
    return

def GetDependentVars():
    '''Return a list of dependent variables: e.g. variables that are
    constrained in terms of other variables

    :returns: a list of variable names

    '''
    return dependentVars

def GetIndependentVars():
    '''Return a list of independent variables: e.g. variables that are
    created by constraints of other variables

    :returns: a list of variable names

    '''
    return independentVars

def PrintIndependentVars(parmDict,varyList,sigDict,PrintAll=False,pFile=None):
    '''Print the values and uncertainties on the independent variables'''
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
    print >>pFile,130*'-'
    print >>pFile,"Variables generated by constraints"
    printlist.append(3*[None])
    for name,val,esd in printlist:
        if len(s1) > 120 or name is None:
            print >>pFile,''
            print >>pFile,s1
            print >>pFile,s2
            print >>pFile,s3
            s1 = ''
            if name is None: break
        if s1 == "":
            s1 = ' name  :'
            s2 = ' value :'
            s3 = ' sig   :'
        s1 += '%15s' % (name)
        s2 += '%15.4f' % (val)
        if esd is None:
            s3 += '%15s' % ('n/a')
        else:    
            s3 += '%15.4f' % (esd)

def ComputeDepESD(covMatrix,varyList,parmDict):
    '''Compute uncertainties for dependent parameters from independent ones
    returns a dictionary containing the esd values for dependent parameters
    '''
    sigmaDict = {}
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        #if invmultarr is None: continue # probably not needed
        try: 
            valuelist = [parmDict[var] for var in mapvars]
        except KeyError:
            continue
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
        s += 'Fixed Variables:\n'
        for v in fixedVarList:
            s += '    ' + v + '\n'
    s += 'Variable mapping relations:\n'
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    for varlist,mapvars,multarr,invmultarr in zip(
        dependentParmList,indParmList,arrayList,invarrayList):
        for i,mv in enumerate(mapvars):
            if multarr is None:
                s += '  ' + str(mv) + ' is equivalent to parameter(s): '
                j = 0
                for v,m in zip(varlist,invmultarr):
                    if debug: print 'v,m[0]: ',v,m[0]
                    if j > 0: s += '  & '
                    j += 1
                    s += str(v)
                    if m != 1:
                        s += " / " + str(m[0])                        
                s += '\n'
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
    if inputOnly: return s
    s += 'Inverse variable mapping relations:\n'
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
            # grouped variables: need to add in the derv. w/r
            # dependent variables to the independent ones
            if name not in varyList: continue # skip if independent var not varied
            if multarr is None:
                for v,m in zip(varlist,invmultarr):
                    if debug: print 'start dMdv',dMdv[varyList.index(name)]
                    if debug: print 'add derv',v,'/',m[0],'to derv',name,'add=',derivDict[v] / m[0]
                    if m == 0: continue
                    dMdv[varyList.index(name)] += derivDict[v] / m[0]
            else:
                for v,m in zip(varlist,multarr[i,:]):
                    if debug: print 'start dMdv',dMdv[varyList.index(name)]
                    if debug: print 'add derv',v,'*',m,'to derv',name,'add=',m * derivDict[v]
                    if m == 0: continue
                    dMdv[varyList.index(name)] += m * derivDict[v]

def Map2Dict(parmDict,varyList):
    '''Create (or update) the Independent Parameters from the original
    set of Parameters

    Removes dependent variables from the varyList

    This should be done once, after the constraints have been
    defined using :func:`StoreEquivalence`,
    :func:`GroupConstraints` and :func:`GenerateConstraints` and
    before any variable refinement is done. This completes the parameter
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
    # now remove fixed variables from the varyList
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
        parmDict.update(zip(varlist,
                            np.dot(invmultarr,np.array(valuelist)))
                        )

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
            raise Exception,"Singular input to GramSchmidtOrtho"
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
        raise Exception,'Singular input'

def _RowEchelon(m,arr,collist):
    '''Convert the first m rows in Matrix arr to row-echelon form
    exchanging columns in the matrix and collist as needed.

    throws an exception if the matrix is singular because
    the first m rows are not linearly independent
    '''
    n = len(collist)
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
#    varyList = ['0::A0', '0::AUiso:0', '0::Afrac:1', '0::Afrac:2', '0::Afrac:3', '0::Afrac:4', '0::dAx:5', '0::dAy:5', '0::dAz:5', '0::AUiso:5', ':0:Back;0', ':0:Back;1', ':0:Back;2', ':0:Back;3', ':0:Back;4', ':0:Back;5', ':0:Back;6', ':0:Back;7', ':0:Back;8', ':0:Back;9', ':0:Back;10', ':0:Back;11', ':0:U', ':0:V', ':0:W', ':0:X', ':0:Y', ':0:Scale', ':0:DisplaceX', ':0:DisplaceY']
#    constrDict = [
#        {'0::Afrac:4': 24.0, '0::Afrac:1': 16.0, '0::Afrac:3': 24.0, '0::Afrac:2': 16.0},
#        {'0::Afrac:1': 1.0, '0::Afrac:2': 1.0},
#        {'0::Afrac:4': 1.0, '0::Afrac:3': 1.0}]
#    fixedList = ['40.0', '1.0', '1.0']

    errmsg, warnmsg = CheckConstraints(varylist,constrDict,fixedList)
    if errmsg:
        print "*** Error ********************"
        print errmsg
    if warnmsg:
        print "*** Warning ********************"
        print warnmsg
    if errmsg or warnmsg:
        sys.exit()
    groups,parmlist = GroupConstraints(constrDict)
    GenerateConstraints(groups,parmlist,varylist,constrDict,fixedList)
    print VarRemapShow(varylist)
    parmdict.update( {
        '0:12:Scale': 1.0, '0:11:Scale': 1.0, '0:14:Scale': 1.0, '0:13:Scale': 1.0, '0:0:Scale': 2.0,
        '0:0:eA': 0.0,
        '2::C(10,6,1)': 0.2, '1::C(10,6,1)': 0.3,
        '1::C(10,0,1)': 0.2, '2::C(10,0,1)': 0.3,
        '1::AUiso:0': 0.02, '0::AUiso:0': 0.03,
        '0::A0': 0.0,
        '2::atomx:3':0.23,'2::atomy:3':-.23, '2::atomz:3':-0.11,
        })
    print 'parmdict start',parmdict
    print 'varylist start',varylist
    before = parmdict.copy()
    Map2Dict(parmdict,varylist)
    print 'parmdict before and after Map2Dict'
    print '  key / before / after'
    for key in sorted(parmdict.keys()):
        print '  '+key,'\t',before.get(key),'\t',parmdict[key]
    print 'varylist after',varylist
    before = parmdict.copy()
    Dict2Map(parmdict,varylist)
    print 'after Dict2Map'
    print '  key / before / after'
    for key in sorted(parmdict.keys()):
        print '  '+key,'\t',before.get(key),'\t',parmdict[key]
#    dMdv = len(varylist)*[0]
#    deriv = {}
#    for i,v in enumerate(parmdict.keys()): deriv[v]=i
#    Dict2Deriv(varylist,deriv,dMdv)
