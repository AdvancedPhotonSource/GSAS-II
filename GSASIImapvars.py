########### SVN repository information ###################
# $Date: 2011-01-07 13:27:24 -0600 (Fri, 07 Jan 2011) $
# $Author: vondreele & toby $
# $Revision: 230 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIImapvars.py $
# $Id: GSASIImapvars.py 230 2011-01-07 19:27:24Z vondreele $
########### SVN repository information ###################
"""Module that implements algebraic contraints, parameter redefinition
and parameter simplification contraints.

Parameter redefinition (new vars) is done by creating one or more relationships
between a set of parameters
   Mx1 * Px + My1 * Py +...
   Mx2 * Px + Mz2 * Pz + ...
where Pj is a parameter name and Mjk is a constant.

Constant constraint Relations can also be supplied in the form of an equation:
  nx1 * Px + ny1 * Py +... = C1
where Cn is a constant. These equations define an algebraic
constrant.

Parameters can also be "fixed" (held), which prevents them from being refined.

All of the above three cases is supplied as input using routines
GroupConstraints and GenerateConstraints. The input consists of a list of
relationship dictionaries:
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
parameter from being varied. Note that all parameters in a relationship must specified as
varied (appear in varyList) or none can be varied. This is checked in GenerateConstraints
(as well as for generated relationships in SetVaryFlags).
* If all parameters in a parameter redefinition (new var) relationship are varied, the
  parameter assigned to this expression (::constr:n, see paramPrefix) newly generated
  parameter is varied. Note that any generated "missing" relations are not varied. Only
  the input relations are varied.
* If all parameters in a fixed constraint equation are varied, the generated "missing"
  relations in the group are all varied. This provides the N-C degrees of freedom. 

External Routines:
   To define a set of constrained and unconstrained relations, one
     calls InitVars, GroupConstraints and GenerateConstraints. It
     is best to supply all relations/equations in one call to these
     routines so that grouping is done correctly.

   To implement parameter redefinition, one calls
     StoreEquivalence. This should be called for every set of
     equivalence relationships. There is no harm in using
     StoreEquivalence with the same independent variable:
       StoreEquivalence('x',('y',))
       StoreEquivalence('x',('z',))
     (though StoreEquivalence('x',('y','z')) will run more
     efficiently) but mixing independent and dependent variables is
     problematic. This is not allowed:
        StoreEquivalence('x',('y',))
        StoreEquivalence('y',('z',))
   Use StoreEquivalence before calling GenerateConstraints or
      CheckConstraints

   To check that input in internally consistent, use CheckConstraints

   To show a summary of the parameter remapping, one calls
      VarRemapShow

   To determine values for the parameters created in this module, one
      calls Map2Dict. This will not apply contraints.

   To take values from the new independent parameters and constraints,
      one calls Dict2Map. This will apply contraints.

   Use Dict2Deriv to determine derivatives on independent parameters
      from those on dependent ones
      
   Use ComputeDepESD to compute uncertainties on dependent variables

Global Variables:
   dependentParmList: contains a list by group of lists of
     parameters used in the group. Note that parameters listed in
     dependentParmList should not be refined as they will not affect
     the model
   indParmList: a list of groups of Independent parameters defined in
     each group. This contains both parameters used in parameter
     redefinitions as well as names of generated new parameters.

   fixedVarList: a list of variables that have been 'fixed'
     by defining them as equal to a constant (::var: = 0). Note that
     the constant value is ignored at present. These variables are
     later removed from varyList which prevents them from being refined. 
     Unlikely to be used externally.
   arrayList: a list by group of relationship matrices to relate
     parameters in dependentParmList to those in indParmList. Unlikely
     to be used externally.
   invarrayList: a list by group of relationship matrices to relate
     parameters in indParmList to those in dependentParmList. Unlikely
     to be used externally.
   fixedDict: a dictionary containing the fixed values corresponding
     to parameter equations.  The dict key is an ascii string, but the
     dict value is a float.  Unlikely to be used externally.
"""

import numpy as np
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
paramPrefix = "::constr:"
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

def GroupConstraints(constrDict):
    """divide the constraints into groups that share no parameters.
    Input
       constrDict: a list of dicts defining relationships/constraints
       constrDict = [{<constr1>}, {<constr2>}, ...]
       where {<constr1>} is {'param1': mult1, 'param2': mult2,...}
    Returns two lists of lists:
      a list of relationship list indices for each group pointing to
        elements in constrDict and
      a list containing the parameters used in each group
      """
    assignedlist = [] # relationships that have been used
    groups = [] # contains a list of grouplists
    ParmList = []
    for i,consi in enumerate(constrDict):
        if i in assignedlist: continue # already in a group, skip
        # starting a new group
        grouplist = [i,]
        assignedlist.append(i)
        groupset = set(consi.keys())
        changes = True # always loop at least once
        while(changes): # loop until we can't find anything to add to the current group
            changes = False # but don't loop again unless we find something
            for j,consj in enumerate(constrDict):
                if j in assignedlist: continue # already in a group, skip
                if len(set(consj.keys()) & groupset) > 0: # true if this needs to be added
                    changes = True
                    grouplist.append(j)
                    assignedlist.append(j)
                    groupset = groupset | set(consj.keys())
        group = sorted(grouplist)
        varlist = sorted(list(groupset))
        groups.append(group)
        ParmList.append(varlist)
    return groups,ParmList

def CheckConstraints(varyList,constrDict,fixedList):
    '''Takes a list of relationship entries comprising a group of
    constraints and checks for inconsistencies such as conflicts in
    parameter/variable definitions and or inconsistently varied parameters.
    Input: see GenerateConstraints
    Output: returns two strings
      the first lists conflicts internal to the specified constraints
      the second lists conflicts where the varyList specifies some
        parameters in a constraint, but not all
      If there are no errors, both strings will be empty
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    errmsg = ''
    warnmsg = ''
    fixVlist = []
    # process fixed (held) variables
    for cdict in constrDict:
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

    #print 'indepVarList',indepVarList
    #print 'depVarList',depVarList
    # check for errors:
    inboth = set(indepVarList).intersection(set(depVarList))
    if len(inboth) > 0:
        errmsg += "\nThe following parameters(s) are used as both dependent and independent variables in Equivalence relations:\n"
        s = ''
        for var in inboth:
            if s != "": s+= ", "
            s += str(v)
        errmsg += '\t'+ s + '\n'
    if len(multdepVarList) > 0:
        errmsg += "\nThe following parameters(s) are used in multiple Equivalence relations as dependent variables:\n"
        s = ''
        for var in multdepVarList:
            if s != "": s+= ", "
            s += str(v)            
        errmsg += '\t'+ s + '\n'
    equivVarList = list(set(indepVarList).union(set(depVarList)))
    #print 'equivVarList',equivVarList
    inboth = set(fixVlist).intersection(set(equivVarList))
    if len(inboth) > 0:
        errmsg += "\nThe following parameter(s) are used in both Equivalence and Fixed constraints:\n"
        s = ''
        for var in inboth:
            if s != "": s+= ", "
            s += str(v)
        errmsg += '\t'+ s + '\n'

    groups,parmlist = GroupConstraints(constrDict)
    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create a warning message.
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        VaryFree = False
        for rel in group:
            varied = 0
            notvaried = ''
            for var in constrDict[rel]:
                if var in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += var
                if var in fixVlist:
                    errmsg += '\nParameter '+var+" is Fixed and used in a constraint:\n\t"
                    errmsg += FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
                if var in equivVarList:
                    errmsg += '\nParameter '+var+" is Equivalenced and used in a constraint:\n\t"
                    errmsg += FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
            if varied > 0 and varied != len(constrDict[rel]):
                warnmsg += "\nNot all variables refined in constraint:\n\t"
                warnmsg += FormatConstraint(constrDict[rel],fixedList[rel])
                warnmsg += '\nNot refined: ' + notvaried + '\n'
    if errmsg or warnmsg: return errmsg,warnmsg

    # now look for process each group and create the relations that are needed to form
    # non-singular square matrix
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        try: 
            arr = MakeArray(constrDict, group, varlist)
        except:
            errmsg += "\nOver-constrained input. "
            errmsg += "There are more constraints" + str(len(group))
            errmsg += "\n\tthan variables" + str(len(varlist)) + "\n"
            for rel in group:
                errmsg += FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
        try:
            constrArr = FillInMissingRelations(arr,len(group))
        except Exception,msg:
            if msg.message.startswith('VectorProduct'):
                errmsg += "\nSingular input. "
                errmsg += "This set of constraints is not linearly independent:\n\n"
            else:
                errmsg += "\nInconsistent input. "
                errmsg += "Cound not generate a set of non-singular constraints"
                errmsg += "\n\tfor this set of input:\n\n"
            for rel in group:
                errmsg += FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
            #import traceback
            #print traceback.format_exc()
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
                for var in constrDict[rel]:
                    if var in varyList:
                        varied += 1
                    else:
                        if notvaried: notvaried += ', '
                        notvaried += var
            else:
                fixedval = None
            if fixedval is None:
                varname = paramPrefix + str(consNum)
                mapvar.append(varname)
                consNum += 1
                if VaryFree or varied > 0:
                    varyList.append(varname)
            else:
                mapvar.append(fixedval)
            if varied > 0 and notvaried != '':
                warnmsg += "\nNot all variables refined in generated constraint"
                warnmsg += '\nPlease report this unexpected error\n'
                for rel in group:
                    warnmsg += FormatConstraint(constrDict[rel],fixedList[rel])
                    warnmsg += "\n"
                warnmsg += '\n\tNot refined: ' + notvaried + '\n'
        try:
            np.linalg.inv(constrArr)
        except:
            errmsg += "\nSingular input. "
            errmsg += "The following constraints are not "
            errmsg += "linearly independent\n\tor do not "
            errmsg += "allow for generation of a non-singular set\n"
            errmsg += 'Please report this unexpected error\n'
            for rel in group:
                errmsg += FormatConstraint(constrDict[rel],fixedList[rel])
                errmsg += "\n"
    return errmsg,warnmsg

def GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList):
    '''Takes a list of relationship entries comprising a group of
    constraints and builds the relationship lists and their inverse
    and stores them in global variables Also checks for internal
    conflicts or inconsistencies in parameter/variable definitions.
    Input:
       groups,parmlist: see GroupConstraints
       varyList: a list of parameters that will be varied
       constrDict: a list of dicts defining relationships/constraints
         constrDict = [{<constr1>}, {<constr2>}, ...]
         where {<constr1>} is {'param1': mult1, 'param2': mult2,...}
       fixedList: a list of values for each constraint/variable in
          constrDict, value is either a float (contraint) or None (for
          a new variable).
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    msg = ''

    # process fixed (held) variables
    for cdict in constrDict:
        if len(cdict) == 1:
            fixedVarList.append(cdict.keys()[0])
    #print 'fixedVarList',fixedVarList
    
    # process equivalences: make a list of dependent and independent vars
    #    and check for repeated uses (repetition of a parameter as an
    #    independent var is OK)
    indepVarList = []
    depVarList = []
    multdepVarList = []
    for varlist,mapvars,multarr,invmultarr in zip(
        dependentParmList,indParmList,arrayList,invarrayList):
        if multarr is None:
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

    #print 'indepVarList',indepVarList
    #print 'depVarList',depVarList
    # check for errors:
    inboth = set(indepVarList).intersection(set(depVarList))
    if len(inboth) > 0:
        msg += "\nThe following parameters(s) are used as both dependent and independent variables in Equivalence relations:\n"
        s = ''
        for var in inboth:
            if s != "": s+= ", "
            s += str(v)
        msg += '\t'+ s + '\n'
    if len(multdepVarList) > 0:
        msg += "\nThe following parameters(s) are used in multiple Equivalence relations as dependent variables:\n"
        s = ''
        for var in multdepVarList:
            if s != "": s+= ", "
            s += str(v)            
        msg += '\t'+ s + '\n'
    equivVarList = list(set(indepVarList).union(set(depVarList)))
    #print 'equivVarList',equivVarList
    inboth = set(fixedVarList).intersection(set(equivVarList))
    if len(inboth) > 0:
        msg += "\nError! The following variables are used in both Equivalence and Fixed constraints:\n"
        s = ''
        for var in inboth:
            if s != "": s+= ", "
            s += str(v)
        msg += '\t'+ s + '\n'

    # scan through parameters in each relationship. Are all varied? If only some are
    # varied, create an error message. 
    # If all are varied and this is a constraint equation, then set VaryFree flag
    # so that newly created relationships will be varied
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        VaryFree = False
        for rel in group:
            varied = 0
            notvaried = ''
            for var in constrDict[rel]:
                if var in varyList:
                    varied += 1
                else:
                    if notvaried: notvaried += ', '
                    notvaried += var
                if var in fixedVarList:
                    msg += '\nError: parameter '+var+" is Fixed and used in a constraint:\n\t"
                    msg += FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
                if var in equivVarList:
                    msg += '\nError: parameter '+var+" is Equivalenced and used in a constraint:\n\t"
                    msg += FormatConstraint(constrDict[rel],fixedList[rel])+"\n"
            if varied > 0 and varied != len(constrDict[rel]):
                msg += "\nNot all variables refined in constraint:\n\t"
                msg += FormatConstraint(constrDict[rel],fixedList[rel])
                msg += '\nNot refined: ' + notvaried + '\n'
            if fixedList[rel] is not None and varied > 0:
                VaryFree = True
    # if there were errors found, go no farther
    if msg:
        print ' *** ERROR in constraint definitions! ***'
        print msg
        raise Exception
                
    # now process each group and create the relations that are needed to form
    # non-singular square matrix
    for group,varlist in zip(groups,parmlist):
        if len(varlist) == 1: continue
        arr = MakeArray(constrDict, group, varlist)
        constrArr = FillInMissingRelations(arr,len(group))
        mapvar = []
        group = group[:]
        # scan through all generated and input variables, add to the varied list
        # all the new parameters where VaryFree has been set or where all the
        # dependent parameters are varied. Check again for inconsistent variable use
        # for new variables -- where varied and unvaried parameters get grouped
        # together. I don't think this can happen when not flagged before, but
        # it does not hurt to check again. 
        for i in range(len(varlist)):
            varied = 0
            notvaried = ''
            if len(group) > 0:
                rel = group.pop(0)
                fixedval = fixedList[rel]
                for var in constrDict[rel]:
                    if var in varyList:
                        varied += 1
                    else:
                        if notvaried: notvaried += ', '
                        notvaried += var
            else:
                fixedval = None
            if fixedval is None:
                varname = paramPrefix + str(consNum)
                mapvar.append(varname)
                consNum += 1
                if VaryFree or varied > 0:
                    varyList.append(varname)
            else:
                mapvar.append(fixedval)
            if varied > 0 and notvaried != '':
                msg += "\nNot all variables refined in generated constraint\n\t"
                msg += '\nNot refined: ' + notvaried + '\n'
        dependentParmList.append(varlist)
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

    if debug: # on debug, show what is parsed & generated, semi-readable
        print 50*'-'
        for group,varlist,multarr,inv,mapvar in zip(groups,parmlist,arrayList,invarrayList,indParmList):
            print '\n*** relation(s) in group:',group,'\tvars in group:',varlist
            print 'new parameters:', mapvar
            print 'Input relationship matrix'
            print multarr[:len(group)]
            added = len(group) - len(varlist)
            if added < 0:
                print 'added relationship rows'
                print multarr[added:]
            print 'Inverse relationship matrix'
            print inv

def StoreEquivalence(independentVar,dependentList):
    '''Takes a list of dependent parameter(s) and stores their
    relationship to a single independent parameter (independentVar)
    input:
       independentVar -- name of parameter that will set others (may
         be varied)
       dependentList -- list of parameter that will set from
         independentVar (may not be varied). Each item can be a parameter
         name or a tuple containing a name and multiplier:
         ('parm1',('parm2',.5),)
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
            raise Exception, "Cannot parse "+repr(var) + " as var or (var,multiplier)"
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
    constrained in terms of other variables'''
    dependentVars = []
    global dependentParmList
    for lst in dependentParmList:
        for itm in lst: dependentVars.append(itm)
    return dependentVars

def GetIndependentVars():
    '''Return a list of independent variables: e.g. variables that are
    created by constraints of other variables'''
    independentVars = []
    global indParmList,fixedDict
    for lst in indParmList:
        for name in lst:
            if name in fixedDict: continue
            independentVars.append(name)
    return independentVars

def PrintIndependentVars(parmDict,varyList,sigDict,PrintAll=False):
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
    print 130*'-'
    print "Variables generated by constraints"
    printlist.append(3*[None])
    for name,val,esd in printlist:
        if len(s1) > 110 or name is None:
            print
            print s1
            print s2
            print s3
            s1 = ''
            if name is None: break
        if s1 == "":
            s1 = ' name  :'
            s2 = ' value :'
            s3 = ' sig   :'
        s1 += '%12s' % (name)
        s2 += '%12.6f' % (val)
        if esd is None:
            s3 += '%12s' % ('n/a')
        else:    
            s3 += '%12.6f' % (esd)

def ComputeDepESD(covMatrix,varyList,parmDict):
    '''Compute uncertainties for dependent parameters from independent ones
    returns a dictionary containing the esd values for dependent parameters
    '''
    sigmaDict = {}
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        #if invmultarr is None: continue # probably not needed
        valuelist = [parmDict[var] for var in mapvars]
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

def FormatConstraint(RelDict,RelVal):
    '''Formats a Constraint or Function for use in a convenient way'''
    linelen = 45
    s = [""]
    for var,val in RelDict.items():
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

def VarRemapShow(varyList):
    '''List out the saved relationships.
    Returns a string containing the details.
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
                    #print v,m[0]
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
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,invarrayList
    for varlist,mapvars,multarr,invmultarr in zip(dependentParmList,indParmList,arrayList,invarrayList):
        for i,name in enumerate(mapvars):
            # grouped variables: need to add in the derv. w/r
            # dependent variables to the independent ones
            if name not in varyList: continue # skip if independent var not varied
            if multarr is None:
                for v,m in zip(varlist,invmultarr):
                    #print 'add derv',v,'/',m[0],'to derv',name
                    if m == 0: continue
                    dMdv[varyList.index(name)] += derivDict[v] / m[0]
            else:
                for v,m in zip(varlist,multarr[i,:]):
                    #print 'add derv',v,'*',m,'to derv',name
                    if m == 0: continue
                    dMdv[varyList.index(name)] += m * derivDict[v]

def Map2Dict(parmDict,varyList):
    '''Create (or update) the Independent Parameters from the original
    set of Parameters

    Removes dependent variables from the varyList

    This should be done once, before any variable refinement is done
    to complete the parameter dictionary with the Independent Parameters
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
    '''Convert the remapped values back to the original parameters
    
    This should be done to apply constraints to parameter values (use
    Map2Dict once first). It should be done as part of the Model function
    evaluation, before any computation is done
    '''
    #I think this needs fixing to update parmDict with new values 
    #   from the constraints based on what is in varyList - RVD. Don't think so -- BHT
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    # reset fixed values (should not be needed, but very quick) 
    # - this seems to update parmDict with {'0':0.0} & {'1.0':1.0} - probably not what was intended
    # not needed, but also not a problem - BHT
    parmDict.update(fixedDict)
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        #if invmultarr is None: continue
        valuelist = [parmDict[var] for var in mapvars]
        parmDict.update(zip(varlist,
                            np.dot(invmultarr,np.array(valuelist)))
                        )

#======================================================================
# internal routines follow (these routines are unlikely to be called
# from outside the module)
def FillInMissingRelations(arrayin,nvars):
    '''Fill in unspecified rows in array with non-colinear unit vectors
    arrayin is a square array, where only the first nvars rows are defined.
    
    Returns a new array where all rows are defined

    Throws an exception if there is no non-signular result
    (likely because two or more input rows are co-linear)
    '''
    a = arrayin.copy()
    n = nvars
    # fill in the matrix with basis vectors not colinear with any of the starting set
    for i in range(n,len(a)):
        try:
            a[n] = VectorProduct(a[:n])
        except:
            raise Exception,"VectorProduct failed, singular input?"
        n += 1
    # use the G-S algorithm to compute a complete set of unit vectors orthogonal
    # to the first in array
    gs = GramSchmidtOrtho(a) 
    n = nvars
    # now replace the created vectors with the ones least correlated to the
    # initial set
    vlist = [v for v in gs[1:]] # drop the first row
    for j in range(n,len(a)):
        minavg = None
        droplist = []
        for k in range(len(vlist)):
            v = vlist[k]
            avgcor = 0
            for i in range(n):
                cor = np.dot(a[i],vlist[k])/(np.linalg.norm(a[i]) * np.linalg.norm(vlist[k]) )
                if np.allclose(cor, 1.0):
                    droplist.append(k) # found a vector co-linear to one we already have
                                       # don't need it again
                    #vlist.pop(k) 
                    break 
                avgcor += cor
            else:
                if minavg == None:
                    minavg = abs(avgcor/n)
                    minnum = k
                elif abs(avgcor/n) < minavg:
                    minavg = abs(avgcor/n)
                    minnum = k
        if minavg == None:
            raise Exception,"Failed to find a non-colinear G-S vector for row %d. Should not happen!" % n
        a[j] = vlist[minnum]
        droplist.append(minnum) # don't need to use this vector again
        for i in sorted(droplist,reverse=True):
            vlist.pop(i) 
        n += 1
    return a


def MakeArray(constrDict, group, varlist):
    """Convert the multipliers in a constraint group to an array of
    coefficients and place in a square numpy array, adding empty rows as
    needed.

    Throws an exception if some sanity checks fail.
    """
    # do some error checks
    if len(varlist) < len(group): # too many relationships -- no can do
        raise Exception,'The number of relationships (%d) exceeds the number of parameters (%d):\n\t%s\n\t%s'% (
            len(group),len(varlist),group,varlist)
    # put the coefficients into an array
    multarr = np.zeros(2*[len(varlist),])
    i = 0
    for cnum in group:
        cdict = constrDict[cnum]
        j = 0
        for var in varlist:
            m = cdict.get(var)
            if m != None:
                multarr[i,j] = m
            j += 1
        i += 1
    return multarr

def GramSchmidtOrtho(arrayin):
    '''Use the Gram-Schmidt process (http://en.wikipedia.org/wiki/Gram-Schmidt) to
    find orthonormal unit vectors relative to first row.
    input: 
       arrayin: a 2-D non-signular square array
    returns:
       a orthonormal set of unit vectors as a square array
    '''
    def proj(a,b):
        'Projection operator'
        return a*(np.dot(a,b)/np.dot(a,a))
    a = arrayin.copy()
    for j in range (len(a)):
        for i in range(j):
            a[j] = a[j] - proj(a[i],a[j])
        if np.allclose(np.linalg.norm(a[j]),0.0):
            raise Exception,"Singular input to GramSchmidtOrtho"
        a[j] = a[j]/np.linalg.norm(a[j])
    return a

def VectorProduct(vl):
    '''Evaluate the n-dimensional "cross" product of the list of vectors in vl
    vl can be any length. The length of each vector is arbitrary, but
    all must be the same length. See http://en.wikipedia.org/wiki/Levi-Civita_symbol

    This appears to return a vector perpendicular to the supplied set.

    Throws an exception if a vector can not be obtained because the input set of
    vectors is already complete or contains any redundancies.
    
    Uses LeviCitvitaMult recursively to obtain all permutations of vector elements
    '''
    i = 0
    newvec = []
    for e in vl[0]:
        i += 1
        tl = [([i,],1),]
        seq = LeviCitvitaMult(tl ,vl)
        tsum = 0
        for terms,prod in seq:
            tsum += EvalLCterm(terms) * prod
        newvec.append(tsum)
    if np.allclose(newvec,0.0):
        raise Exception,"VectorProduct failed, singular or already complete input"
    return newvec

def LeviCitvitaMult(tin,vl):
    '''Recursion formula to compute all permutations of vector element products
    The first term in each list is the indices of the Levi-Civita symbol, note
    that if an index is repeated, the value is zero, so the element is not included.
    The second term is the product of the vector elements
    '''
    v = vl[0]
    vl1 = vl[1:]        

    j = 0
    tl = []
    for e in vl[0]:
        j += 1
        for ind,vals in tin:
            if j in ind: continue
            tl.append((ind+[j,],vals*e))
    if len(vl1):
        return LeviCitvitaMult(tl,vl1)
    else:
        return tl

def EvalLCterm(term):
    '''Evaluate the Levi-Civita symbol Epsilon(i,j,k,...)'''
    p = 1
    for i in range(len(term)):
        for j in range(i+1,len(term)):
            p *= (term[j] - term[i])
            if p == 0: return 0
    return p/abs(p)


if __name__ == "__main__":
    import sys
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
    varyList = ['0::A0', '0::AUiso:0', '0::Afrac:1', '0::Afrac:2', '0::Afrac:3', '0::Afrac:4', '0::dAx:5', '0::dAy:5', '0::dAz:5', '0::AUiso:5', ':0:Back:0', ':0:Back:1', ':0:Back:2', ':0:Back:3', ':0:Back:4', ':0:Back:5', ':0:Back:6', ':0:Back:7', ':0:Back:8', ':0:Back:9', ':0:Back:10', ':0:Back:11', ':0:U', ':0:V', ':0:W', ':0:X', ':0:Y', ':0:Scale', ':0:DisplaceX', ':0:DisplaceY']
    constrDict = [
        {'0::Afrac:4': 24.0, '0::Afrac:1': 16.0, '0::Afrac:3': 24.0, '0::Afrac:2': 16.0},
        {'0::Afrac:1': 1.0, '0::Afrac:2': 1.0},
        {'0::Afrac:4': 1.0, '0::Afrac:3': 1.0}]
    fixedList = ['40.0', '1.0', '1.0']

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
    dMdv = len(varylist)*[0]
    deriv = {}
    for i,v in enumerate(parmdict.keys()): deriv[v]=i
    Dict2Deriv(varylist,deriv,dMdv)
