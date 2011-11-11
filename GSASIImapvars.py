# $Date: 2011-01-07 13:27:24 -0600 (Fri, 07 Jan 2011) $
# $Author: toby $
# $Revision: 230 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIImapvars.py $
# $Id: GSASIImapvars.py 230 2011-01-07 19:27:24Z vondreele $
########### SVN repository information ###################
"""Module that implements algebraic contraints, parameter redefinition
and parameter simplification contraints.

Parameter redefinition is done by creating one or more relationships
between a set of parameters
   Mx1 * Px + My1 * Py +...
   Mx2 * Px + Mz2 * Pz + ...
where Pj is a parameter name and Mjk is a constant.

Relations are typically supplied as input to InputParse, where either
of two keywords can be attached to a relationship:
 * VARY which means the generated parameter name will be included in
  the list of variables to be refined (varyList)
 * VARYFREE which means that all generated parameter names for a
  group, including the ones for generated relationships (see below)
  will be included in the list of variables to be refined (varyList)
  The case of VARY and VARYFREE is immaterial.

Relations can also be supplied in the form of an equation:
  nx1 * Px + ny1 * Py +... = C1
where Cn is a constant. These equations define an algebraic
constrant. The keyword VARY makes no sense when used with a constraint
equation, but VARYFREE can be used.  #RVD??? is VARYFREE required or default??

Unconstrained relations describe a new, independent, parameter, which
is defined in terms of dependent parameters that are defined in the
Model, while constrained relations effectively reduce the complexity
of the Model by removing a degree of freedom.

Relationships are grouped so that a set of dependent parameters appear
in only one group (done in routine GroupConstraints.) Note that if a
group contains relationships/equations that involve N dependent
parameters, there must exist N-C new parameters, where C is the number
of contraint equations in the group. Routine GenerateConstraints takes
the output from InputParse and GroupConstraints generates the
"missing" relationships and saves that information in the module's
global variables. Each generated parameter is named
sequentially using paramPrefix.

Parameter redefinition is done by equating one (independent) parameter to several
(now dependent) parameters, in algebraic form:
  P1 = n2 * P2 = n3 * P3 ...
Each equality in the relationship reduces the complexity of the model
by one potentially variable parameter. Input is provided to routine
StoreEquivalence in the form of an independent parameter and a list of
dependent parameters, optionally with a multiplier.

Note that none of the dependent parameters in any constraint or
reformulation can be refined (see dependentParmList, below).

External Routines:
   To define a set of constrained and unconstrained relations, one
     calls InputParse, GroupConstraints and GenerateConstraints. It
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
     problematic. This may not work properly:
        StoreEquivalence('x',('y',))
        StoreEquivalence('y',('z',))

   To show a summary of the parameter remapping, one calls
      VarRemapShow

   To determine values for the parameters created in this module, one
      calls Map2Dict. This will not apply contraints.

   To take values from the new independent parameters and constraints,
      one calls Dict2Map. This will apply contraints.

Global Variables:
   dependentParmList: contains a list by group of lists of
     parameters used in the group. Note that parameters listed in
     dependentParmList should not be refined as they will not affect
     the model
   indParmList: a list of groups of Independent parameters defined in
     each group. This contains both parameters used in parameter
     redefinitions as well as names of generated new parameters.
   
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
import re
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

# compile regular expressions used for parsing input
rex_mult = re.compile('[+-]?[0-9.]+[eE]?[+-]?[0-9]*')
rex_star = re.compile('\s*\*\s*')
rex_var = re.compile('\w+')
rex_plusminus = re.compile('\s*[+-=]\s*')
rex_vary = re.compile('\s*Vary\s*', re.IGNORECASE)
rex_varyfree = re.compile('(.*)\s*VaryFree\s*', re.IGNORECASE)

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

def InputParse(mapList):
    '''Converts a set relationships used to remap parameters into composite
    parameters or to one or more constants.
        
    Input:
      mapList: a list of character strings where each string defines
    a relationship in form:
    ('<const11>*<prm1> [+-] <const12>*<prm2> [+-] ... = <value> [VaryFree]',
    '<const21>*<prm1> [+-] <const22>*<prm2> [+-] ... [Vary/VaryFree]',
    ...)
    these define either a constraint or a new independent parameter,
    where <constXX> is a constant containing an optional sign, digits,
    optionally a decimal point and/or an exponent, prefixed by e or E,
    and <prmN> is the name of a parameter defined in the Model.
    Parameters can be included in any order. Multiplying a parameter
    by 0 causes that parameter to be included in a group of
    relationships.  Note that if the same parameter is included twice
    in a relationship/equation, the coefficients are added.
    
    When an equality is defined, a constant, <value>, is
    included. This then describes a constraint equation.
    
    Vary is an optional keyword that indicates the indicated remapped
    parameter should be varied (note, case insensitive). Should be the
    last item on the line. Makes no sense with an equation.

    VaryFree is an optional keyword that indicates all possible
    remapped parameters should be varied (note, case
    insensitive). Makes most sense with a constraint equation.  Should
    be the last item on the line.

    returns
      constrDict: a list with a dict for each item in mapList that
        defines the relationship. The keys are parameter names and the
        values are the multiplier for the parameter name
      constrFlag: a list for each item in mapList with a list contains
        'Vary' and/or 'VaryFree'
      fixedList: a list for each item in mapList. Contains the value
        (as a string) for each contraint equation or None for a
        constraint relationship.
    '''
    i = 0
    constrDict = []
    fixedList = []
    constrFlag = []
    if debug: print 50*'-','\n(Input)'
    for line in mapList:
        line = line.strip()
        if len(line) == 0: continue # empty lines are ignored
        constrDict.append({})
        constrFlag.append([])
        fixedList.append(None)
        i += 1
        fixedval = None
        if debug: print '\t',line
        try: 
            j = 0
            sign = 1.0
            while len(line):
                j += 1
                m = sign * float(rex_mult.match(line).group())
                line = line[rex_mult.match(line).end():]
                j += 1
                line = line[rex_star.match(line).end():]
                j += 1
                v = rex_var.match(line).group()
                line = line[rex_var.match(line).end():]
                #print m,'times',v
                #if v not in varlist: varlist.append(v)
                if constrDict[i-1].get(v) == None:
                    constrDict[i-1][v] = m
                else:
                    constrDict[i-1][v] += m
                if len(line.strip()) == 0: break
                j += 1
                # advance to next separator (+-=)
                pm = rex_plusminus.match(line)
                if pm != None:
                    line = line[pm.end():]
                    pm = pm.group()
                else:
                    pm = ''
                if pm.strip() == '+':
                    sign = 1.0
                elif pm.strip() == '-':
                    sign = -1.0
                elif pm.strip() == '=':
                    # found a fixed value, also check for a VaryFree flag
                    if rex_varyfree.match(line):
                        constrFlag[-1].append('VaryFree')
                        #fixedval = float(rex_varyfree.split(line)[1])
                        fixedval = rex_varyfree.split(line)[1].strip()
                    else:
                        #fixedval = float(line.strip())
                        fixedval = line.strip()
                    fixedList[-1] = fixedval
                    line = ""
                elif rex_varyfree.match(line):
                    constrFlag[-1].append('VaryFree')
                    line = line[rex_varyfree.match(line).end():]
                elif rex_vary.match(line):
                    constrFlag[-1].append('Vary')
                    line = line[rex_vary.match(line).end():]
                else:
                    # something else is on the line, but not a keyword
                    raise Exception
        except SyntaxError:
            if debug: print 'Error in line',i,'token',j,'@','"'+line+'"'
            raise Exception,'Error in line %d token %d, beginning with %s'% (i,j, line)

    if debug: # on debug, show what is parsed in a semi-readable
        print 50*'-','\n(parsed relationship/equation & flag)'
        for i in range(len(constrDict)):
            flags = ''
            for f in constrFlag[i]:
                if flags != '':
                    flags = flags + ' + ' + f
                else:
                    flags = f
            if fixedList[i] is None:
                print '#',i+1,constrDict[i],'\t',constrFlag[i]
            else:
                print '#',i+1,constrDict[i],'=',fixedList[i],'\t',constrFlag[i]

    return constrDict,constrFlag,fixedList

def GroupConstraints(constrDict):
    """divide the constraints into groups that share no parameters.
    Input
       constrDict: a list of dicts defining relationships/constraints
       (see InputParse)
    Returns two lists of lists:
      a list of relationship list indices for each group and
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

def GenerateConstraints(groups,parmlist,varyList,constrDict,constrFlag,fixedList):
    '''Takes a list of relationship entries comprising a group of constraints and
    builds the relationship lists and their inverse and stores them in global variables
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,consNum
    for group,varlist in zip(groups,parmlist):
        VaryFree = False
        for row in group:
            if 'VaryFree' in constrFlag[row]: VaryFree = True
        arr = MakeArray(constrDict, group, varlist)
        constrArr = FillInMissingRelations(arr,len(group))
        mapvar = []
        group = group[:]
        for i in range(len(varlist)):
            vary = VaryFree
            if len(group) > 0:
                rel = group.pop(0)
                fixedval = fixedList[rel]
                if 'Vary' in constrFlag[rel]: vary = True
            else:
                fixedval = None
            if fixedval is None:
                varname = paramPrefix + str(consNum)
                mapvar.append(varname)
                consNum += 1
                if vary: varyList.append(varname)
            else:
                mapvar.append(fixedval)
        dependentParmList.append(varlist)
        arrayList.append(constrArr)
        invarrayList.append(np.linalg.inv(constrArr))
        indParmList.append(mapvar)
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

def VarRemapShow(varyList):
    '''List out the saved relationships.
    Returns a string containing the details.
    '''
    s = 'Mapping relations:\n'
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    for varlist,mapvars,multarr in zip(dependentParmList,indParmList,arrayList):
        i = 0
        for mv in mapvars:
            if multarr is None:
                s += '  ' + str(mv) + ' defines parameter(s): '
                j = 0
                for v in varlist:
                    if j > 0: s += '  & '
                    j += 1
                    s += str(v)
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
            i += 1
    s += 'Inverse mapping relations:\n'
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        i = 0
        for mv in varlist:
            s += '  %s = ' % mv
            j = 0
            for m,v in zip(invmultarr[i,:],mapvars):
                if m == 0: continue
                if j > 0: s += ' + '
                j += 1
                s += '(%s * %s)' % (m,v)
            s += '\n'
            i += 1
    return s

def Map2Dict(parmDict,varyList):
    '''Create (or update) the Independent Parameters from the original
    set of Parameters

    This should be done once, before any variable refinement is done
    to complete the parameter dictionary with the Independent Parameters
    '''
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
    # overwrite dict with constraints - why not parmDict.update(fixDict)?
    parmDict.update(fixedDict)

def Dict2Map(parmDict,varyList):
    #I think this needs fixing to update parmDict with new values 
    #   from the constraints based on what is in varyList - RVD 
    '''Convert the remapped values back to the original parameters
    
    This should be done to apply constraints to parameter values (use
    Map2Dict first). It should be done as part of the Model function
    evaluation, before any computation is done
    '''
    global dependentParmList,arrayList,invarrayList,indParmList,fixedDict
    # reset fixed values (should not be needed, but very quick) 
    # - this seems to update parmDict with {'0':0.0} & {'1.0':1.0} - probably not what was intended
    parmDict.update(fixedDict)
    for varlist,mapvars,invmultarr in zip(dependentParmList,indParmList,invarrayList):
        if invmultarr is None: continue
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
    if len(varlist) == 1: # one to one mapping -- something is probably wrong
        raise Exception,'Why remap a single parameter? (%s)'% (varlist[0])
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
    #remap = MapParameters() # create an object (perhaps better as a module)
    #remap.debug = 1
    #debug = 1
    parmdict = {}
    StoreEquivalence('2::atomx:3',('2::atomy:3', ('2::atomz:3',2,), ))
    print VarRemapShow()

    parmdict = {
        '2::atomx:3':.1 ,
        '2::atomy:3':.2 , # conflicts with constraint
        '2::atomz:3':.3 ,
    }

    mapList = [
        "1*p1 + 2e0*p2 - 1.0*p3",
        "1*p4 + 1*p7",
        "1*p1+2*p2-3.0*p2 VARY",
        "1*p21 + 0*p22 - 0*p23 + 0*p24 varyfree",
        "1*p5 + 2*p6 = 1 varyfree",
        "-1 * p6 + 1*p5",
        "-10e-1 * p1 - -2*p2 + 3.0*p4",
        ]
    constrDict,constrFlag,fixedList = InputParse(mapList)
    print constrDict
    print constrFlag
    print fixedList
    groups,parmlist = GroupConstraints(constrDict)
    GenerateConstraints(groups,parmlist,constrDict,constrFlag,fixedList)
    print VarRemapShow()
    parmdict.update( {'p1':1,'p2':2,'p3':3,'p4':4,
                      'p6':6,'p5':5,  # conflicts with constraint
                      'p7':7,
                      "p21":.1,"p22":.2,"p23":.3,"p24":.4,
                      })
    #print 'parmdict start',parmdict
    before = parmdict.copy()
    Map2Dict(parmdict,[])
    print 'parmdict before and after Map2Dict'
    print '  key / before / after'
    for key in sorted(parmdict.keys()):
        print '  '+key,'\t',before.get(key),'\t',parmdict[key]

    before = parmdict.copy()
    Dict2Map(parmdict,[])
    print 'after Dict2Map'
    print '  key / before / after'
    for key in sorted(parmdict.keys()):
        print '  '+key,'\t',before.get(key),'\t',parmdict[key]
