# -*- coding: utf-8 -*-
#GSASIIobj - data objects for GSAS-II
'''
Classes and routines defined in :mod:`GSASIIobj` follow. 
'''
# Note that documentation for GSASIIobj.py has been moved
# to file docs/source/GSASIIobj.rst

from __future__ import division, print_function
import platform
import re
import random as ran
import sys
import os.path
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import GSASIIpath
import GSASIImath as G2mth
import GSASIIspc as G2spc
import numpy as np

DefaultControls = {
    'deriv type':'analytic Hessian',
    'min dM/M':0.001,'shift factor':1.,'max cyc':3,'F**2':False,'SVDtol':1.e-6,
    'UsrReject':{'minF/sig':0.,'MinExt':0.01,'MaxDF/F':100.,'MaxD':500.,'MinD':0.05},
    'Copy2Next':False,'Reverse Seq':False,'HatomFix':False,
    'Author':'no name','newLeBail':False,
    'FreePrm1':'Sample humidity (%)',
    'FreePrm2':'Sample voltage (V)',
    'FreePrm3':'Applied load (MN)',
    'ShowCell':False,
    }
'''Values to be used as defaults for the initial contents of the ``Controls``
data tree item.
'''

restraintNames = [['Bond','Bonds'],['Angle','Angles'],['Plane','Planes'],
                  ['Chiral','Volumes'],['Torsion','Torsions'],['Rama','Ramas'],
                  ['ChemComp','Sites'],['Texture','HKLs'],['Moments','Moments'],
                  ['General','General']]
'''Names of restraint keys for the restraint dict and the location of 
the restraints in each dict
'''
    
def StripUnicode(string,subs='.'):
    '''Strip non-ASCII characters from strings

    :param str string: string to strip Unicode characters from
    :param str subs: character(s) to place into string in place of each
      Unicode character. Defaults to '.'

    :returns: a new string with only ASCII characters
    '''
    s = ''
    for c in string:
        if ord(c) < 128:
            s += c
        else:
            s += subs
    return s
#    return s.encode('ascii','replace')

def MakeUniqueLabel(lbl,labellist):
    '''Make sure that every a label is unique against a list by adding
    digits at the end until it is not found in list.

    :param str lbl: the input label
    :param list labellist: the labels that have already been encountered
    :returns: lbl if not found in labellist or lbl with ``_1-9`` (or
      ``_10-99``, etc.) appended at the end
    '''
    lbl = StripUnicode(lbl.strip(),'_')
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

#VarDesc = {}  # removed 1/30/19 BHT as no longer needed (I think)
#''' This dictionary lists descriptions for GSAS-II variables,
#as set in :func:`CompileVarDesc`. See that function for a description
#for how keys and values are written.
#'''

reVarDesc = {}
''' This dictionary lists descriptions for GSAS-II variables where
keys are compiled regular expressions that will match the name portion
of a parameter name. Initialized in :func:`CompileVarDesc`.
'''

reVarStep = {}
''' This dictionary lists the preferred step size for numerical 
derivative computation w/r to a GSAS-II variable. Keys are compiled 
regular expressions and values are the step size for that parameter. 
Initialized in :func:`CompileVarDesc`.
'''
# create a default space group object for P1; N.B. fails when building documentation
try:
    P1SGData = G2spc.SpcGroup('P 1')[1] # data structure for default space group
except:
    pass

def GetPhaseNames(fl):
    ''' Returns a list of phase names found under 'Phases' in GSASII gpx file
    NB: there is another one of these in GSASIIstrIO.py that uses the gpx filename

    :param file fl: opened .gpx file
    :return: list of phase names
    '''
    PhaseNames = []
    while True:
        try:
            data = cPickle.load(fl)
        except EOFError:
            break
        datum = data[0]
        if 'Phases' == datum[0]:
            for datus in data[1:]:
                PhaseNames.append(datus[0])
    fl.seek(0)          #reposition file
    return PhaseNames

def SetNewPhase(Name='New Phase',SGData=None,cell=None,Super=None):
    '''Create a new phase dict with default values for various parameters

    :param str Name: Name for new Phase

    :param dict SGData: space group data from :func:`GSASIIspc:SpcGroup`;
      defaults to data for P 1

    :param list cell: unit cell parameter list; defaults to
      [1.0,1.0,1.0,90.,90,90.,1.]

    '''
    if SGData is None: SGData = P1SGData
    if cell is None: cell=[1.0,1.0,1.0,90.,90.,90.,1.]
    phaseData = {
        'ranId':ran.randint(0,sys.maxsize),
        'General':{
            'Name':Name,
            'Type':'nuclear',
            'Modulated':False,
            'AtomPtrs':[3,1,7,9],
            'SGData':SGData,
            'Cell':[False,]+cell,
            'Pawley dmin':1.0,
            'Data plot type':'None',
            'SH Texture':{
                'Order':0,
                'Model':'cylindrical',
                'Sample omega':[False,0.0],
                'Sample chi':[False,0.0],
                'Sample phi':[False,0.0],
                'SH Coeff':[False,{}],
                'SHShow':False,
                'PFhkl':[0,0,1],
                'PFxyz':[0,0,1],
                'PlotType':'Pole figure',
                'Penalty':[['',],0.1,False,1.0]}},
        'Atoms':[],
        'Drawing':{},
        'Histograms':{},
        'Pawley ref':[],
        'RBModels':{},
        }
    if Super and Super.get('Use',False):
        phaseData['General'].update({'Modulated':True,'Super':True,'SuperSg':Super['ssSymb']})
        phaseData['General']['SSGData'] = G2spc.SSpcGroup(SGData,Super['ssSymb'])[1]
        phaseData['General']['SuperVec'] = [Super['ModVec'],False,Super['maxH']]

    return phaseData

def ReadCIF(URLorFile):
    '''Open a CIF, which may be specified as a file name or as a URL using PyCifRW
    (from James Hester).
    The open routine gets confused with DOS names that begin with a letter and colon
    "C:\\dir\" so this routine will try to open the passed name as a file and if that
    fails, try it as a URL

    :param str URLorFile: string containing a URL or a file name. Code will try first
      to open it as a file and then as a URL.

    :returns: a PyCifRW CIF object.
    '''
    import CifFile as cif # PyCifRW from James Hester

    # alternate approach:
    #import urllib
    #ciffile = 'file:'+urllib.pathname2url(filename)

    try:
        fp = open(URLorFile,'r')
        cf = cif.ReadCif(fp)
        fp.close()
        return cf
    except IOError:
        return cif.ReadCif(URLorFile)

def TestIndexAll():
    '''Test if :func:`IndexAllIds` has been called to index all phases and 
    histograms (this is needed before :func:`G2VarObj` can be used. 

    :returns: Returns True if indexing is needed.
    '''
    if PhaseIdLookup or AtomIdLookup or HistIdLookup:
        return False
    return True
        
def IndexAllIds(Histograms,Phases):
    '''Scan through the used phases & histograms and create an index
    to the random numbers of phases, histograms and atoms. While doing this,
    confirm that assigned random numbers are unique -- just in case lightning
    strikes twice in the same place.

    Note: this code assumes that the atom random Id (ranId) is the last
    element each atom record.

    This is called when phases & histograms are looked up 
    in these places (only): 
        
     * :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` (which loads the histograms and phases from a GPX file),
     * :meth:`~GSASIIdataGUI.GSASII.GetUsedHistogramsAndPhasesfromTree` (which does the same thing but from the data tree.) 
     * :meth:`~GSASIIdataGUI.GSASII.OnFileClose` (clears out an old project)
    
    Note that globals :data:`PhaseIdLookup` and :data:`PhaseRanIdLookup` are 
    also set in :func:`AddPhase2Index` to temporarily assign a phase number
    as a phase is being imported.
 
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
            Phases[ph]['ranId'] = ranId = ran.randint(0,sys.maxsize)
        pId = str(Phases[ph]['pId'])
        PhaseIdLookup[pId] = (ph,ranId)
        PhaseRanIdLookup[ranId] = pId
        shortname = ph  #[:10]
        while shortname in ShortPhaseNames.values():
            shortname = ph[:8] + ' ('+ pId + ')'
        ShortPhaseNames[pId] = shortname
        AtomIdLookup[pId] = {}
        AtomRanIdLookup[pId] = {}
        for iatm,at in enumerate(Phases[ph]['Atoms']):
            ranId = at[cia+8]
            while ranId in AtomRanIdLookup[pId]: # check for dups
                print ("\n\n*** Phase "+str(ph)+" atom "+str(iatm)+" has repeated ranId. Fixing.\n")
                at[cia+8] = ranId = ran.randint(0,sys.maxsize)
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
            Histograms[hist]['ranId'] = ranId = ran.randint(0,sys.maxsize)
        hId = str(Histograms[hist]['hId'])
        HistIdLookup[hId] = (hist,ranId)
        HistRanIdLookup[ranId] = hId
        shortname = hist[:15]
        while shortname in ShortHistNames.values():
            shortname = hist[:11] + ' ('+ hId + ')'
        ShortHistNames[hId] = shortname

def AddPhase2Index(rdObj,filename):
    '''Add a phase to the index during reading
    Used where constraints are generated during import (ISODISTORT CIFs)        
    '''
    ranId = rdObj.Phase['ranId']
    ph = 'from  '+filename #  phase is not named yet
    if ranId in PhaseRanIdLookup: return
    maxph = -1
    for r in PhaseRanIdLookup:
        maxph = max(maxph,int(PhaseRanIdLookup[r]))
    PhaseRanIdLookup[ranId] = pId = str(maxph + 1)
    PhaseIdLookup[pId] = (ph,ranId)
    shortname = 'from '+ os.path.splitext((os.path.split(filename))[1])[0]
    while shortname in ShortPhaseNames.values():
        shortname = ph[:8] + ' ('+ pId + ')'
    ShortPhaseNames[pId] = shortname
    AtomIdLookup[pId] = {}
    AtomRanIdLookup[pId] = {}
    for iatm,at in enumerate(rdObj.Phase['Atoms']):
        ranId = at[-1]
        while ranId in AtomRanIdLookup[pId]: # check for dups
            print ("\n\n*** Phase "+str(ph)+" atom "+str(iatm)+" has repeated ranId. Fixing.\n")
            at[-1] = ranId = ran.randint(0,sys.maxsize)
        AtomRanIdLookup[pId][ranId] = str(iatm)
        #if Phases[ph]['General']['Type'] == 'macromolecular':
        #    label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
        #else:
        #    label = at[ct-1]
        label = at[0]
        AtomIdLookup[pId][str(iatm)] = (label,ranId)

def LookupAtomId(pId,ranId):
    '''Get the atom number from a phase and atom random Id

    :param int/str pId: the sequential number of the phase
    :param int ranId: the random Id assigned to an atom

    :returns: the index number of the atom (str)
    '''
    if not AtomRanIdLookup:
        raise Exception('Error: LookupAtomId called before IndexAllIds was run')
    if pId is None or pId == '':
        raise KeyError('Error: phase is invalid (None or blank)')
    pId = str(pId)
    if pId not in AtomRanIdLookup:
        raise KeyError('Error: LookupAtomId does not have phase '+pId)
    if ranId not in AtomRanIdLookup[pId]:
        raise KeyError('Error: LookupAtomId, ranId '+str(ranId)+' not in AtomRanIdLookup['+pId+']')
    return AtomRanIdLookup[pId][ranId]

def LookupAtomLabel(pId,index):
    '''Get the atom label from a phase and atom index number

    :param int/str pId: the sequential number of the phase
    :param int index: the index of the atom in the list of atoms

    :returns: the label for the atom (str) and the random Id of the atom (int)
    '''
    if not AtomIdLookup:
        raise Exception('Error: LookupAtomLabel called before IndexAllIds was run')
    if pId is None or pId == '':
        raise KeyError('Error: phase is invalid (None or blank)')
    pId = str(pId)
    if pId not in AtomIdLookup:
        raise KeyError('Error: LookupAtomLabel does not have phase '+pId)
    if index not in AtomIdLookup[pId]:
        raise KeyError('Error: LookupAtomLabel, ranId '+str(index)+' not in AtomRanIdLookup['+pId+']')
    return AtomIdLookup[pId][index]

def LookupPhaseId(ranId):
    '''Get the phase number and name from a phase random Id

    :param int ranId: the random Id assigned to a phase
    :returns: the sequential Id (pId) number for the phase (str)
    '''
    if not PhaseRanIdLookup:
        raise Exception('Error: LookupPhaseId called before IndexAllIds was run')
    if ranId not in PhaseRanIdLookup:
        raise KeyError('Error: LookupPhaseId does not have ranId '+str(ranId))
    return PhaseRanIdLookup[ranId]

def LookupPhaseName(pId):
    '''Get the phase number and name from a phase Id

    :param int/str pId: the sequential assigned to a phase
    :returns:  (phase,ranId) where phase is the name of the phase (str)
      and ranId is the random # id for the phase (int)
    '''
    if not PhaseIdLookup:
        raise Exception('Error: LookupPhaseName called before IndexAllIds was run')
    if pId is None or pId == '':
        raise KeyError('Error: phase is invalid (None or blank)')
    pId = str(pId)
    if pId not in PhaseIdLookup:
        raise KeyError('Error: LookupPhaseName does not have index '+pId)
    return PhaseIdLookup[pId]

def LookupHistId(ranId):
    '''Get the histogram number and name from a histogram random Id

    :param int ranId: the random Id assigned to a histogram
    :returns: the sequential Id (hId) number for the histogram (str)
    '''
    if not HistRanIdLookup:
        raise Exception('Error: LookupHistId called before IndexAllIds was run')
    if ranId not in HistRanIdLookup:
        raise KeyError('Error: LookupHistId does not have ranId '+str(ranId))
    return HistRanIdLookup[ranId]

def LookupHistName(hId):
    '''Get the histogram number and name from a histogram Id

    :param int/str hId: the sequential assigned to a histogram
    :returns:  (hist,ranId) where hist is the name of the histogram (str)
      and ranId is the random # id for the histogram (int)
    '''
    if not HistIdLookup:
        raise Exception('Error: LookupHistName called before IndexAllIds was run')
    if hId is None or hId == '':
        raise KeyError('Error: histogram is invalid (None or blank)')
    hId = str(hId)
    if hId not in HistIdLookup:
        raise KeyError('Error: LookupHistName does not have index '+hId)
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

    # special handling for parameter names without a colon
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
            lbl = 'Seq. ref.'
        else:
            lbl = ShortPhaseNames.get(l[0],'? #'+str(l[0]))
        if l[1] == '*':
            hlbl = 'Seq. ref.'
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
        elif l[4] is not None: # rigid body parameter or modulation parm
            lbl = ShortPhaseNames.get(l[0],'phase?')
            if 'RB' in l[2]:    #rigid body parm
                s = "RB body #"+str(l[3])+" (type "+str(l[4])+") in "+str(lbl) + ','
            else: #modulation parm
                s = 'Atom %s wave %s in %s'%(LookupAtomLabel(l[0],l[3])[0],l[4],lbl)
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
    '''Set the values in the variable lookup tables 
    (:attr:`reVarDesc` and :attr:`reVarStep`).
    This is called in :func:`getDescr` and :func:`getVarStep` so this
    initialization is always done before use. These variables are 
    also used in script `makeVarTbl.py` which creates the table in section 3.2
    of the Sphinx docs (:ref:`VarNames_table`).

    Note that keys may contain regular expressions, where '[xyz]'
    matches 'x' 'y' or 'z' (equivalently '[x-z]' describes this as range 
    of values). '.*' matches any string. For example::

    'AUiso':'Atomic isotropic displacement parameter',

    will match variable ``'p::AUiso:a'``.
    If parentheses are used in the key, the contents of those parentheses can be
    used in the value, such as::

    'AU([123][123])':'Atomic anisotropic displacement parameter U\\1',

    will match ``AU11``, ``AU23``,... and `U11`, `U23` etc will be displayed
    in the value when used.

    '''
    if reVarDesc: return # already done
    for key,value in {
        # derived or other sequential vars
        '([abc])$' : 'Lattice parameter, \\1, from Ai and Djk', # N.B. '$' prevents match if any characters follow
        u'\u03B1' : u'Lattice parameter, \u03B1, computed from both Ai and Djk',
        u'\u03B2' : u'Lattice parameter, \u03B2, computed from both Ai and Djk',
        u'\u03B3' : u'Lattice parameter, \u03B3, computed from both Ai and Djk',
        # ambiguous, alas:
        'Scale' : 'Phase fraction (as p:h:Scale) or Histogram scale factor (as :h:Scale)',
        # Phase vars (p::<var>)
        'A([0-5])' : ('Reciprocal metric tensor component \\1',1e-5),
        '([vV]ol)' : 'Unit cell volume', # probably an error that both upper and lower case are used
        # Atom vars (p::<var>:a)
        'dA([xyz])$' : ('Refined change to atomic coordinate, \\1',1e-6),
        'A([xyz])$' : 'Fractional atomic coordinate, \\1',
        'AUiso':('Atomic isotropic displacement parameter',1e-4),
        'AU([123][123])':('Atomic anisotropic displacement parameter U\\1',1e-4),
        'Afrac': ('Atomic site fraction parameter',1e-5),
        'Amul': 'Atomic site multiplicity value',
        'AM([xyz])$' : 'Atomic magnetic moment parameter, \\1',
        # Atom deformation parameters
        'Akappa([0-6])'  : ' Atomic orbital softness for orbital, \\1',
        'ANe([01])' : ' Atomic <j0> orbital population for orbital, \\1',
        'AD\\([0-6],[0-6]\\)([0-6])' : ' Atomic sp. harm. coeff for orbital, \\1',
        'AD\\([0-6],-[0-6]\\)([0-6])' : ' Atomic sp. harm. coeff for orbital, \\1',     #need both!
        # Hist (:h:<var>) & Phase (HAP) vars (p:h:<var>)
        'Back(.*)': 'Background term #\\1',
        'BkPkint;(.*)':'Background peak #\\1 intensity',
        'BkPkpos;(.*)':'Background peak #\\1 position',
        'BkPksig;(.*)':'Background peak #\\1 Gaussian width',
        'BkPkgam;(.*)':'Background peak #\\1 Cauchy width',
#        'Back File' : 'Background file name',
        'BF mult' : 'Background file multiplier',
        'Bab([AU])': 'Babinet solvent scattering coef. \\1',
        'D([123][123])' : 'Anisotropic strain coef. \\1',
        'Extinction' : 'Extinction coef.',
        'MD' : 'March-Dollase coef.',
        'Mustrain;.*' : 'Microstrain coefficient (delta Q/Q x 10**6)',
        'Size;.*' : 'Crystallite size value (in microns)',
        'eA$' : 'Cubic mustrain value',
        'Ep$' : 'Primary extinction',
        'Es$' : 'Secondary type II extinction',
        'Eg$' : 'Secondary type I extinction',
        'Flack' : 'Flack parameter',
        'TwinFr' : 'Twin fraction',
        'Layer Disp'  : 'Layer displacement along beam',
        #Histogram vars (:h:<var>)
        'Absorption' : 'Absorption coef.',
        'LayerDisp'  : 'Bragg-Brentano Layer displacement',
        'Displace([XY])' : ('Debye-Scherrer sample displacement \\1',0.1),
        'Lam' : ('Wavelength',1e-6),
        'I\\(L2\\)\\/I\\(L1\\)' : ('Ka2/Ka1 intensity ratio',0.001),
        'Polariz.' : ('Polarization correction',1e-3),
        'SH/L' : ('FCJ peak asymmetry correction',1e-4),
        '([UVW])$' : ('Gaussian instrument broadening \\1',1e-5),
        '([XYZ])$' : ('Cauchy instrument broadening \\1',1e-5),
        'Zero' : 'Debye-Scherrer zero correction',
        'Shift' : 'Bragg-Brentano sample displ.',
        'SurfRoughA' : 'Bragg-Brenano surface roughness A',
        'SurfRoughB' : 'Bragg-Brenano surface roughness B',
        'Transparency' : 'Bragg-Brentano sample tranparency',
        'DebyeA' : 'Debye model amplitude',
        'DebyeR' : 'Debye model radius',
        'DebyeU' : 'Debye model Uiso',
        'RBV.*' : 'Vector rigid body parameter',
        'RBVO([aijk])' : 'Vector rigid body orientation parameter \\1',
        'RBVP([xyz])' : 'Vector rigid body \\1 position parameter',
        'RBVf' : 'Vector rigid body site fraction',
        'RBV([TLS])([123AB][123AB])' : 'Residue rigid body group disp. param.',
        'RBVU' : 'Residue rigid body group Uiso param.',
        'RBRO([aijk])' : 'Residue rigid body orientation parameter \\1',
        'RBRP([xyz])' : 'Residue rigid body \\1 position parameter',
        'RBRTr;.*' : 'Residue rigid body torsion parameter',
        'RBRf' : 'Residue rigid body site fraction',
        'RBR([TLS])([123AB][123AB])' : 'Residue rigid body group disp. param.',
        'RBRU' : 'Residue rigid body group Uiso param.',
        'RBSAtNo' : 'Atom number for spinning rigid body',
        'RBSO([aijk])' : 'Spinning rigid body orientation parameter \\1',
        'RBSP([xyz])' : 'Spinning rigid body \\1 position parameter',
        'RBSShRadius' : 'Spinning rigid body shell radius',
        'RBSShC([1-20,1-20])'  : 'Spinning rigid body sph. harmonics term',
        'constr([0-9]*)' : 'Generated degree of freedom from constraint',
        'nv-(.+)' : 'New variable assignment with name \\1',
        # supersymmetry parameters  p::<var>:a:o 'Flen','Fcent'?
        'mV([0-2])$' : 'Modulation vector component \\1',
        'Fsin'  :   'Sin site fraction modulation',
        'Fcos'  :   'Cos site fraction modulation',
        'Fzero'  :   'Crenel function offset',      #may go away
        'Fwid'   :   'Crenel function width',
        'Tmin'   :   'ZigZag/Block min location',
        'Tmax'   :   'ZigZag/Block max location',
        '([XYZ])max': 'ZigZag/Block max value for \\1',
        '([XYZ])sin'  : 'Sin position wave for \\1',
        '([XYZ])cos'  : 'Cos position wave for \\1',
        'U([123][123])sin$' :  'Sin thermal wave for U\\1',
        'U([123][123])cos$' :  'Cos thermal wave for U\\1',
        'M([XYZ])sin$' :  'Sin mag. moment wave for \\1',
        'M([XYZ])cos$' :  'Cos mag. moment wave for \\1',
        # PDF peak parms (l:<var>;l = peak no.)
        'PDFpos'  : 'PDF peak position',
        'PDFmag'  : 'PDF peak magnitude',
        'PDFsig'  : 'PDF peak std. dev.',
        # SASD vars (l:<var>;l = component)
        'Aspect ratio' : 'Particle aspect ratio',
        'Length' : 'Cylinder length',
        'Diameter' : 'Cylinder/disk diameter',
        'Thickness' : 'Disk thickness',
        'Shell thickness' : 'Multiplier to get inner(<1) or outer(>1) sphere radius',
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
        'e([12][12])' : 'strain tensor e\\1',   # strain vars e11, e22, e12
        'Dcalc': 'Calc. d-spacing',
        'Back$': 'background parameter',
        'pos$': 'peak position',
        'int$': 'peak intensity',
        'WgtFrac':'phase weight fraction',
        'alpha':'TOF profile term',
        'alpha-([01])':'Pink profile term',
        'beta-([01q])':'TOF/Pink profile term',
        'sig-([012q])':'TOF profile term',
        'dif([ABC])':'TOF to d-space calibration',
        'C\\([0-9]*,[0-9]*\\)' : 'spherical harmonics preferred orientation coef.',
        'Pressure': 'Pressure level for measurement in MPa',
        'Temperature': 'T value for measurement, K',
        'FreePrm([123])': 'User defined measurement parameter \\1',
        'Gonio. radius': 'Distance from sample to detector, mm',
        }.items():
        # Needs documentation: HAP: LeBail, newLeBail
        # hist: Azimuth, Chi, Omega, Phi, Bank, nDebye, nPeaks
        
        if len(value) == 2:
            #VarDesc[key] = value[0]
            reVarDesc[re.compile(key)] = value[0]
            reVarStep[re.compile(key)] = value[1]
        else:
            #VarDesc[key] = value
            reVarDesc[re.compile(key)] = value

def removeNonRefined(parmList):
    '''Remove items from variable list that are not refined and should not 
    appear as options for constraints

    :param list parmList: a list of strings of form "p:h:VAR:a" where
      VAR is the variable name

    :returns: a list after removing variables where VAR matches a 
      entry in local variable NonRefinedList
    '''
    NonRefinedList = ['Omega','Type','Chi','Phi', 'Azimuth','Gonio. radius',
                          'Lam1','Lam2','Back','Temperature','Pressure',
                          'FreePrm1','FreePrm2','FreePrm3',
                          'Source','nPeaks','LeBail','newLeBail','Bank',
                          'nDebye', #'',
                    ]
    return [prm for prm in parmList if prm.split(':')[2] not in NonRefinedList]
        
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
            try:
                return m.expand(reVarDesc[key])
            except:
                print('Error in key: %s'%key)
    return None

def getVarStep(name,parmDict=None):
    '''Return a step size for computing the derivative of a GSAS-II variable

    :param str name: A complete variable name (with colons, :)
    :param dict parmDict: A dict with parameter values or None (default)

    :returns: a float that should be an appropriate step size, either from 
      the value supplied in :func:`CompileVarDesc` or based on the value for 
      name in parmDict, if supplied. If not found or the value is zero, 
      a default value of 1e-5 is used. If parmDict is None (default) and 
      no value is provided in :func:`CompileVarDesc`, then None is returned.
    '''
    CompileVarDesc() # compile the regular expressions, if needed
    for key in reVarStep:
        m = key.match(name)
        if m:
            return reVarStep[key]
    if parmDict is None: return None
    val = parmDict.get(key,0.0)
    if abs(val) > 0.05:
        return abs(val)/1000.
    else:
        return 1e-5

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
      (such as \\*::var)
    :param list varlist: the list of all variable names used in
      the current project
    :returns: a list of matching GSAS-II variables (may be empty)
    '''
    rexp = re.compile(varname.replace('*','[0-9]+'))
    return sorted([var for var in varlist if rexp.match(var)])

def prmLookup(name,prmDict):
    '''Looks for a parameter in a min/max dictionary, optionally
    considering a wild card for histogram or atom number (use of
    both will never occur at the same time).

    :param name: a GSAS-II parameter name (str, see :func:`getVarDescr` 
      and :func:`CompileVarDesc`) or a :class:`G2VarObj` object. 
    :param dict prmDict: a min/max dictionary, (parmMinDict 
      or parmMaxDict in Controls) where keys are :class:`G2VarObj`
      objects. 
    :returns: Two values, (**matchname**, **value**), are returned where:

       * **matchname** *(str)* is the :class:`G2VarObj` object 
         corresponding to the actual matched name, 
         which could contain a wildcard even if **name** does not; and 
       * **value** *(float)* which contains the parameter limit.
    '''
    for key,value in prmDict.items():
        if str(key) == str(name): return key,value
        if key == name: return key,value
    return None,None
        

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

def SortVariables(varlist):
    '''Sorts variable names in a sensible manner
    '''
    def cvnnums(var):
        v = []
        for i in var.split(':'):
#            if i == '' or i == '*':
#                v.append(-1)
#                continue
            try:
                v.append(int(i))
            except:
                v.append(-1)
        return v
    return sorted(varlist,key=cvnnums)

class G2VarObj(object):
    '''Defines a GSAS-II variable either using the phase/atom/histogram
    unique Id numbers or using a character string that specifies
    variables by phase/atom/histogram number (which can change).
    Note that :func:`GSASIIstrIO.GetUsedHistogramsAndPhases`, 
    which calls :func:`IndexAllIds` (or 
    :func:`GSASIIscriptable.G2Project.index_ids`) should be used to 
    (re)load the current Ids
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
                    raise Exception("Incorrect number of colons in var name "+str(args[0]))
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
            raise Exception("Incorrectly called GSAS-II parameter name")

        #print "DEBUG: created ",self.phase,self.histogram,self.name,self.atom

    def __str__(self):
        return self.varname()

    def __hash__(self):
        'Allow G2VarObj to be a dict key by implementing hashing'
        return hash(self.varname())

    def varname(self,hist=None):
        '''Formats the GSAS-II variable name as a "traditional" GSAS-II variable
        string (p:h:<var>:a) or (p:h:<var>)

        :param str/int hist: if specified, overrides the histogram number
          with the specified value
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
        if hist is not None and self.histogram:
            hist = str(hist)
        elif self.histogram == "*":
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
            elif ":" in str(self.atom):
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
        '''Allow comparison of G2VarObj to other G2VarObj objects or strings.
        If any field is a wildcard ('*') that field matches.
        '''
        if type(other) is str:
            other = G2VarObj(other)
        elif type(other) is not G2VarObj:
            raise Exception("Invalid type ({}) for G2VarObj comparison with {}"
                            .format(type(other),other))
        if self.phase != other.phase and self.phase != '*' and other.phase != '*':
            return False
        if self.histogram != other.histogram and self.histogram != '*' and other.histogram != '*':
            return False
        if self.atom != other.atom and self.atom != '*' and other.atom != '*':
            return False
        if self.name != other.name:
            return False
        return True
    
    def fmtVarByMode(self, seqmode, note, warnmsg):
        '''Format a parameter object for display. Note that these changes 
        are only temporary and are only shown only when the Constraints 
        data tree is selected.

        * In a non-sequential refinement or where the mode is 'use-all', the 
          name is converted unchanged to a str
        * In a sequential refinement when the mode is 'wildcards-only' the 
          name is converted unchanged to a str but a warning is added 
          for non-wildcarded HAP or Histogram parameters
        * In a sequential refinement or where the mode is 'auto-wildcard', 
          a histogram number is converted to a wildcard (*) and then 
          converted to str

        :param str mode: the sequential mode (see above)
        :param str note: value displayed on the line of the constraint/equiv.
        :param str warnmsg: a message saying the constraint is not used

        :returns: varname, explain, note, warnmsg (all str values) where:

          * varname is the parameter expressed as a string,
          * explain is blank unless there is a warning explanation about 
            the parameter or blank
          * note is the previous value unless overridden 
          * warnmsg is the previous value unless overridden 
        '''
        explain = ''
        s = self.varname()
        if seqmode == 'auto-wildcard':
            if self.histogram: s = self.varname('*')
        elif seqmode == 'wildcards-only' and self.histogram:
            if self.histogram != '*':
                warnmsg = 'Ignored due to use of a non-wildcarded histogram number'
                note = 'Ignored'
                explain = '\nIgnoring: '+self.varname()+' does not contain a wildcard.\n'
        elif seqmode != 'use-all' and seqmode != 'wildcards-only':
            print('Unexpected mode',seqmode,' in fmtVarByMode')
        return s,explain,note,warnmsg

    def _show(self):
        'For testing, shows the current lookup table'
        print ('phases'+ self.IDdict['phases'])
        print ('hists'+ self.IDdict['hists'])
        print ('atomDict'+ self.IDdict['atoms'])

#==========================================================================
def SetDefaultSample():
    'Fills in default items for the Sample dictionary for Debye-Scherrer & SASD'
    return {
        'InstrName':'',
        'ranId':ran.randint(0,sys.maxsize),
        'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],
        'DisplaceX':[0.0,False],'DisplaceY':[0.0,False],
        'Temperature':300.,'Pressure':0.1,'Time':0.0,
        'FreePrm1':0.,'FreePrm2':0.,'FreePrm3':0.,
        'Gonio. radius':200.0,
        'Omega':0.0,'Chi':0.0,'Phi':0.0,'Azimuth':0.0,
#SASD items
        'Materials':[{'Name':'vacuum','VolFrac':1.0,},{'Name':'vacuum','VolFrac':0.0,}],
        'Thick':1.0,'Contrast':[0.0,0.0],       #contrast & anomalous contrast
        'Trans':1.0,                            #measured transmission
        'SlitLen':0.0,                          #Slit length - in Q(A-1)
        }
######################################################################
class ImportBaseclass(object):
    '''Defines a base class for the reading of input files (diffraction
    data, coordinates,...). See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use a subclass of this class.
    '''
    class ImportException(Exception):
        '''Defines an Exception that is used when an import routine hits an expected error,
        usually in .Reader.

        Good practice is that the Reader should define a value in self.errors that
        tells the user some information about what is wrong with their file.
        '''
        pass

    UseReader = True  # in __init__ set value of self.UseReader to False to skip use of current importer
    def __init__(self,formatName,longFormatName=None,
                 extensionlist=[],strictExtension=False,):
        self.formatName = formatName # short string naming file type
        if longFormatName: # longer string naming file type
            self.longFormatName = longFormatName
        else:
            self.longFormatName = formatName
        # define extensions that are allowed for the file type
        # for windows, remove any extensions that are duplicate, as case is ignored
        if sys.platform == 'windows' and extensionlist:
            extensionlist = list(set([s.lower() for s in extensionlist]))
        self.extensionlist = extensionlist
        # If strictExtension is True, the file will not be read, unless
        # the extension matches one in the extensionlist
        self.strictExtension = strictExtension
        self.errors = ''
        self.warnings = ''
        self.SciPy = False          #image reader needed scipy
        # used for readers that will use multiple passes to read
        # more than one data block
        self.repeat = False
        self.selections = []
        self.repeatcount = 0
        self.readfilename = '?'
        self.scriptable = False
        #print 'created',self.__class__

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        self.errors = ''
        self.warnings = ''
        self.SciPy = False          #image reader needed scipy
        self.repeat = False
        self.repeatcount = 0
        self.readfilename = '?'


#    def Reader(self, filename, filepointer, ParentFrame=None, **unused):
#        '''This method must be supplied in the child class to read the file.
#        if the read fails either return False or raise an Exception
#        preferably of type ImportException.
#        '''
#        #start reading
#        raise ImportException("Error occurred while...")
#        self.errors += "Hint for user on why the error occur
#        return False # if an error occurs
#        return True # if read OK

    def ExtensionValidator(self, filename):
        '''This methods checks if the file has the correct extension
        
        :returns:
        
          * False if this filename will not be supported by this reader (only
            when strictExtension is True)
          * True if the extension matches the list supplied by the reader
          * None if the reader allows un-registered extensions
          
        '''
        if filename:
            ext = os.path.splitext(filename)[1]
            if not ext and self.strictExtension: return False
            for ext in self.extensionlist:                
                if sys.platform == 'windows':
                    if filename.lower().endswith(ext): return True
                else:
                    if filename.endswith(ext): return True
        if self.strictExtension:
            return False
        else:
            return None

    def ContentsValidator(self, filename):
        '''This routine will attempt to determine if the file can be read
        with the current format.
        This will typically be overridden with a method that
        takes a quick scan of [some of]
        the file contents to do a "sanity" check if the file
        appears to match the selected format.
        the file must be opened here with the correct format (binary/text)
        '''
        #filepointer.seek(0) # rewind the file pointer
        return True

    def CIFValidator(self, filepointer):
        '''A :meth:`ContentsValidator` for use to validate CIF files.
        '''
        filepointer.seek(0)
        for i,l in enumerate(filepointer):
            if i >= 1000: return True
            '''Encountered only blank lines or comments in first 1000
            lines. This is unlikely, but assume it is CIF anyway, since we are
            even less likely to find a file with nothing but hashes and
            blank lines'''
            line = l.strip()
            if len(line) == 0: # ignore blank lines
                continue
            elif line.startswith('#'): # ignore comments
                continue
            elif line.startswith('data_'): # on the right track, accept this file
                return True
            else: # found something invalid
                self.errors = 'line '+str(i+1)+' contains unexpected data:\n'
                if all([ord(c) < 128 and ord(c) != 0 for c in str(l)]): # show only if ASCII
                    self.errors += '  '+str(l)
                else:
                    self.errors += '  (binary)'
                self.errors += '\n  Note: a CIF should only have blank lines or comments before'
                self.errors += '\n        a data_ statement begins a block.'
                return False

######################################################################
class ImportPhase(ImportBaseclass):
    '''Defines a base class for the reading of files with coordinates

    Objects constructed that subclass this (in import/G2phase_*.py etc.) will be used
    in :meth:`GSASIIdataGUI.GSASII.OnImportPhase` and in 
    :func:`GSASIIscriptable.import_generic`.
    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use this class.

    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):
        # call parent __init__
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)
        self.Phase = None # a phase must be created with G2IO.SetNewPhase in the Reader
        self.SymOps = {} # specified when symmetry ops are in file (e.g. CIF)
        self.Constraints = None

######################################################################
class ImportStructFactor(ImportBaseclass):
    '''Defines a base class for the reading of files with tables
    of structure factors.

    Structure factors are read with a call to :meth:`GSASIIdataGUI.GSASII.OnImportSfact`
    which in turn calls :meth:`GSASIIdataGUI.GSASII.OnImportGeneric`, which calls
    methods :meth:`ExtensionValidator`, :meth:`ContentsValidator` and
    :meth:`Reader`.

    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use import classes in general. The specifics
    for reading a structure factor histogram require that
    the ``Reader()`` routine in the import
    class need to do only a few things: It
    should load :attr:`RefDict` item ``'RefList'`` with the reflection list,
    and set :attr:`Parameters` with the instrument parameters
    (initialized with :meth:`InitParameters` and set with :meth:`UpdateParameters`).
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)

        # define contents of Structure Factor entry
        self.Parameters = []
        'self.Parameters is a list with two dicts for data parameter settings'
        self.InitParameters()
        self.RefDict = {'RefList':[],'FF':{},'Super':0}
        self.Banks = []             #for multi bank data (usually TOF)
        '''self.RefDict is a dict containing the reflection information, as read from the file.
        Item 'RefList' contains the reflection information. See the
        :ref:`Single Crystal Reflection Data Structure<XtalRefl_table>`
        for the contents of each row. Dict element 'FF'
        contains the form factor values for each element type; if this entry
        is left as initialized (an empty list) it will be initialized as needed later.
        '''
    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        ImportBaseclass.ReInitialize(self)
        self.InitParameters()
        self.Banks = []             #for multi bank data (usually TOF)
        self.RefDict = {'RefList':[],'FF':{},'Super':0}

    def InitParameters(self):
        'initialize the instrument parameters structure'
        Lambda = 0.70926
        HistType = 'SXC'
        self.Parameters = [{'Type':[HistType,HistType], # create the structure
                            'Lam':[Lambda,Lambda]
                            }, {}]
        'Parameters is a list with two dicts for data parameter settings'

    def UpdateParameters(self,Type=None,Wave=None):
        'Revise the instrument parameters'
        if Type is not None:
            self.Parameters[0]['Type'] = [Type,Type]
        if Wave is not None:
            self.Parameters[0]['Lam'] = [Wave,Wave]

######################################################################
class ImportPowderData(ImportBaseclass):
    '''Defines a base class for the reading of files with powder data.

    Objects constructed that subclass this (in import/G2pwd_*.py etc.) will be used
    in :meth:`GSASIIdataGUI.GSASII.OnImportPowder` and in 
    :func:`GSASIIscriptable.import_generic`.
    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use this class.
    '''
    def __init__(self,formatName,longFormatName=None,
        extensionlist=[],strictExtension=False,):
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)
        self.clockWd = None  # used in TOF
        self.ReInitialize()

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        ImportBaseclass.ReInitialize(self)
        self.powderentry = ['',None,None] #  (filename,Pos,Bank)
        self.powderdata = [] # Powder dataset
        '''A powder data set is a list with items [x,y,w,yc,yb,yd]:
                np.array(x), # x-axis values
                np.array(y), # powder pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.array(yc), # calc. intensities (zero)
                np.array(yb), # calc. background (zero)
                np.array(yd), # obs-calc profiles
        '''
        self.comments = []
        self.idstring = ''
        self.Sample = SetDefaultSample() # default sample parameters
        self.Controls = {}  # items to be placed in top-level Controls
        self.GSAS = None     # used in TOF
        self.repeat_instparm = True # Should a parm file be
        #                             used for multiple histograms?
        self.instparm = None # name hint from file of instparm to use
        self.instfile = '' # full path name to instrument parameter file
        self.instbank = '' # inst parm bank number
        self.instmsg = ''  # a label that gets printed to show
                           # where instrument parameters are from
        self.numbanks = 1
        self.instdict = {} # place items here that will be transferred to the instrument parameters
        self.pwdparms = {} # place parameters that are transferred directly to the tree
                           # here (typically from an existing GPX file)
######################################################################
class ImportSmallAngleData(ImportBaseclass):
    '''Defines a base class for the reading of files with small angle data.
    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use this class.
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):

        ImportBaseclass.__init__(self,formatName,longFormatName,extensionlist,
            strictExtension)
        self.ReInitialize()

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        ImportBaseclass.ReInitialize(self)
        self.smallangleentry = ['',None,None] #  (filename,Pos,Bank)
        self.smallangledata = [] # SASD dataset
        '''A small angle data set is a list with items [x,y,w,yc,yd]:
                np.array(x), # x-axis values
                np.array(y), # powder pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.array(yc), # calc. intensities (zero)
                np.array(yd), # obs-calc profiles
                np.array(yb), # preset bkg
        '''
        self.comments = []
        self.idstring = ''
        self.Sample = SetDefaultSample()
        self.GSAS = None     # used in TOF
        self.clockWd = None  # used in TOF
        self.numbanks = 1
        self.instdict = {} # place items here that will be transferred to the instrument parameters

######################################################################
class ImportReflectometryData(ImportBaseclass):
    '''Defines a base class for the reading of files with reflectometry data.
    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use this class.
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):

        ImportBaseclass.__init__(self,formatName,longFormatName,extensionlist,
            strictExtension)
        self.ReInitialize()

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        ImportBaseclass.ReInitialize(self)
        self.reflectometryentry = ['',None,None] #  (filename,Pos,Bank)
        self.reflectometrydata = [] # SASD dataset
        '''A small angle data set is a list with items [x,y,w,yc,yd]:
                np.array(x), # x-axis values
                np.array(y), # powder pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.array(yc), # calc. intensities (zero)
                np.array(yd), # obs-calc profiles
                np.array(yb), # preset bkg
        '''
        self.comments = []
        self.idstring = ''
        self.Sample = SetDefaultSample()
        self.GSAS = None     # used in TOF
        self.clockWd = None  # used in TOF
        self.numbanks = 1
        self.instdict = {} # place items here that will be transferred to the instrument parameters

######################################################################
class ImportPDFData(ImportBaseclass):
    '''Defines a base class for the reading of files with PDF G(R) data.
    See :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use this class.
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):

        ImportBaseclass.__init__(self,formatName,longFormatName,extensionlist,
            strictExtension)
        self.ReInitialize()

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings'
        ImportBaseclass.ReInitialize(self)
        self.pdfentry = ['',None,None] #  (filename,Pos,Bank)
        self.pdfdata = [] # PDF G(R) dataset
        '''A pdf g(r) data set is a list with items [x,y]:
                np.array(x), # r-axis values
                np.array(y), # pdf g(r)
        '''
        self.comments = []
        self.idstring = ''
        self.numbanks = 1

######################################################################
class ImportImage(ImportBaseclass):
    '''Defines a base class for the reading of images

    Images are read in only these places:

      * Initial reading is typically done from a menu item
        with a call to :meth:`GSASIIdataGUI.GSASII.OnImportImage`
        which in turn calls :meth:`GSASIIdataGUI.GSASII.OnImportGeneric`. That calls
        methods :meth:`ExtensionValidator`, :meth:`ContentsValidator` and
        :meth:`Reader`. This returns a list of reader objects for each read image.
        Also used in :func:`GSASIIscriptable.import_generic`.

      * Images are read alternatively in :func:`GSASIIIO.ReadImages`, which puts image info
        directly into the data tree.

      * Unlike all other data types read by GSAS-II, images are only kept in memory as 
        they are used and function :func:`GSASIIIO.GetImageData` is used to reread images
        if they are reloaded. For quick retrieval of previously read images, it may be useful to 
        save sums of images or save a keyword (see ``ImageTag``, below

    When reading an image, the ``Reader()`` routine in the ImportImage class
    should set:

      * :attr:`Comments`: a list of strings (str),
      * :attr:`Npix`: the number of pixels in the image (int),
      * :attr:`Image`: the actual image as a numpy array (np.array)
      * :attr:`Data`: a dict defining image parameters (dict). Within this dict the following
        data items are used:

         * ``pixelSize``: size of each pixel (x,y) in microns (such as ``[200.,200.]``.
         * ``wavelength``: wavelength in :math:`\\AA`.
         * ``distance``: distance of detector from sample in cm.
         * ``center``: uncalibrated center of beam on detector (such as ``[204.8,204.8]``, in mm
           measured from top left corner of the detector
         * ``size``: size of image in pixels (x,y) (such as ``[2048,2048]``).
         * ``ImageTag``: image number or other keyword used to retrieve image from
           a multi-image data file (defaults to ``1`` if not specified).
         * ``sumfile``: holds sum image file name if a sum was produced from a multi image file
         * ``PolaVal``: has two values, the polarization fraction (typically 0.95-0.99 
           for synchrotrons, 0.5 for lab instruments) and a refinement flag 
           (such as ``[0.99, False]``).
         * ``setdist``: nominal distance from sample to detector. Note that ``distance`` may 
           be changed during calibration, but ``setdist`` will not be, so that calibration may be 
           repeated. 

    optional data items:

      * :attr:`repeat`: set to True if there are additional images to
        read in the file, False otherwise
      * :attr:`repeatcount`: set to the number of the image.

    Note that the above is initialized with :meth:`InitParameters`.
    (Also see :ref:`Writing a Import Routine<import_routines>`
    for an explanation on how to use import classes in general.)
    '''
    def __init__(self,formatName,longFormatName=None,extensionlist=[],
        strictExtension=False,):
        ImportBaseclass.__init__(self,formatName,longFormatName,
            extensionlist,strictExtension)
        self.InitParameters()

    def ReInitialize(self):
        'Reinitialize the Reader to initial settings -- not used at present'
        ImportBaseclass.ReInitialize(self)
        self.InitParameters()

    def InitParameters(self):
        'initialize the instrument parameters structure'
        self.Comments = ['No comments']
        self.Data = {'samplechangerpos':0.0,'det2theta':0.0,'Gain map':''}
        self.Npix = 0
        self.Image = None
        self.repeat = False
        self.repeatcount = 1
        self.sumfile = ''

    def LoadImage(self,ParentFrame,imagefile,imagetag=None):
        '''Optionally, call this after reading in an image to load it into the tree.
        This saves time by preventing a reread of the same information.
        '''
        if ParentFrame:
            ParentFrame.ImageZ = self.Image   # store the image for plotting
            ParentFrame.oldImagefile = imagefile # save the name of the last image file read
            ParentFrame.oldImageTag = imagetag   # save the tag of the last image file read

#################################################################################################
# shortcut routines
exp = np.exp
sind = sin = s = lambda x: np.sin(x*np.pi/180.)
cosd = cos = c = lambda x: np.cos(x*np.pi/180.)
tand = tan = t = lambda x: np.tan(x*np.pi/180.)
sqrt = sq = lambda x: np.sqrt(x)
pi = lambda: np.pi

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
    # no listed module name, try in current namespace
    if len(df) == 1:
        try:
            fxnobj = eval(f)
            return pkgdict,fxnobj
        except (AttributeError, NameError):
            return None,None

    # includes a package, see if package is already imported
    pkgnam = '.'.join(df[:-1])
    try:
        fxnobj = eval(f)
        pkgdict[pkgnam] = eval(pkgnam)
        return pkgdict,fxnobj
    except (AttributeError, NameError):
        pass
    # package not yet imported, so let's try
    if '.' not in sys.path: sys.path.append('.')
    pkgnam = '.'.join(df[:-1])
    #for pkg in f.split('.')[:-1]: # if needed, descend down the tree
    #    if pkgname:
    #        pkgname += '.' + pkg
    #    else:
    #        pkgname = pkg
    try:
        exec('import '+pkgnam)
        pkgdict[pkgnam] = eval(pkgnam)
        fxnobj = eval(f)
    except Exception as msg:
        print('load of '+pkgnam+' failed with error='+str(msg))
        return {},None
    # can we access the function? I am not exactly sure what
    #    I intended this to test originally (BHT)
    try:
        fxnobj = eval(f,globals(),pkgdict)
        return pkgdict,fxnobj
    except Exception as msg:
        print('call to',f,' failed with error=',str(msg))
        return None,None # not found
                
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

        self.distance_dict  = None  # to be used for defining atom phase/symmetry info
        self.distance_atoms = None  # to be used for defining atom distances

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
            indx = list(varyList).index("::"+self.freeVars[v][0])
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
        defined = list(self.assgnVars.keys()) + list(self.freeVars.keys())
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
        except SyntaxError:
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
        self.su = None
        '''Standard error evaluation where supplied by the evaluator
        '''
        # Patch: for old-style expressions with a (now removed step size)
        for v in self.eObj.assgnVars:
            if not isinstance(self.eObj.assgnVars[v], str):
                self.eObj.assgnVars[v] = self.eObj.assgnVars[v][0]
        self.parmDict = {}
        '''A copy of the parameter dictionary, for distance and angle computation
        '''

    def SetupCalc(self,parmDict):
        '''Do all preparations to use the expression for computation.
        Adds the free parameter values to the parameter dict (parmDict).
        '''
        if self.eObj.expression.startswith('Dist') or self.eObj.expression.startswith('Angle'):
            return
        self.fxnpkgdict = self.eObj.CheckVars()
        # all is OK, compile the expression
        self.compiledExpr = compile(self.eObj.expression,'','eval')

        # look at first value in parmDict to determine its type
        parmsInList = True
        for key in parmDict:
            val = parmDict[key]
            if isinstance(val, str):
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
            if varname in parmDict:
                self.lblLookup[varname] = v
                self.varLookup[v] = varname
                if parmsInList:
                    self.exprDict[v] = parmDict[varname][0]
                else:
                    self.exprDict[v] = parmDict[varname]
            elif '*' in varname:
                varlist = LookupWildCard(varname,list(parmDict.keys()))
                if len(varlist) == 0:
                    self.exprDict[v] = None
                    self.lblLookup[v] = varname # needed?
                    self.exprDict.update(self.fxnpkgdict) # needed?
                    return
                for var in varlist:
                    self.lblLookup[var] = v
                if parmsInList:
                    self.exprDict[v] = np.array([parmDict[var][0] for var in varlist])
                else:
                    self.exprDict[v] = np.array([parmDict[var] for var in varlist])
                self.varLookup[v] = [var for var in varlist]
            else:
                self.exprDict[v] = None
#                raise Exception,"No value for variable "+str(v)
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
        :param dict parmDict: a dict of values, items not in use are ignored
        '''
        if self.eObj.expression.startswith('Dist') or self.eObj.expression.startswith('Angle'):
            self.parmDict = parmDict
            return
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
        self.su = None
        if self.eObj.expression.startswith('Dist'):
#            GSASIIpath.IPyBreak()
            dist = G2mth.CalcDist(self.eObj.distance_dict, self.eObj.distance_atoms, self.parmDict)
            return dist
        elif self.eObj.expression.startswith('Angle'):
            angle = G2mth.CalcAngle(self.eObj.angle_dict, self.eObj.angle_atoms, self.parmDict)
            return angle
        if self.compiledExpr is None:
            raise Exception("EvalExpression called before SetupCalc")
        try:
            val = eval(self.compiledExpr,globals(),self.exprDict)
        except TypeError:
            val = None
        except NameError:
            val = None
        if not np.isscalar(val):
            val = np.sum(val)
        return val

def makeAngleObj(Phase,Oatom,Tatoms):
    General = Phase['General']
    cx,ct = General['AtomPtrs'][:2]
    pId = Phase['pId']
    SGData = General['SGData']
    Atoms = Phase['Atoms']
    aNames = [atom[ct-1] for atom in Atoms]
    tIds = []
    symNos = []
    cellNos = []
    oId = aNames.index(Oatom)
    for Tatom in Tatoms.split(';'):
        sB = Tatom.find('(')+1
        symNo = 0
        if sB:
            sF = Tatom.find(')')
            symNo = int(Tatom[sB:sF])
        symNos.append(symNo)
        cellNo = [0,0,0]
        cB = Tatom.find('[')
        if cB>0:
            cF = Tatom.find(']')+1
            cellNo = eval(Tatom[cB:cF])
        cellNos.append(cellNo)
        tIds.append(aNames.index(Tatom.split('+')[0]))
    # create an expression object
    obj = ExpressionObj()
    obj.expression = 'Angle(%s,%s,\n%s)'%(Tatoms[0],Oatom,Tatoms[1])
    obj.angle_dict = {'pId':pId,'SGData':SGData,'symNo':symNos,'cellNo':cellNos}
    obj.angle_atoms = [oId,tIds]
    return obj

class G2Exception(Exception):
    'A generic GSAS-II exception class'
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class G2RefineCancel(Exception):
    'Raised when Cancel is pressed in a refinement dialog'
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)
    
def HowDidIgetHere(wherecalledonly=False):
    '''Show a traceback with calls that brought us to the current location.
    Used for debugging.

    :param bool wherecalledonly: When True, the entire calling stack is 
      shown. When False (default), only the 2nd to last stack entry (the
      routine that called the calling routine is shown.
    '''
    import traceback
    if wherecalledonly:
        i = traceback.format_list(traceback.extract_stack()[:-1])[-2]
        print(i.strip().rstrip())
    else:
        print (70*'*')
        for i in traceback.format_list(traceback.extract_stack()[:-1]): print(i.strip().rstrip())
        print (70*'*')

# Note that this is GUI code and should be moved at somepoint
def CreatePDFitems(G2frame,PWDRtree,ElList,Qlimits,numAtm=1,FltBkg=0,PDFnames=[]):
    '''Create and initialize a new set of PDF tree entries

    :param Frame G2frame: main GSAS-II tree frame object
    :param str PWDRtree: name of PWDR to be used to create PDF item
    :param dict ElList: data structure with composition
    :param list Qlimits: Q limits to be used for computing the PDF
    :param float numAtm: no. atom in chemical formula
    :param float FltBkg: flat background value
    :param list PDFnames: previously used PDF names

    :returns: the Id of the newly created PDF entry
    '''
    PDFname = 'PDF '+PWDRtree[4:] # this places two spaces after PDF, which is needed is some places
    if PDFname in PDFnames:
        print('Skipping, entry already exists: '+PDFname)
        return None
    #PDFname = MakeUniqueLabel(PDFname,PDFnames)
    Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=PDFname)
    Data = {
        'Sample':{'Name':PWDRtree,'Mult':1.0},
        'Sample Bkg.':{'Name':'','Mult':-1.0,'Refine':False},
        'Container':{'Name':'','Mult':-1.0,'Refine':False},
        'Container Bkg.':{'Name':'','Mult':-1.0},'ElList':ElList,
        'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':10.0*numAtm,'Flat Bkg':FltBkg,
        'DetType':'Area detector','ObliqCoeff':0.3,'Ruland':0.025,'QScaleLim':Qlimits,
        'Lorch':False,'BackRatio':0.0,'Rmax':100.,'noRing':False,'IofQmin':1.0,'Rmin':1.0,
        'I(Q)':[],'S(Q)':[],'F(Q)':[],'G(R)':[],
        #items for sequential PDFfit
        'Datarange':[0.,30.],'Fitrange':[0.,30.],'qdamp':[0.03,False],'qbroad':[0,False],'Temp':300}
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='PDF Controls'),Data)
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='PDF Peaks'),
        {'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]})
    return Id

class ShowTiming(object):
    '''An object to use for timing repeated sections of code.

    Create the object with::
       tim0 = ShowTiming()

    Tag sections of code to be timed with::
       tim0.start('start')
       tim0.start('in section 1')
       tim0.start('in section 2')
       
    etc. (Note that each section should have a unique label.)

    After the last section, end timing with::
       tim0.end()

    Show timing results with::
       tim0.show()
       
    '''
    def __init__(self):
        self.timeSum =  []
        self.timeStart = []
        self.label = []
        self.prev = None
    def start(self,label):
        import time
        if label in self.label:
            i = self.label.index(label)
            self.timeStart[i] = time.time()
        else:
            i = len(self.label)
            self.timeSum.append(0.0)
            self.timeStart.append(time.time())
            self.label.append(label)
        if self.prev is not None:
            self.timeSum[self.prev] += self.timeStart[i] - self.timeStart[self.prev]
        self.prev = i
    def end(self):
        import time
        if self.prev is not None:
            self.timeSum[self.prev] += time.time() - self.timeStart[self.prev]
        self.prev = None
    def show(self):
        sumT = sum(self.timeSum)
        print('Timing results (total={:.2f} sec)'.format(sumT))
        for i,(lbl,val) in enumerate(zip(self.label,self.timeSum)):
            print('{} {:20} {:8.2f} ms {:5.2f}%'.format(i,lbl,1000.*val,100*val/sumT))

def validateAtomDrawType(typ,generalData={}):
    '''Confirm that the selected Atom drawing type is valid for the current 
    phase. If not, use 'vdW balls'. This is currently used only for setting a 
    default when atoms are added to the atoms draw list.
    '''
    if typ in ('lines','vdW balls','sticks','balls & sticks','ellipsoids'):
        return typ
    # elif generalData.get('Type','') == 'macromolecular':
    #     if typ in ('backbone',):
    #         return typ
    return 'vdW balls'

if __name__ == "__main__":
    # test variable descriptions
    for var in '0::Afrac:*',':1:Scale','1::dAx:0','::undefined':
        v = var.split(':')[2]
        print(var+':\t', getDescr(v),getVarStep(v))
    import sys; sys.exit()
    # test equation evaluation
    def showEQ(calcobj):
        print (50*'=')
        print (calcobj.eObj.expression+'='+calcobj.EvalExpression())
        for v in sorted(calcobj.varLookup):
            print ("  "+v+'='+calcobj.exprDict[v]+'='+calcobj.varLookup[v])
        # print '  Derivatives'
        # for v in calcobj.derivStep.keys():
        #     print '    d(Expr)/d('+v+') =',calcobj.EvalDeriv(v)

    obj = ExpressionObj()

    obj.expression = "A*np.exp(B)"
    obj.assgnVars =  {'B': '0::Afrac:1'}
    obj.freeVars =  {'A': [u'A', 0.5, True]}
    #obj.CheckVars()
    calcobj = ExpressionCalcObj(obj)

    obj1 = ExpressionObj()
    obj1.expression = "A*np.exp(B)"
    obj1.assgnVars =  {'B': '0::Afrac:*'}
    obj1.freeVars =  {'A': [u'Free Prm A', 0.5, True]}
    #obj.CheckVars()
    calcobj1 = ExpressionCalcObj(obj1)

    obj2 = ExpressionObj()
    obj2.distance_stuff = np.array([[0,1],[1,-1]])
    obj2.expression = "Dist(1,2)"
    GSASIIpath.InvokeDebugOpts()
    parmDict2 = {'0::Afrac:0':[0.0,True], '0::Afrac:1': [1.0,False]}
    calcobj2 = ExpressionCalcObj(obj2)
    calcobj2.SetupCalc(parmDict2)
    showEQ(calcobj2)

    parmDict1 = {'0::Afrac:0':1.0, '0::Afrac:1': 1.0}
    print ('\nDict = '+parmDict1)
    calcobj.SetupCalc(parmDict1)
    showEQ(calcobj)
    calcobj1.SetupCalc(parmDict1)
    showEQ(calcobj1)

    parmDict2 = {'0::Afrac:0':[0.0,True], '0::Afrac:1': [1.0,False]}
    print ('Dict = '+parmDict2)
    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)
    calcobj1.SetupCalc(parmDict2)
    showEQ(calcobj1)
    calcobj2.SetupCalc(parmDict2)
    showEQ(calcobj2)
