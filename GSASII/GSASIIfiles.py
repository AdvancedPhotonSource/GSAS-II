# -*- coding: utf-8 -*-
'''
Code for accessing files, including support for reading and writing
instrument parameter files and exporting various types of data files.

This module has some routines that require wxPython, but imports
for wx and GSAS-II GUI routines is done on a per-function basis so
that this module can be imported for GSASIIscriptable use when
wx is not installed.
'''
from __future__ import division, print_function
import platform
import os
import sys
import glob
#import inspect
import re

import numpy as np

from . import GSASIIpath
from . import GSASIIlattice as G2lat
from . import GSASIIstrIO as G2stIO
from . import GSASIImapvars as G2mv
from . import GSASIImath as G2mth
#from . import GSASIIlattice as G2lat

#if not sys.platform.startswith('win'):
#    try:
#        from dmp import dump2tmp,undumptmp
#    except:
#        print('Note: Import of dmp skipped')

# declare symbol (pi) and functions allowed in expressions
sind = sin = s = lambda x: np.sin(x*np.pi/180.)
cosd = cos = c = lambda x: np.cos(x*np.pi/180.)
tand = tan = t = lambda x: np.tan(x*np.pi/180.)
sqrt = sq = lambda x: np.sqrt(x)
pi = np.pi

# N.B. This is duplicated in G2IO
def sfloat(S):
    'Convert a string to float. An empty field or a unconvertable value is treated as zero'
    if S.strip():
        try:
            return float(S)
        except ValueError:
            pass
    return 0.0

G2printLevel = 'all'
'''This defines the level of output from calls to :func:`GSASIIfiles.G2Print`,
which should  be used in place of print() within GSASII where possible.
Settings for this are 'all', 'warn', 'error' or 'none'. Best to change this
with :func:`G2SetPrintLevel`.

.. seealso::
    :func:`G2Print`
    :func:`G2SetPrintLevel`.
'''

def G2SetPrintLevel(level):
    '''Set the level of output from calls to :func:`G2Print`, which should
    be used in place of print() within GSASII. Settings for the mode are
    'all', 'warn', 'error' or 'none'

    :param str level: a string used to set the print level, which may be
      'all', 'warn', 'error' or 'none'.
      Note that capitalization and extra letters in level are ignored, so
      'Warn', 'warnings', etc. will all set the mode to 'warn'
    '''
    global G2printLevel
    for mode in  'all', 'warn', 'error', 'none':
        if mode in level.lower():
            G2printLevel = mode
            return
    else:
        G2Print('G2SetPrintLevel Error: level={} cannot be interpreted.',
                    'Use all, warn, error or none.')

def find(name, path):
    '''find 1st occurance of file in path
    '''
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def G2Print(*args,**kwargs):
    '''Print with filtering based level of output (see :func:`G2SetPrintLevel`).
    Use G2Print() as replacement for print().

    :param str mode: if specified, this should contain the mode for printing
      ('error', 'warn' or anything else). If not specified, the first argument
      of the print command (args[0]) should contain the string 'error' for
      error messages and 'warn' for warning messages
      (capitalization and additional letters ignored.)
    '''
    if G2printLevel == 'none': return
    if kwargs.get('mode') is None:
        testStr = str(args[0]).lower()
    else:
        testStr = kwargs['mode'][:].lower()
        del kwargs['mode']
    level = 2
    for i,mode in enumerate(('error', 'warn')):
        if mode in testStr:
            level = i
            break
    if G2printLevel == 'error' and level > 0: return
    if G2printLevel == 'warn' and level > 1: return
    print(*args,**kwargs)

def get_python_versions(packagelist):
    versions = [['Python', sys.version.split()[0]]]
    for pack in packagelist:
        try:
            versions.append([pack.__name__, pack.__version__])
        except:
            pass
    versions.append(['Platform',
                     sys.platform + ' ' + platform.architecture()[0] +
                     ' ' + platform.machine()])
    return versions

ImportErrors = []
condaRequestList = {}
def ImportErrorMsg(errormsg=None,pkg={}):
    '''Store error message(s) from loading importers (usually missing
    packages. Or, report back all messages, if called with no argument.

    :param str errormsg: a string containing the error message. If not
      supplied, the function returns the error message(s).
    :param dict pkg: a dict where the key is the name of the importer and the
      value is a list containing the packages that need to be installed to
      allow the importer to be used.

    :returns: the error messages as a list (an empty list if there are none),
      only if errormsg is None (the default).
    '''
    if errormsg is None:
        return ImportErrors
    ImportErrors.append(errormsg)
    if pkg: NeededPackage(pkg)

def NeededPackage(pkgDict):
    '''Store packages that are needed to add functionality to GSAS-II

    :param dict pkgDict: a dict where the key is the describes the routine
      or the unavailable functionality that requires addition of a Python
      package and the value is a list containing the packages that need to 
      be installed. If a specific version of a package is required then 
      indicate that by placing version information into the name (see 
      https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf). Note that the names containing "pkg=x.y"
      will be translated to "pkg==x.y" for pip installs, but more complex 
      package specifications (see 
      https://pip.pypa.io/en/latest/reference/requirement-specifiers/) will 
      probably work only for conda installations.

      Examples::

          {'MIDAS Zarr importer':['zarr=2.18.*']}
          {'HDF5 image importer':['h5py','hdf5']}
    '''
    condaRequestList.update(pkgDict)

def makeInstDict(names,data,codes):
    inst = dict(zip(names,zip(data,data,codes)))
    for item in inst:
        inst[item] = list(inst[item])
    return inst

def SetPowderInstParms(Iparm, rd):
    '''extracts values from instrument parameters in rd.instdict
    or in array Iparm.
    Create and return the contents of the instrument parameter tree entry.
    '''
    Irads = {0:' ',1:'CrKa',2:'FeKa',3:'CuKa',4:'MoKa',5:'AgKa',6:'TiKa',7:'CoKa'}
    DataType = Iparm['INS   HTYPE '].strip()[:3]  # take 1st 3 chars
    # override inst values with values read from data file
    Bank = rd.powderentry[2]    #should be used in multibank iparm files
    if rd.instdict.get('type'):
        DataType = rd.instdict.get('type')
    data = [DataType,]
    instname = Iparm.get('INS  1INAME ')
    irad = int(Iparm.get('INS  1 IRAD ','0'))
    if instname:
        rd.Sample['InstrName'] = instname.strip()
    if 'C' in DataType:
        wave1 = None
        wave2 = 0.0
        if rd.instdict.get('wave'):
            wl = rd.instdict.get('wave')
            wave1 = wl[0]
            if len(wl) > 1: wave2 = wl[1]
        s = Iparm['INS  1 ICONS']
        if not wave1:
            wave1 = sfloat(s[:10])
            wave2 = sfloat(s[10:20])
        v = (wave1,wave2,
             sfloat(s[20:30])/100.,sfloat(s[55:65]),sfloat(s[40:50])) #get lam1, lam2, zero, pola & ratio
        if not v[1]:
            names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
            v = (v[0],v[2],v[4])
            codes = [0,0,0,0,0]
            rd.Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],'DisplaceY':[0.,False]})
        else:
            names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
            codes = [0,0,0,0,0,0,0]
            rd.Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
                'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
        data.extend(v)
        if 'INS  1PRCF  ' in Iparm:
            v1 = Iparm['INS  1PRCF  '].split()
            v = Iparm['INS  1PRCF 1'].split()
            data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
            azm = float(Iparm.get('INS  1DETAZM','0.0'))
            v = Iparm['INS  1PRCF 2'].split()
            if v1[0] == 3:
                data.extend([float(v[0]),float(v[1]),0.0,float(v[2])+float(v[3],azm)])  #get LX, LY, Z, S+H/L & azimuth
            else:
                data.extend([0.0,0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
        else:
            v1 = Iparm['INS  1PRCF1 '].split()
            v = Iparm['INS  1PRCF11'].split()
            data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
            azm = float(Iparm.get('INS  1DETAZM','0.0'))
            v = Iparm['INS  1PRCF12'].split()
            if v1[0] == 3:
                data.extend([float(v[0]),float(v[1]),0.0,float(v[2])+float(v[3],azm)])  #get LX, LY, Z, S+H/L & azimuth
            else:
                data.extend([0.0,0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
        codes.extend([0,0,0,0,0,0,0])
        Iparm1 = makeInstDict(names,data,codes)
        Iparm1['Source'] = [Irads[irad],Irads[irad]]
        Iparm1['Bank'] = [Bank,Bank,0]
        return [Iparm1,{}]
    elif 'T' in DataType:
        names = ['Type','fltPath','2-theta','difC','difA', 'difB','Zero','alpha','beta-0','beta-1',
            'beta-q','sig-0','sig-1','sig-2','sig-q', 'X','Y','Z','Azimuth',]
        codes = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,]
        azm = 0.
        if 'INS  1DETAZM' in Iparm:
            azm = float(Iparm['INS  1DETAZM'])
        rd.Sample['Azimuth'] = azm
        fltPath0 = 20.                      #arbitrary
        if 'INS   FPATH1' in Iparm:
            s = Iparm['INS   FPATH1'].split()
            fltPath0 = sfloat(s[0])
        if 'INS  1BNKPAR' not in Iparm:     #bank missing from Iparm file
            return []
        s = Iparm['INS  1BNKPAR'].split()
        fltPath1 = sfloat(s[0])
        data.extend([fltPath0+fltPath1,])               #Flight path source-sample-detector
        data.extend([sfloat(s[1]),])               #2-theta for bank
        s = Iparm['INS  1 ICONS'].split()
        data.extend([sfloat(s[0]),sfloat(s[1]),0.0,sfloat(s[2])])    #difC,difA,difB,Zero
        if 'INS  1PRCF  ' in Iparm:
            s = Iparm['INS  1PRCF  '].split()
            pfType = int(s[0])
            s = Iparm['INS  1PRCF 1'].split()
            if abs(pfType) == 1:
                data.extend([sfloat(s[1]),sfloat(s[2]),sfloat(s[3])]) #alpha, beta-0, beta-1
                s = Iparm['INS  1PRCF 2'].split()
                data.extend([0.0,0.0,sfloat(s[1]),sfloat(s[2]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
            elif abs(pfType) in [3,4,5]:
                data.extend([sfloat(s[0]),sfloat(s[1]),sfloat(s[2])]) #alpha, beta-0, beta-1
                if abs(pfType) == 4:
                    data.extend([0.0,0.0,sfloat(s[3]),0.0,0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
                else:
                    s = Iparm['INS  1PRCF 2'].split()
                    data.extend([0.0,0.0,sfloat(s[0]),sfloat(s[1]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
            elif abs(pfType) == 2:
                print('''***WARNING gsas profile function #2 does not give valid GSAS-II diffractometer/profile coefficients ***
                you should request GSAS-II instparm file from Instrument responsible''')
                data.extend([sfloat(s[1]),0.0,1./sfloat(s[3])]) #alpha, beta-0, beta-1
                data.extend([0.0,0.0,sfloat(s[1]),0.0,0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
        else:
            s = Iparm['INS  1PRCF1 '].split()
            pfType = int(s[0])
            s = Iparm['INS  1PRCF11'].split()
            if abs(pfType) == 1:
                data.extend([sfloat(s[1]),sfloat(s[2]),sfloat(s[3])]) #alpha, beta-0, beta-1
                s = Iparm['INS  1PRCF12'].split()
                data.extend([0.0,0.0,sfloat(s[1]),sfloat(s[2]),0.0,0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
            elif abs(pfType) in [3,4,5]:
                data.extend([sfloat(s[0]),sfloat(s[1]),sfloat(s[2])]) #alpha, beta-0, beta-1
                if abs(pfType) == 4:
                    data.extend([0.0,0.0,sfloat(s[3]),0.0,0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
                else:
                    s = Iparm['INS  1PRCF12'].split()
                    data.extend([0.0,0.0,sfloat(s[0]),sfloat(s[1]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y, Z
        Inst1 = makeInstDict(names,data,codes)
        Inst1['Bank'] = [Bank,Bank,0]
        Inst2 = {}
        if pfType < 0:
            Ipab = 'INS  1PAB'+str(-pfType)
            Npab = int(Iparm[Ipab+'  '].strip())
            Inst2['Pdabc'] = []
            for i in range(Npab):
                k = Ipab+str(i+1).rjust(2)
                s = Iparm[k].split()
                Inst2['Pdabc'].append([float(t) for t in s])
            Inst2['Pdabc'] = np.array(Inst2['Pdabc'])
            Inst2['Pdabc'].T[3] += Inst2['Pdabc'].T[0]*Inst1['difC'][0] #turn 3rd col into TOF
        if 'INS  1I ITYP' in Iparm:
            s = Iparm['INS  1I ITYP'].split()
            Ityp = int(s[0])
            Tminmax = [float(s[1])*1000.,float(s[2])*1000.]
            Itypes = ['Exponential','Maxwell/Exponential','','Maxwell/Chebyschev','']
            if Ityp in [1,2,4]:
                Inst2['Itype'] = Itypes[Ityp-1]
                Inst2['Tminmax'] = Tminmax
                Icoeff = []
                Iesd = []
                Icovar = []
                for i in range(3):
                    s = Iparm['INS  1ICOFF'+str(i+1)].split()
                    Icoeff += [float(S) for S in s]
                    s = Iparm['INS  1IECOF'+str(i+1)].split()
                    Iesd += [float(S) for S in s]
                NT = 10
                for i in range(8):
                    s = Iparm['INS  1IECOR'+str(i+1)]
                    if i == 7:
                        NT = 8
                    Icovar += [float(s[6*j:6*j+6]) for j in range(NT)]
                Inst2['Icoeff'] = Icoeff
                Inst2['Iesd'] = Iesd
                Inst2['Icovar'] = Icovar
        return [Inst1,Inst2]
    elif 'E' in DataType:
        tth = float(Iparm['INS  1 ICONS'])
        s = Iparm['INS  1PRCF11'].split()
        names = ['Type','2-theta','XE','YE','ZE','A','B', 'C']
        codes = [0,0,0,0,0,0,0,0]
        data.extend([tth,0.0,0.0,0.0,sfloat(s[0]),sfloat(s[1]),sfloat(s[2])])    #A,B,C
        Iparm1 = makeInstDict(names,data,codes)
        Iparm1['Bank'] = [Bank,Bank,0]
        return [Iparm1,{}]

def ReadInstprm(instLines, bank, Sample={}):
    '''Read contents of a GSAS-II (new) .instprm instrument parameter file

    :param list instLines: contents of GSAS-II parameter file as a
          list of str; N.B. lines can be concatenated with ';'
    :param int bank: bank number to use when instprm file has values
          for multiple banks (noted by headers of '#BANK n:...'.). This
          is ignored for instprm files without those headers.
          If bank is None with multiple banks, the first bank is used.
          Note that multibank .instprm files are made by
          a "Save all profile" command in Instrument Parameters.
    :param dict Sample: A dict containing sample parameters,
           typically corresponding to rd.Sample,
           where rd is a reader object that
           is being read from. Sample parameters
           determined by instrument settings or information
           from the instprm file are placed here.
    :returns: bank,instdict where bank is the sample parameter set
           number and instdict is the instrument parameter dict

    Note if 'Type' is set as Debye-Scherrer or Bragg-Brentano this will be used and
    will set defaults in the sample parameters. Otherwise, a single-wavelength file
    will set Debye-Scherrer mode and dual wavelength will set Bragg-Brentano.
    '''
    if 'GSAS-II' not in instLines[0]:
        raise ValueError("Not a valid GSAS-II instprm file")

    newItems = []
    newVals = []
    NewSample = {}
    Found = False
    il = 0
    if bank is None:
        banklist = set()
        for S in instLines:
            if S[0] == '#' and 'Bank' in S:
                banklist.add(int(S.split(':')[0].split()[1]))
        # Picks the first bank by default
        if len(banklist) > 1:
            bank = sorted(banklist)[0]
        else:
            bank = 1
#        rd.powderentry[2] = bank
    while il < len(instLines):
        S = instLines[il]
        if S[0] == '#':
            if Found:
                break
            if 'Bank' in S:
                if bank == int(S.split(':')[0].split()[1]):
                    il += 1
                    S = instLines[il]
                else:
                    il += 1
                    S = instLines[il]
                    while il < len(instLines) and '#Bank' not in S:
                        il += 1
                        if il == len(instLines):
                            raise ValueError("Bank {} not found in instprm file".format(bank))
                        S = instLines[il]
                    continue
            else:
                il += 1
                S = instLines[il]
        Found = True
        if '"""' in S:
            delim = '"""'
        elif "'''" in S:
            delim = "'''"
        else:
            S = S.replace(' ', '')
            SS = S.strip().split(';')
            for s in SS:
                item, val = s.split(':', 1)
                try:
                    val = float(val)
                except:
                    pass
                if item == 'Gonio.radius':
                    NewSample.update({'Gonio. radius':float(val)})
                elif item == 'Diff-type':
                    NewSample.update({'Type':val})
                elif item == 'InstrName':
                    NewSample.update({item:val})
                else:
                    newItems.append(item)
                    newVals.append(val)
            il += 1
            continue
        # read multiline values, delimited by ''' or """
        item, val = S.strip().split(':', 1)
        if item in ['XE','YE','ZE','WE']:   #skip cn to keV conversion factors
            il += 1
            continue
        val = val.replace(delim, '').rstrip()
        val += '\n'
        while True:
            il += 1
            if il >= len(instLines):
                break
            S = instLines[il]
            if delim in S:
                val += S.replace(delim, '').rstrip()
                val += '\n'
                break
            else:
                val += S.rstrip()
                val += '\n'
        newItems.append(item)
        newVals.append(val)
        il += 1
    if NewSample.get('Type','') == 'Bragg-Brentano':
        Sample.update({'Shift':[0.,False],'Transparency':[0.,False],
                           'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
    elif NewSample.get('Type','') == 'Debye-Scherrer':
        Sample.update({'Absorption':[0.,False],'DisplaceX':[0.,False],'DisplaceY':[0.,False]})
    elif 'Lam1' in newItems:
        Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
                           'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
    else:
        Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],
                           'DisplaceY':[0.,False]})
    Sample.update(NewSample)
    return bank,[makeInstDict(newItems, newVals, len(newVals)*[False]), {}]

def WriteInstprm(fp, InstPrm, Sample={}, bank=None):
    '''Write the contents of a GSAS-II (new) .instprm instrument parameter file
    ToDo: use this inside G2frame.OnSave and G2frame.OnSaveAll

    :param file fp: Pointer to open file to be written.
    :param dict InstPrm: Instrument parameters
    :param dict Sample: Sample parameters (optional)
    :param int bank: Bank number. If not None (default), this causes
      a "#Bank" heading to be placed in the file before the
      parameters are written.
    '''
    if bank is not None:
        fp.write(f"#Bank {bank}: GSAS-II instrument parameter file; do not add/delete items!\n")
        indent = '  '
    else:
        fp.write("#GSAS-II instrument parameter file; do not add/delete items!\n")
        indent = ''
    for item in InstPrm:
        if item not in  ['XE','YE','ZE','WE']:  #skip cn to keV conversion factors
            fp.write(f"{indent}{item}:{InstPrm[item][1]}\n")
    # stick in some instrumental things that are listed in Sample
    if "Type" in Sample:
        fp.write(f"{indent}Diff-type:{Sample['Type']}\n")
    for item in ('Gonio. radius','InstrName'):
        if not Sample.get(item): continue
        fp.write(f"{indent}{item}:{Sample[item]}\n")

# version of LoadImportRoutines from before switch to "main"
# def _old_LoadImportRoutines(prefix, errprefix=None, traceback=False):
#     '''Routine to locate GSASII importers matching a prefix string.
#
#     Warns if more than one file with the same name is in the path
#     or if a file is found that is not in the main directory tree.
#     '''
#     if errprefix is None:
#         errprefix = prefix
#
#     readerlist = []
#     import_files = {}
#     if '.' not in sys.path: sys.path.append('.')
#     for path in sys.path:
#         for filename in glob.iglob(os.path.join(path, 'G2'+prefix+'*.py')):
#             pkg = os.path.splitext(os.path.split(filename)[1])[0]
#             if pkg in import_files:
#                 G2Print('Warning: importer {} overrides {}'.format(import_files[pkg],os.path.abspath(filename)))
#             elif not filename.startswith(GSASIIpath.path2GSAS2):
#                 G2Print('Note: found importer in non-standard location:'+
#                             f'\n\t{os.path.abspath(filename)}')
#                 import_files[pkg] = filename
#             else:
#                 import_files[pkg] = filename

#     for pkg in sorted(import_files.keys()):
#         try:
#             exec('import '+pkg)
#             #print(eval(pkg+'.__file__'))
#             for name, value in inspect.getmembers(eval(pkg)):
#                 if name.startswith('_'):
#                     continue
#                 if inspect.isclass(value):
#                     for method in 'Reader', 'ExtensionValidator', 'ContentsValidator':
#                         if not hasattr(value, method):
#                             break
#                         if not callable(getattr(value, method)):
#                             break
#                     else:
#                         reader = value()
#                         if reader.UseReader:
#                             readerlist.append(reader)
#         except AttributeError:
#             G2Print ('Import_' + errprefix + ': Attribute Error ' + import_files[pkg])
#             if traceback:
#                 traceback.print_exc(file=sys.stdout)
#         except Exception as exc:
#             G2Print ('\nImport_' + errprefix + ': Error importing file ' + import_files[pkg])
#             G2Print (u'Error message: {}\n'.format(exc))
#             if traceback:
#                 traceback.print_exc(file=sys.stdout)
#     return readerlist

def LoadImportRoutines(prefix, errprefix=None, traceback=False):
    from . import imports
    readerlist = []
    # how to import directly from a file for extra search magic if we need
    # https://stackoverflow.com/questions/67631/how-can-i-import-a-module-dynamically-given-the-full-path
    # TODO handle non-bundled readers
    for mod_name in (_ for _ in dir(imports) if _.startswith(f'G2{prefix}')):
        mod = getattr(imports, mod_name)
        for member_name in dir(mod):
            if member_name.startswith('_'):
                continue
            member = getattr(mod, member_name)
            if all(
                hasattr(member, meth)
                for meth in ('Reader', 'ExtensionValidator', 'ContentsValidator')
            ):
                reader = member()
                if reader.UseReader:
                    readerlist.append(reader)
    return readerlist

# version of LoadExportRoutines from before switch to "main" 
# def _LoadExportRoutines(parent, usetraceback=False):
#     '''Routine to locate GSASII exporters. Warns if more than one file
#     with the same name is in the path or if a file is found that is not
#     in the main directory tree.
#     '''
#     exporterlist = []
#     export_files = {}
#     if '.' not in sys.path: sys.path.append('.')
#     for path in sys.path:
#         for filename in glob.iglob(os.path.join(path,"G2export*.py")):
#             pkg = os.path.splitext(os.path.split(filename)[1])[0]
#             if pkg in export_files:
#                 G2Print('Warning: exporter {} overrides {}'.format(export_files[pkg],os.path.abspath(filename)))
#             elif not filename.startswith(GSASIIpath.path2GSAS2):
#                 G2Print('Note, found non-standard exporter: {}'.format(os.path.abspath(filename)))
#                 export_files[pkg] = filename
#             else:
#                 export_files[pkg] = filename
#     # go through the routines and import them, saving objects that
#     # have export routines (method Exporter)
#     for pkg in sorted(export_files.keys()):
#         try:
#             exec('import '+pkg)
#             for clss in inspect.getmembers(eval(pkg)): # find classes defined in package
#                 if clss[0].startswith('_'): continue
#                 if not inspect.isclass(clss[1]): continue
#                 # check if we have the required methods
#                 if not hasattr(clss[1],'Exporter'): continue
#                 if not callable(getattr(clss[1],'Exporter')): continue
#                 if parent is None:
#                     if not hasattr(clss[1],'Writer'): continue
#                 else:
#                     if not hasattr(clss[1],'loadParmDict'): continue
#                     if not callable(getattr(clss[1],'loadParmDict')): continue
#                 try:
#                     exporter = clss[1](parent) # create an export instance
#                 except AttributeError:
#                     pass
#                 except Exception as exc:
#                     G2Print ('\nExport init: Error substantiating class ' + clss[0])
#                     G2Print (u'Error message: {}\n'.format(exc))
#                     if usetraceback:
#                         import traceback
#                         traceback.print_exc(file=sys.stdout)
#                     continue
#                 exporterlist.append(exporter)
#         except AttributeError:
#             G2Print ('Export Attribute Error ' + export_files[pkg])
#             if usetraceback:
#                 import traceback
#                 traceback.print_exc(file=sys.stdout)
#         except Exception as exc:
#             G2Print ('\nExport init: Error importing file ' + export_files[pkg])
#             G2Print (u'Error message: {}\n'.format(exc))
#             if usetraceback:
#                 import traceback
#                 traceback.print_exc(file=sys.stdout)
#     return exporterlist

def LoadExportRoutines(parent, usetraceback=False):
    from . import exports
    exporterlist = []
    # how to import directly from a file for extra search magic if we need
    # https://stackoverflow.com/questions/67631/how-can-i-import-a-module-dynamically-given-the-full-path
    # TODO handle non-bundled readers
    for mod_name in (_ for _ in dir(exports) if _.startswith('G2export')):
        mod = getattr(exports, mod_name)
        for member_name in dir(mod):
            if member_name.startswith('_'):
                continue
            member = getattr(mod, member_name)
            if not hasattr(member, 'Exporter'):
                continue
            if parent is None:
                if not hasattr(member, 'Writer'):
                    continue
            else:
                if not hasattr(member, 'loadParmDict'):
                    continue
            try:
                exporter = member(parent)
                exporterlist.append(exporter)
            except Exception as exc:
                G2Print ('\nExport init: Error substantiating class ' + member_name)
                G2Print ('Error message: {}\n'.format(exc))
                if usetraceback:
                    import traceback
                    traceback.print_exc(file=sys.stdout)
                continue
    return exporterlist


def readColMetadata(imagefile):
    '''Reads image metadata from a column-oriented metadata table
    (1-ID style .par file). Called by :func:`GetColumnMetadata`

    The .par file has any number of columns separated by spaces.
    The directory for the file must be specified in
    Config variable :data:`config_example.Column_Metadata_directory`.
    As an index to the .par file a second "label file" must be specified with the
    same file root name as the .par file but the extension must be .XXX_lbls (where
    .XXX is the extension of the image) or if that is not present extension
    .lbls.

    :param str imagefile: the full name of the image file (with extension, directory optional)

    :returns: a dict with parameter values. Named parameters will have the type based on
       the specified Python function, named columns will be character strings

    The contents of the label file will look like this::

        # define keywords
        filename:lambda x,y: "{}_{:0>6}".format(x,y)|33,34
        distance: float | 75
        wavelength:lambda keV: 12.398425/float(keV)|9
        pixelSize:lambda x: [74.8, 74.8]|0
        ISOlikeDate: lambda dow,m,d,t,y:"{}-{}-{}T{} ({})".format(y,m,d,t,dow)|0,1,2,3,4
        Temperature: float|53
        FreePrm2: int | 34 | Free Parm2 Label
        # define other variables
        0:day
        1:month
        2:date
        3:time
        4:year
        7:I_ring

    This file contains three types of lines in any order.
     * Named parameters are evaluated with user-supplied Python code (see
       subsequent information). Specific named parameters are used
       to determine values that are used for image interpretation (see table,
       below). Any others are copied to the Comments subsection of the Image
       tree item.
     * Column labels are defined with a column number (integer) followed by
       a colon (:) and a label to be assigned to that column. All labeled
       columns are copied to the Image's Comments subsection.
     * Comments are any line that does not contain a colon.

    Note that columns are numbered starting at zero.

    Any named parameter may be defined provided it is not a valid integer,
    but the named parameters in the table have special meanings, as descibed.
    The parameter name is followed by a colon. After the colon, specify
    Python code that defines or specifies a function that will be called to
    generate a value for that parameter.

    Note that several keywords, if defined in the Comments, will be found and
    placed in the appropriate section of the powder histogram(s)'s Sample
    Parameters after an integration: ``Temperature``, ``Pressure``, ``Time``,
    ``FreePrm1``, ``FreePrm2``, ``FreePrm3``, ``Omega``, ``Chi``, and ``Phi``.

    After the Python code, supply a vertical bar (|) and then a list of one
    more more columns that will be supplied as arguments to that function.

    Note that the labels for the three FreePrm items can be changed by
    including that label as a third item with an additional vertical bar. Labels
    will be ignored for any other named parameters.

    The examples above are discussed here:

    ``filename:lambda x,y: "{}_{:0>6}".format(x,y)|33,34``
        Here the function to be used is defined with a lambda statement::

          lambda x,y: "{}_{:0>6}".format(x,y)

        This function will use the format function to create a file name from the
        contents of columns 33 and 34. The first parameter (x, col. 33) is inserted directly into
        the file name, followed by a underscore (_), followed by the second parameter (y, col. 34),
        which will be left-padded with zeros to six characters (format directive ``:0>6``).

        When there will be more than one image generated per line in the .par file, an alternate way to
        generate list of file names takes into account the number of images generated::

          lambda x,y,z: ["{}_{:0>6}".format(x,int(y)+i) for i in range(int(z))]

        Here a third parameter is used to specify the number of images generated, where
        the image number is incremented for each image.

    ``distance: float | 75``
        Here the contents of column 75 will be converted to a floating point number
        by calling float on it. Note that the spaces here are ignored.

    ``wavelength:lambda keV: 12.398425/float(keV)|9``
        Here we define an algebraic expression to convert an energy in keV to a
        wavelength and pass the contents of column 9 as that input energy

    ``pixelSize:lambda x: [74.8, 74.8]|0``
        In this case the pixel size is a constant (a list of two numbers). The first
        column is passed as an argument as at least one argument is required, but that
        value is not used in the expression.

    ``ISOlikeDate: lambda dow,m,d,t,y:"{}-{}-{}T{} ({})".format(y,m,d,t,dow)|0,1,2,3,4``
        This example defines a parameter that takes items in the first five columns
        and formats them in a different way. This parameter is not one of the pre-defined
        parameter names below. Some external code could be used to change the month string
        (argument ``m``) to a integer from 1 to 12.

    ``FreePrm2: int | 34 | Free Parm2 Label``
        In this example, the contents of column 34 will be converted to an integer and
        placed as the second free-named parameter in the Sample Parameters after an
        integration. The label for this parameter will be changed to "Free Parm2 Label".

    **Pre-defined parameter names**

    =============  =========  ========  =====================================================
     keyword       required    type      Description
    =============  =========  ========  =====================================================
       filename    yes         str or   generates the file name prefix for the matching image
                               list     file (MyImage001 for file /tmp/MyImage001.tif) or
                                        a list of file names.
     polarization  no         float     generates the polarization expected based on the
                                        monochromator angle, defaults to 0.99.
       center      no         list of   generates the approximate beam center on the detector
                              2 floats  in mm, such as [204.8, 204.8].
       distance    yes        float     generates the distance from the sample to the detector
                                        in mm
       pixelSize   no         list of   generates the size of the pixels in microns such as
                              2 floats  [200.0, 200.0].
       wavelength  yes        float     generates the wavelength in Angstroms
    =============  =========  ========  =====================================================

    '''
    dir,fil = os.path.split(os.path.abspath(imagefile))
    imageName,ext = os.path.splitext(fil)
    if not GSASIIpath.GetConfigValue('Column_Metadata_directory'): return
    parfiles = glob.glob(os.path.join(GSASIIpath.GetConfigValue('Column_Metadata_directory'),'*.par'))
    if len(parfiles) == 0:
        G2Print('Sorry, No Column metadata (.par) file found in '+
              GSASIIpath.GetConfigValue('Column_Metadata_directory'))
        return {}
    for parFil in parfiles: # loop over all .par files (hope just 1) in image dir until image is found
        parRoot = os.path.splitext(parFil)[0]
        for e in (ext+'_lbls','.lbls'):
            if os.path.exists(parRoot+e):
                lblFil = parRoot+e
                break
        else:
            G2Print('Warning: No labels definitions found for '+parFil)
            continue
        labels,lbldict,keyCols,keyExp,errors = readColMetadataLabels(lblFil)
        if errors:
            print('Errors in labels file '+lblFil)
            for i in errors: print('  '+i)
            continue
        else:
            G2Print('Read '+lblFil)
        # scan through each line in this .par file, looking for the matching image rootname
        fp = open(parFil,'r')
        for iline,line in enumerate(fp):
            items = line.strip().split(' ')
            nameList = keyExp['filename'](*[items[j] for j in keyCols['filename']])
            if type(nameList) is str:
                if nameList != imageName: continue
                name = nameList
            else:
                for name in nameList:
                    if name == imageName: break # got a match
                else:
                    continue
            # parse the line and finish
            metadata = evalColMetadataDicts(items,labels,lbldict,keyCols,keyExp)
            metadata['par file'] = parFil
            metadata['lbls file'] = lblFil
            G2Print("Metadata read from {} line {}".format(parFil,iline+1))
            fp.close()
            return metadata
        else:
            G2Print("Image {} not found in {}".format(imageName,parFil))
            fp.close()
            continue
        fp.close()
    else:
        G2Print("Warning: No .par metadata for image {}".format(imageName))
        return {}

def readColMetadataLabels(lblFil):
    '''Read the .*lbls file and setup for metadata assignments
    '''
    lbldict = {}
    keyExp = {}
    keyCols = {}
    labels = {}
    errors = []
    fp = open(lblFil,'r')         # read column labels
    for iline,line in enumerate(fp): # read label definitions
        line = line.strip()
        if not line or line[0] == '#': continue # comments
        items = line.split(':')
        if len(items) < 2: continue # lines with no colon are also comments
        # does this line a definition for a named parameter?
        key = items[0]
        try:
            int(key)
        except ValueError: # try as named parameter since not a valid number
            items = line.split(':',1)[1].split('|')
            try:
                f = eval(items[0]) # compile the expression
                if not callable(f):
                    errors += ['Expression "{}" for key {} is not a function (line {})'.
                           format(items[0],key,iline)]
                    continue
                keyExp[key] = f
            except Exception as msg:
                errors += ['Expression "{}" for key {} is not valid (line {})'.
                           format(items[0],key,iline)]
                errors += [str(msg)]
                continue
            keyCols[key] = [int(i) for i in items[1].strip().split(',')]
            if key.lower().startswith('freeprm') and len(items) > 2:
                labels[key] = items[2]
            continue
        if len(items) == 2: # simple column definition
            lbldict[int(items[0])] = items[1]
    fp.close()
    if 'filename' not in keyExp:
        errors += ["File {} is invalid. No valid filename expression.".format(lblFil)]
    return labels,lbldict,keyCols,keyExp,errors

def evalColMetadataDicts(items,labels,lbldict,keyCols,keyExp,ShowError=False):
    '''Evaluate the metadata for a line in the .par file
    '''
    metadata = {lbldict[j]:items[j] for j in lbldict}
    named = {}
    for key in keyExp:
        try:
            res = keyExp[key](*[items[j] for j in keyCols[key]])
        except:
            if ShowError:
                res = "*** error ***"
            else:
                continue
        named[key] = res
    metadata.update(named)
    for lbl in labels: # add labels for FreePrm's
        metadata['label_'+lbl[4:].lower()] = labels[lbl]
    return metadata

def GetColumnMetadata(reader):
    '''Add metadata to an image from a column-type metadata file
    using :func:`readColMetadata`

    :param reader: a reader object from reading an image

    '''
    if not GSASIIpath.GetConfigValue('Column_Metadata_directory'): return
    parParms = readColMetadata(reader.readfilename)
    if not parParms: return # check for read failure
    specialKeys = ('filename',"polarization", "center", "distance", "pixelSize", "wavelength",)
    reader.Comments = ['Metadata from {} assigned by {}'.format(parParms['par file'],parParms['lbls file'])]
    for key in parParms:
        if key in specialKeys+('par file','lbls file'): continue
        reader.Comments += ["{} = {}".format(key,parParms[key])]
    if "polarization" in parParms:
        reader.Data['PolaVal'][0] = parParms["polarization"]
    else:
        reader.Data['PolaVal'][0] = 0.99
    if "center" in parParms:
        reader.Data['center'] = parParms["center"]
    if "pixelSize" in parParms:
        reader.Data['pixelSize'] = parParms["pixelSize"]
    if "wavelength" in parParms:
        reader.Data['wavelength'] = parParms['wavelength']
    else:
        G2Print('Error: wavelength not defined in {}'.format(parParms['lbls file']))
    if "distance" in parParms:
        reader.Data['distance'] = parParms['distance']
        reader.Data['setdist'] = parParms['distance']
    else:
        G2Print('Error: distance not defined in {}'.format(parParms['lbls file']))

def LoadControls(Slines,data):
    'Read values from a .imctrl (Image Controls) file'
    cntlList = ['color','wavelength','distance','tilt','invert_x','invert_y','type','Oblique',
        'fullIntegrate','outChannels','outAzimuths','LRazimuth','IOtth','azmthOff','DetDepth',
        'calibskip','pixLimit','cutoff','calibdmin','Flat Bkg','varyList','setdist',
        'PolaVal','SampleAbs','dark image','background image','twoth']
    save = {}
    for S in Slines:
        if S[0] == '#':
            continue
        [key,val] = S.strip().split(':',1)
        if key in ['type','calibrant','binType','SampleShape','color',]:    #strings
            save[key] = val
        elif key in ['varyList',]:
            save[key] = eval(val)   #dictionary
        elif key in ['rotation']:
            save[key] = float(val)
        elif key in ['center',]:
            if ',' in val:
                save[key] = eval(val)
            else:
                vals = val.strip('[] ').split()
                save[key] = [float(vals[0]),float(vals[1])]
        elif key in cntlList:
            save[key] = eval(val)
    data.update(save)

def WriteControls(filename,data):
    'Write current values to a .imctrl (Image Controls) file'
    File = open(filename,'w')
    keys = ['type','color','wavelength','calibrant','distance','center','Oblique',
            'tilt','rotation','azmthOff','fullIntegrate','LRazimuth','setdist',
            'IOtth','outChannels','outAzimuths','invert_x','invert_y','DetDepth',
            'calibskip','pixLimit','cutoff','calibdmin','Flat Bkg','varyList',
            'binType','SampleShape','PolaVal','SampleAbs','dark image','background image',
            'twoth']
    for key in keys:
        if key not in data:     #uncalibrated!
            continue
        File.write(key+':'+str(data[key])+'\n')
    File.close()

def GetImageData(G2frame,imagefile,imageOnly=False,ImageTag=None,FormatName=''):
    '''Read a single image with an image importer. This is called to reread an image
    after it has already been imported with :meth:`GSASIIdataGUI.GSASII.OnImportGeneric`
    (or :func:`GSASIImiscGUI.ReadImages` in Auto Integration) so it is not necessary to reload metadata.

    :param wx.Frame G2frame: main GSAS-II Frame and data object.
    :param str imagefile: name of image file
    :param bool imageOnly: If True return only the image. Formerly, if False
      return more information (see below). Now must be True
    :param int/str ImageTag: specifies a particular image to be read from a file.
      First image is read if None (default).
    :param str formatName: the image reader formatName

    :returns: an image as a numpy array.
      Formerly if imageOnly=False this would return a list of four items: Comments, Data, Npix and the Image,
      but now an exception occurs when imageOnly=False.
    '''
    # determine which formats are compatible with this file
    primaryReaders = []
    secondaryReaders = []
    for rd in G2frame.ImportImageReaderlist:
        flag = rd.ExtensionValidator(imagefile)
        if flag is None:
            secondaryReaders.append(rd)
        elif flag:
            if not FormatName:
                primaryReaders.append(rd)
            elif FormatName in rd.formatName:       #This is a kluge because the rd.formatName was changed!
                primaryReaders.append(rd)
    if len(secondaryReaders) + len(primaryReaders) == 0:
        print('Error: No matching format for file '+imagefile)
        raise Exception('No image read')
    errorReport = ''
    if not imagefile:
        return
    for rd in primaryReaders+secondaryReaders:
        rd.ReInitialize() # purge anything from a previous read
        rd.errors = "" # clear out any old errors
        if not rd.ContentsValidator(imagefile): # rejected on cursory check
            errorReport += "\n  "+rd.formatName + ' validator error'
            if rd.errors:
                errorReport += ': '+rd.errors
                continue
        if imageOnly:
            ParentFrame = None # prevent GUI access on reread
        else:
            ParentFrame = G2frame
        if GSASIIpath.GetConfigValue('debug'):
            flag = rd.Reader(imagefile,ParentFrame,blocknum=ImageTag)
        else:
            flag = False
            try:
                flag = rd.Reader(imagefile,ParentFrame,blocknum=ImageTag)
            except rd.ImportException as detail:
                rd.errors += "\n  Read exception: "+str(detail)
            except Exception as detail:
                import traceback
                rd.errors += "\n  Unhandled read exception: "+str(detail)
                rd.errors += "\n  Traceback info:\n"+str(traceback.format_exc())
        if flag: # this read succeeded
            if rd.Image is None:
                raise Exception('No image read. Strange!')
            if GSASIIpath.GetConfigValue('Transpose'):
                print ('Transposing Image!')
                rd.Image = rd.Image.T
            #rd.readfilename = imagefile
            if imageOnly:
                return rd.Image
            else:
                raise Exception('GetImageData must be called with imageOnly=True')
                #return rd.Comments,rd.Data,rd.Npix,rd.Image
    else:
        print('Error reading file '+imagefile)
        print('Error messages(s)\n'+errorReport)
        raise Exception('No image read')

def RereadImageData(ImageReaderlist,imagefile,ImageTag=None,FormatName=''):
    '''Read a single image with an image importer. This is called to
    reread an image after it has already been imported, so it is not
    necessary to reload metadata.

    Based on :func:`GetImageData` which this can replace
    where imageOnly=True

    :param list ImageReaderlist: list of Reader objects for images
    :param str imagefile: name of image file
    :param int/str ImageTag: specifies a particular image to be read from a file.
      First image is read if None (default).
    :param str formatName: the image reader formatName

    :returns: an image as a numpy array
    '''
    # determine which formats are compatible with this file
    primaryReaders = []
    secondaryReaders = []
    for rd in ImageReaderlist:
        flag = rd.ExtensionValidator(imagefile)
        if flag is None:
            secondaryReaders.append(rd)
        elif flag:
            if not FormatName:
                primaryReaders.append(rd)
            elif FormatName == rd.formatName:
                primaryReaders.append(rd)
    if len(secondaryReaders) + len(primaryReaders) == 0:
        G2Print('Error: No matching format for file '+imagefile)
        raise Exception('No image read')
    errorReport = ''
    if not imagefile:
        return
    for rd in primaryReaders+secondaryReaders:
        rd.ReInitialize() # purge anything from a previous read
        rd.errors = "" # clear out any old errors
        if not rd.ContentsValidator(imagefile): # rejected on cursory check
            errorReport += "\n  "+rd.formatName + ' validator error'
            if rd.errors:
                errorReport += ': '+rd.errors
                continue
        flag = rd.Reader(imagefile,None,blocknum=ImageTag)
        if flag: # this read succeeded
            if rd.Image is None:
                raise Exception('No image read. Strange!')
            if GSASIIpath.GetConfigValue('Transpose'):
                G2Print ('Warning: Transposing Image!')
                rd.Image = rd.Image.T
            #rd.readfilename = imagefile
            return rd.Image
    else:
        G2Print('Error reading file '+imagefile)
        G2Print('Error messages(s)\n'+errorReport)
        raise Exception('No image read')

def readMasks(filename,masks,ignoreThreshold):
    '''Read a GSAS-II masks file'''
    File = open(filename,'r')
    save = {}
    oldThreshold = masks['Thresholds'][0]
    S = File.readline()
    while S:
        if S[0] == '#':
            S = File.readline()
            continue
        [key,val] = S.strip().split(':',1)
        if key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds','Xlines','Ylines']:
            if ignoreThreshold and key == 'Thresholds':
                S = File.readline()
                continue
            save[key] = eval(val)
            if key == 'Thresholds':
                save[key][0] = oldThreshold
                save[key][1][1] = min(oldThreshold[1],save[key][1][1])
        S = File.readline()
    File.close()
    masks.update(save)
    # CleanupMasks
    for key in ['Points','Rings','Arcs','Polygons']:
        masks[key] = masks.get(key,[])
        masks[key] = [i for i in masks[key] if len(i)]

def PDFWrite(PDFentry,fileroot,PDFsaves,PDFControls,Inst={},Limits=[]):
    '''Write PDF-related data (G(r), S(Q),...) into files, as
    selected.

    :param str PDFentry: name of the PDF entry in the tree. This is
      used for comments in the file specifying where it came from;
      it can be arbitrary
    :param str fileroot: name of file(s) to be written. The extension
      will be ignored.
    :param list PDFsaves: flags that determine what type of file will be
      written:
      PDFsaves[0], if True writes a I(Q) file with a .iq extension
      PDFsaves[1], if True writes a S(Q) file with a .sq extension
      PDFsaves[2], if True writes a F(Q) file with a .fq extension
      PDFsaves[3], if True writes a G(r) file with a .gr extension
      PDFsaves[4], if True writes G(r) in a pdfGUI input file with
      a .gr extension. Note that if PDFsaves[3] and PDFsaves[4] are
      both True, the pdfGUI overwrites the G(r) file.
      PDFsaves[5], if True writes F(Q) & g(R) with .fq & .gr extensions
      overwrites these if selected by option 2, 3 or 4
    :param dict PDFControls: The PDF parameters and computed results
    :param dict Inst: Instrument parameters from the PDWR entry used
      to compute the PDF. Needed only when PDFsaves[4] is True.
    :param list Limits: Computation limits from the PDWR entry used
      to compute the PDF. Needed only when PDFsaves[4] is True.
    '''
    import scipy.interpolate as scintp
    fileroot = os.path.splitext(fileroot)[0]
    if PDFsaves[0]:     #I(Q)
        iqfilename = fileroot+'.iq'
        iqdata = PDFControls['I(Q)'][1]
        iqfxn = scintp.interp1d(iqdata[0],iqdata[1],kind='linear')
        iqfile = open(iqfilename,'w')
        iqfile.write('#T I(Q) %s\n'%(PDFentry))
        iqfile.write('#L Q     I(Q)\n')
        qnew = np.arange(iqdata[0][0],iqdata[0][-1],0.005)
        iqnew = zip(qnew,iqfxn(qnew))
        for q,iq in iqnew:
            iqfile.write("%15.6g %15.6g\n" % (q,iq))
        iqfile.close()
        G2Print (' I(Q) saved to: '+iqfilename)

    if PDFsaves[1]:     #S(Q)
        sqfilename = fileroot+'.sq'
        sqdata = PDFControls['S(Q)'][1]
        sqfxn = scintp.interp1d(sqdata[0],sqdata[1],kind='linear')
        sqfile = open(sqfilename,'w')
        sqfile.write('#T S(Q) %s\n'%(PDFentry))
        sqfile.write('#L Q     S(Q)\n')
        qnew = np.arange(sqdata[0][0],sqdata[0][-1],0.005)
        sqnew = zip(qnew,sqfxn(qnew))
        for q,sq in sqnew:
            sqfile.write("%15.6g %15.6g\n" % (q,sq))
        sqfile.close()
        G2Print (' S(Q) saved to: '+sqfilename)

    if PDFsaves[2]:     #F(Q)
        fqfilename = fileroot+'.fq'
        fqdata = PDFControls['F(Q)'][1]
        fqfxn = scintp.interp1d(fqdata[0],fqdata[1],kind='linear')
        fqfile = open(fqfilename,'w')
        fqfile.write('#T F(Q) %s\n'%(PDFentry))
        fqfile.write('#L Q     F(Q)\n')
        qnew = np.arange(fqdata[0][0],fqdata[0][-1],0.005)
        fqnew = zip(qnew,fqfxn(qnew))
        for q,fq in fqnew:
            fqfile.write("%15.6g %15.6g\n" % (q,fq))
        fqfile.close()
        G2Print (' F(Q) saved to: '+fqfilename)

    if PDFsaves[3]:     #G(R)
        grfilename = fileroot+'.gr'
        grdata = PDFControls['G(R)'][1]
        grfxn = scintp.interp1d(grdata[0],grdata[1],kind='linear')
        grfile = open(grfilename,'w')
        grfile.write('#T G(R) %s\n'%(PDFentry))
        grfile.write('#L R     G(R)\n')
        rnew = np.arange(grdata[0][0],grdata[0][-1],0.010)
        grnew = zip(rnew,grfxn(rnew))
        for r,gr in grnew:
            grfile.write("%15.6g %15.6g\n" % (r,gr))
        grfile.close()
        G2Print (' G(R) saved to: '+grfilename)

    if PDFsaves[4]: #pdfGUI file for G(R)
        grfilename = fileroot+'.gr'
        grdata = PDFControls['G(R)'][1]
        qdata = PDFControls['I(Q)'][1][0]
        grfxn = scintp.interp1d(grdata[0],grdata[1],kind='linear')
        grfile = open(grfilename,'w')
        rnew = np.arange(grdata[0][0],grdata[0][-1],0.010)
        grnew = zip(rnew,grfxn(rnew))

        grfile.write('[DEFAULT]\n')
        grfile.write('\n')
        grfile.write('version = GSAS-II-v'+str(GSASIIpath.GetVersionNumber())+'\n')
        grfile.write('\n')
        grfile.write('# input and output specifications\n')
        grfile.write('dataformat = Qnm\n')
        grfile.write('inputfile = %s\n'%(PDFControls['Sample']['Name']))
        grfile.write('backgroundfile = %s\n'%(PDFControls['Sample Bkg.']['Name']))
        grfile.write('outputtype = gr\n')
        grfile.write('\n')
        grfile.write('# PDF calculation setup\n')
        if 'x' in Inst['Type']:
            grfile.write('mode = %s\n'%('xray'))
        elif 'N' in Inst['Type']:
            grfile.write('mode = %s\n'%('neutron'))
        wave = G2mth.getMeanWave(Inst)
        grfile.write('wavelength = %.5f\n'%(wave))
        formula = ''
        for el in PDFControls['ElList']:
            formula += el
            num = PDFControls['ElList'][el]['FormulaNo']
            if num == round(num):
                formula += '%d'%(int(num))
            else:
                formula += '%.2f'%(num)
        grfile.write('composition = %s\n'%(formula))
        grfile.write('bgscale = %.3f\n'%(-PDFControls['Sample Bkg.']['Mult']))
        highQ = 2.*np.pi/G2lat.Pos2dsp(Inst,Limits[1][1])
        grfile.write('qmaxinst = %.2f\n'%(highQ))
        grfile.write('qmin = %.5f\n'%(qdata[0]))
        grfile.write('qmax = %.4f\n'%(qdata[-1]))
        grfile.write('rmin = %.2f\n'%(PDFControls['Rmin']))
        grfile.write('rmax = %.2f\n'%(PDFControls['Rmax']))
        grfile.write('rstep = 0.01\n')
        grfile.write('\n')
        grfile.write('# End of config '+63*'-')
        grfile.write('\n')
        grfile.write('#### start data\n')
        grfile.write('#S 1\n')
        grfile.write('#L r($\\AA$)  G($\\AA^{-2}$)\n')
        for r,gr in grnew:
            grfile.write("%15.2F %15.6F\n" % (r,gr))
        grfile.close()
        G2Print (' G(R) saved to: '+grfilename)

    if len(PDFsaves) > 5 and PDFsaves[5]: #RMCProfile files for F(Q) & g(r) overwrites any above

        fqfilename = fileroot+'.fq'
        fqdata = PDFControls['F(Q)'][1]
        fqfxn = scintp.interp1d(fqdata[0],fqdata[1],kind='linear')
        fqfile = open(fqfilename,'w')
        qnew = np.arange(fqdata[0][0],fqdata[0][-1],0.005)
        nq = qnew.shape[0]
        fqfile.write('%20d\n'%(nq-1))
        fqfile.write(fqfilename+'\n')
        fqnew = zip(qnew,fqfxn(qnew))
        for q,fq in fqnew:
            if not q:
                continue
            fqfile.write("%15.6g %15.6g\n" % (q,fq))
        fqfile.close()
        G2Print (' F(Q) saved to: '+fqfilename)

        grfilename = fileroot+'.gr'
        grdata = PDFControls['g(r)'][1]
        grfxn = scintp.interp1d(grdata[0],grdata[1],kind='linear')
        grfile = open(grfilename,'w')
        rnew = np.arange(grdata[0][0],grdata[0][-1],0.010)
        nr = rnew.shape[0]
        grfile.write('%20d\n'%(nr-1))
        grfile.write(grfilename+'\n')
        grnew = zip(rnew,grfxn(rnew))
        for r,gr in grnew:
            if not r:
                continue
            grfile.write("%15.6g %15.6g\n" % (r,gr))
        grfile.close()
        G2Print (' G(R) saved to: '+grfilename)

#===========================================================================
#-- Routines for formatting output originally in module GSASIIpy3
#
# formatting for unique cell parameters by Laue type
cellGUIlist = [
    [['m3','m3m'],4,[" Unit cell: a = "],["{:.5f}"],[True],[0]],
    [['3R','3mR'],6,[" a = ",u" \u03B1 = "],["{:.5f}","{:.3f}"],[True,True],[0,3]],
    [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],6,[" a = "," c = "],["{:.5f}","{:.5f}"],[True,True],[0,2]],
    [['mmm'],8,[" a = "," b = "," c = "],["{:.5f}","{:.5f}","{:.5f}"],
        [True,True,True],[0,1,2]],
    [['2/m'+'a'],10,[" a = "," b = "," c = ",u" \u03B1 = "],
        ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[True,True,True,True,],[0,1,2,3]],
    [['2/m'+'b'],10,[" a = "," b = "," c = ",u" \u03B2 = "],
        ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[True,True,True,True,],[0,1,2,4]],
    [['2/m'+'c'],10,[" a = "," b = "," c = ",u" \u03B3 = "],
        ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[True,True,True,True,],[0,1,2,5]],
    [['-1'],7,[" a = "," b = "," c = ",u" \u03B1 = ",u" \u03B2 = ",u" \u03B3 = "],
         ["{:.5f}","{:.5f}","{:.5f}","{:.3f}","{:.3f}","{:.3f}"],
         [True,True,True,True,True,True],[0,1,2,3,4,5]]
    ]

def FormulaEval(string):
    '''Evaluates a algebraic formula into a float, if possible. Works
    properly on fractions e.g. 2/3 only with python 3.0+ division.

    Expressions such as 2/3, 3*pi, sin(45)/2, 2*sqrt(2), 2**10 can all
    be evaluated.

    :param str string: Character string containing a Python expression
      to be evaluated.

    :returns: the value for the expression as a float or None if the expression does not
      evaluate to a valid number.

    '''
    try:
        val = float(eval(string))
        if np.isnan(val) or np.isinf(val): return None
    except:
        return None
    return val

def FormatPadValue(val,maxdigits=None):
    '''Format a float to fit in ``maxdigits[0]`` spaces with maxdigits[1] after decimal.

    :param float val: number to be formatted.

    :param list maxdigits: the number of digits & places after decimal to be used for display of the
      number (defaults to [10,2]).

    :returns: a string with exactly maxdigits[0] characters (except under error conditions),
      but last character will always be a space
    '''
    if maxdigits is None:
        digits = [10,2]
    else:
        digits = list(maxdigits)
    fmt = '{:'+str(digits[0])+'}'
    s = fmt.format(FormatValue(val,digits))
    if s[-1] == ' ':
        return s
    else:
        return s+' '


def FormatValue(val,maxdigits=None):
    '''Format a float to fit in at most ``maxdigits[0]`` spaces with maxdigits[1] after decimal.
    Note that this code has been hacked from FormatSigFigs and may have unused sections.

    :param float val: number to be formatted.

    :param list maxdigits: the number of digits, places after decimal and 'f' or 'g' to be used for display of the
      number (defaults to [10,2,'f']).

    :returns: a string with <= maxdigits characters (usually).
    '''
    if 'str' in str(type(val)) and (val == '?' or val == '.'):
        return val
    if maxdigits is None:
        digits = [10,2,'f']
    else:
        digits = list(maxdigits)
    if len(digits) == 2:
        digits.append('f')
    if not val:
        digits[2] = 'f'
    fmt="{:"+str(digits[0])+"."+str(digits[1])+digits[2]+"}"
    string = fmt.format(float(val)).strip() # will standard .f formatting work?
    if len(string) <= digits[0]:
        if ':' in string: # deal with weird bug where a colon pops up in a number when formatting (EPD 7.3.2!)
            string = str(val)
        if digits[1] > 0 and not 'e' in string.lower(): # strip off extra zeros on right side
            string = string.rstrip('0')
            if string[-1] == '.': string += "0"
        return string
    if val < 0: digits[0] -= 1 # negative numbers, reserve space for the sign
    decimals = digits[0] - digits[1]
    if abs(val) > 1e99: # for very large numbers, use scientific notation and use all digits
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],digits[0]-6))+"}"
    elif abs(val) > 1e9:
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],digits[0]-5))+"}"
    elif abs(val) < 10**(4-decimals): # make sure at least 4 decimals show
        # this clause is probably no longer needed since the number probably shows as 0.0
        decimals = min(digits[0]-5,digits[1])
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    elif abs(val) >= 10**(decimals-1): # deal with large numbers in smaller spaces
        decimals = max(0,digits[0]-5)
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    elif abs(val) < 1: # use f format for small numbers
        # this clause is probably no longer needed since the number probably shows as 0.0
        decimals = min(digits[0]-3,digits[1])
        fmt = "{" + (":{:d}.{:d}f".format(digits[0],decimals))+"}"
    else: # in range where g formatting should do what I want
        # used?
        decimals = digits[0] - 6
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    try:
        return fmt.format(float(val)).strip()
    except ValueError:
        print ('FormatValue Error with val,maxdigits,fmt= %f %d %s'%(val,maxdigits,fmt))
        return str(val)

def FormatSigFigs(val, maxdigits=10, sigfigs=5, treatAsZero=1e-20):
    '''Format a float to use ``maxdigits`` or fewer digits with ``sigfigs``
    significant digits showing (if room allows).

    :param float val: number to be formatted.

    :param int maxdigits: the number of digits to be used for display of the
       number (defaults to 10).

    :param int sigfigs: the number of significant figures to use, if room allows

    :param float treatAsZero: numbers that are less than this in magnitude
      are treated as zero. Defaults to 1.0e-20, but this can be disabled
      if set to None.

    :returns: a string with <= maxdigits characters (I hope).
    '''
    if 'str' in str(type(val)) and (val == '?' or val == '.'):
        return val
    if treatAsZero is not None:
        if abs(val) < treatAsZero:
            return '0.0'
    # negative numbers, leave room for a sign
    if np.isnan(val):
        return str(val)
    if val < 0: maxdigits -= 1
    if abs(val) < 1e-99 or abs(val) > 9.999e99:
        decimals = min(maxdigits-6,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}" # create format string
    elif abs(val) < 1e-9 or abs(val) > 9.999e9:
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) < 9.9999999*10**(sigfigs-maxdigits):
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) >= 10**sigfigs: # deal with large numbers in smaller spaces
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) < 1: # small numbers, add to decimal places
        decimals = sigfigs - int(np.log10(np.abs(val)))
        fmt = "{" + (":{:d}.{:d}f".format(maxdigits,decimals))+"}"
    else: # larger numbers, remove decimal places
        decimals = sigfigs - 1 - int(np.log10(np.abs(val)))
        if decimals <= 0:
            fmt = "{" + (":{:d}.0f".format(maxdigits))+"}."
        else:
            fmt = "{" + (":{:d}.{:d}f".format(maxdigits,decimals))+"}"
    try:
        return fmt.format(float(val)).strip()
    except ValueError:
        print ('FormatValue Error with val,maxdigits, sigfigs, fmt=%f %d %d %s'%(val, maxdigits,sigfigs, fmt))
        return str(val)

#===========================================================================
def openInNewTerm(project=None,g2script=None,pythonapp=sys.executable):
    '''Open a new and independent GSAS-II session in separate terminal
    or console window and as a separate process that will continue
    even if the calling process exits.
    Intended to work on all platforms.

    This could be used to run other scripts inside python other than GSAS-II

    :param str project: the name of an optional parameter to be
      passed to the script (usually a .gpx file to be opened in
      a new GSAS-II session)
    :param str g2script: the script to be run. If None (default)
      the G2.py file in the same directory as this file will
      be used.
    :param str pythonapp: the Python interpreter to be used.
      Defaults to sys.executable which is usually what is wanted.
    :param str terminal: a name for a preferred terminal emulator
    '''
    import subprocess
    if g2script is None:
        g2script = os.path.join(os.path.dirname(__file__),'G2.py')

    if sys.platform == "darwin":
        if project:
            script = f'''
set python to "{pythonapp}"
set appwithpath to "{g2script}"
set filename to "{project}"
set filename to the quoted form of the POSIX path of filename

tell application "Terminal"
     activate
     do script python & " " & appwithpath & " " & filename & "; exit"
end tell
'''
        else:
            script = f'''
set python to "{pythonapp}"
set appwithpath to "{g2script}"

tell application "Terminal"
     activate
     do script python & " " & appwithpath & " " & "; exit"
end tell
'''
        subprocess.Popen(["osascript","-e",script])
    elif sys.platform.startswith("win"):
        cmds = [pythonapp, g2script]
        if project: cmds += [project]
        subprocess.Popen(cmds,creationflags=subprocess.CREATE_NEW_CONSOLE)
    else:
        import shutil
        script = ''
        # try a bunch of common terminal emulators in Linux
        # there does not appear to be a good way to way to specify this
        # perhaps this should be a GSAS-II config option
        for term in ("lxterminal", "gnome-terminal", 'konsole', "xterm",
                         "terminator", "terminology", "tilix"):
            try:
                found = shutil.which(term)
                if not found: continue
            except AttributeError:
                print(f'shutil.which() failed (why?); assuming {term} present')
                found = True
            if term == "gnome-terminal":
                #terminal = 'gnome-terminal -t "GSAS-II console" --'
                cmds = [term,'--title','"GSAS-II console"','--']
                script = "echo; echo Press Enter to close window; read line"
                break
            elif term == "lxterminal":
               #terminal = 'lxterminal -t "GSAS-II console" -e'
               cmds = [term,'-t','"GSAS-II console"','-e']
               script = "echo;echo Press Enter to close window; read line"
               break
            elif term == "xterm":
                #terminal = 'xterm -title "GSAS-II console" -hold -e'
                cmds = [term,'-title','"GSAS-II console"','-hold','-e']
                script = "echo; echo This window can now be closed"
                break
            elif term == "terminator":
                cmds = [term,'-T','"GSAS-II console"','-x']
                script = "echo;echo Press Enter to close window; read line"
                break
            elif term == "konsole":
                cmds = [term,'-p','tabtitle="GSAS-II console"','--hold','-e']
                script = "echo; echo This window can now be closed"
                break
            elif term == "tilix":
                cmds = [term,'-t','"GSAS-II console"','-e']
                script = "echo;echo Press Enter to close window; read line"
                break
            elif term == "terminology":
                cmds = [term,'-T="GSAS-II console"','--hold','-e']
                script = "echo; echo This window can now be closed"
                break
        else:
            print("No known terminal was found to use, Can't start {}")
            return

        fil = '/tmp/GSAS2-launch.sh'
        cmds += ['/bin/sh',fil]
        fp = open(fil,'w')
        if project:
            fp.write(f"{pythonapp} {g2script} {project}\n")
        else:
            fp.write(f"{pythonapp} {g2script}\n")
        fp.write(f"rm {fil}\n")
        if script:
            fp.write(f"{script}\n")
        fp.close()
        subprocess.Popen(cmds,start_new_session=True)

def CleanupFromZip(label,cleanupList):
    '''Delete files extracted from a zip archive, typically created with
    :func:`GSASIIctrl.ExtractFileFromZip` during the data import process.

    Note that image files should not be deleted as they will be reused
    every time the image is displayed, but everything else will be saved
    in the data tree and the file copy is not needed.

    :param str label: for imports, this is the type of file being read
    :param list cleanupList: a list of files that have been extracted from
      the zip archive and can be deleted.
    '''
    if not cleanupList: return
    if 'image' in label:
        print("images don't get removed. Retaining zip-extracted file(s):")
        print('\t','\n\t'.join(cleanupList))
        return
    else:
        print("Zip-extracted file(s) will be deleted:")
        for f in cleanupList:
            try:
                os.remove(f)
                print('\tdeleted:',f)
            except:
                print('\tdelete failed:',f)

def striphist(var,insChar=''):
    'strip a histogram number from a var name'
    sv = var.split(':')
    if len(sv) <= 1: return var
    if sv[1]:
        sv[1] = insChar
    return ':'.join(sv)

def trim(val):
    '''Simplify a string containing leading and trailing spaces
    as well as newlines, tabs, repeated spaces etc. into a shorter and
    more simple string, by replacing all ranges of whitespace
    characters with a single space.

    :param str val: the string to be simplified

    :returns: the (usually) shortened version of the string
    '''
    return re.sub('\\s+', ' ', val).strip()

######################################################################
# base class for reading various types of data files
#   not used directly, only by subclassing
######################################################################
class ExportBaseclass(object):
    '''Defines a base class for the exporting of GSAS-II results.

    This class is subclassed in the various exports/G2export_*.py files. Those files
    are imported in :meth:`GSASIIdataGUI.GSASII._init_Exports` which defines the
    appropriate menu items for each one and the .Exporter method is called
    directly from the menu item.

    Routines may also define a .Writer method, which is used to write a single
    file without invoking any GUI objects.
    '''
    # TODO: review exporters producing exceptions where .Writer can't be used where G2frame is None (see CIF)
    # TODO: review conflicting uses of .Writer with mode (SeqRef) & elsewhere
    # TODO: move this class to G2fil
    def __init__(self,G2frame,formatName,extension,longFormatName=None,):
        self.G2frame = G2frame
        self.formatName = formatName # short string naming file type
        self.extension = extension
        if longFormatName: # longer string naming file type
            self.longFormatName = longFormatName
        else:
            self.longFormatName = formatName
        self.OverallParms = {}
        self.Phases = {}
        self.Histograms = {}
        self.powderDict = {}
        self.sasdDict = {}
        self.refdDict = {}
        self.xtalDict = {}
        self.parmDict = {}
        self.sigDict = {}
        self.fp = None
        # updated in InitExport:
        self.currentExportType = None # type of export that has been requested
        # updated in ExportSelect (when used):
        self.phasenam = None # a list of selected phases
        self.histnam = None # a list of selected histograms
        self.filename = None # name of file to be written (single export) or template (multiple files)
        self.dirname = '' # name of directory where file(s) will be written
        self.fullpath = '' # name of file being written -- full path

        # items that should be defined in a subclass of this class
        self.exporttype = []  # defines the type(s) of exports that the class can handle.
        # The following types are defined: 'project', "phase", "powder", "single"
        self.multiple = False # set as True if the class can export multiple phases or histograms
        # self.multiple is ignored for "project" exports

    def InitExport(self,event):
        '''Determines the type of menu that called the Exporter and
        misc initialization.
        '''
        self.filename = None # name of file to be written (single export)
        self.dirname = '' # name of file to be written (multiple export)
        if event:
            self.currentExportType = self.G2frame.ExportLookup.get(event.Id)

    def MakePWDRfilename(self,hist):
        '''Make a filename root (no extension) from a PWDR histogram name

        :param str hist: the histogram name in data tree (starts with "PWDR ")
        '''
        file0 = ''
        file1 = hist[5:]
        # replace repeated blanks
        while file1 != file0:
            file0 = file1
            file1 = file0.replace('  ',' ').strip()
        file0 = file1.replace('Azm= ','A')
        # if angle has unneeded decimal places on aziumuth, remove them
        if file0[-3:] == '.00': file0 = file0[:-3]
        file0 = file0.replace('.','_')
        file0 = file0.replace(' ','_')
        return file0

    def ExportSelect(self,AskFile='ask'):
        '''Selects histograms or phases when needed. Sets a default file name when
        requested into self.filename; always sets a default directory in self.dirname.

        :param bool AskFile: Determines how this routine processes getting a
          location to store the current export(s).

          * if AskFile is 'ask' (default option), get the name of the file to be written;
            self.filename and self.dirname are always set. In the case where
            multiple files must be generated, the export routine should do this
            based on self.filename as a template.
          * if AskFile is 'dir', get the name of the directory to be used;
            self.filename is not used, but self.dirname is always set. The export routine
            will always generate the file name.
          * if AskFile is 'single', get only the name of the directory to be used when
            multiple items will be written (as multiple files) are used
            *or* a complete file name is requested when a single file
            name is selected. self.dirname is always set and self.filename used
            only when a single file is selected.
          * if AskFile is 'default', creates a name of the file to be used from
            the name of the project (.gpx) file. If the project has not been saved,
            then the name of file is requested.
            self.filename and self.dirname are always set. In the case where
            multiple file names must be generated, the export routine should do this
            based on self.filename.
          * if AskFile is 'default-dir', sets self.dirname from the project (.gpx)
            file. If the project has not been saved, then a directory is requested.
            self.filename is not used.

        :returns: True in case of an error
        '''

        numselected = 1
        from . import GSASIIctrlGUI as G2G
        if self.currentExportType == 'phase':
            if len(self.Phases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any phases.')
                return True
            elif len(self.Phases) == 1:
                self.phasenam = list(self.Phases.keys())
            elif self.multiple:
                choices = sorted(self.Phases.keys())
                phasenum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if phasenum is None: return True
                self.phasenam = [choices[i] for i in phasenum]
                if not self.phasenam: return True
                numselected = len(self.phasenam)
            else:
                choices = sorted(self.Phases.keys())
                phasenum = G2G.ItemSelector(choices,self.G2frame)
                if phasenum is None: return True
                self.phasenam = [choices[phasenum]]
                numselected = len(self.phasenam)
        elif self.currentExportType == 'single':
            if len(self.xtalDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any single crystal data.')
                return True
            elif len(self.xtalDict) == 1:
                self.histnam = list(self.xtalDict.values())
            elif self.multiple:
                choices = sorted(self.xtalDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.xtalDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'powder':
            if len(self.powderDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any powder data.')
                return True
            elif len(self.powderDict) == 1:
                self.histnam = list(self.powderDict.values())
            elif self.multiple:
                choices = sorted(self.powderDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.powderDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'sasd':
            if len(self.sasdDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any small angle data.')
                return True
            elif len(self.sasdDict) == 1:
                self.histnam = self.sasdDict.values()
            elif self.multiple:
                choices = sorted(self.sasdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.sasdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'refd':
            if len(self.refdDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any reflectivity data.')
                return True
            elif len(self.refdDict) == 1:
                self.histnam = self.refdDict.values()
            elif self.multiple:
                choices = sorted(self.refdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.refdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'image':
            if len(self.Histograms) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any images.')
                return True
            elif len(self.Histograms) == 1:
                self.histnam = list(self.Histograms.keys())
            else:
                choices = sorted(self.Histograms.keys())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not hnum: return True
                    self.histnam = [choices[i] for i in hnum]
                else:
                    if hnum is None: return True
                    self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        if self.currentExportType == 'map':
            # search for phases with maps
            mapPhases = []
            choices = []
            for phasenam in sorted(self.Phases):
                phasedict = self.Phases[phasenam] # pointer to current phase info
                if len(phasedict['General']['Map'].get('rho',[])):
                    mapPhases.append(phasenam)
                    if phasedict['General']['Map'].get('Flip'):
                        choices.append('Charge flip map: '+str(phasenam))
                    elif phasedict['General']['Map'].get('MapType'):
                        choices.append(
                            str(phasedict['General']['Map'].get('MapType'))
                            + ' map: ' + str(phasenam))
                    else:
                        choices.append('unknown map: '+str(phasenam))
            # select a map if needed
            if len(mapPhases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any maps.')
                return True
            elif len(mapPhases) == 1:
                self.phasenam = mapPhases
            else:
                phasenum = G2G.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not phasenum: return True
                    self.phasenam = [mapPhases[i] for i in phasenum]
                else:
                    if phasenum is None: return True
                    self.phasenam = [mapPhases[phasenum]]
            numselected = len(self.phasenam)

        # items selected, now set self.dirname and usually self.filename
        if AskFile == 'ask' or (AskFile == 'single' and numselected == 1) or (
            AskFile == 'default' and not self.G2frame.GSASprojectfile
            ):
            filename = self.askSaveFile()
            if not filename: return True
            self.dirname,self.filename = os.path.split(filename)
        elif AskFile == 'dir' or AskFile == 'single' or (
            AskFile == 'default-dir' and not self.G2frame.GSASprojectfile
            ):
            self.dirname = self.askSaveDirectory()
            if not self.dirname: return True
        elif AskFile == 'default-dir' or AskFile == 'default':
            self.dirname,self.filename = os.path.split(
                os.path.splitext(self.G2frame.GSASprojectfile)[0] + self.extension
                )
        else:
            raise Exception('This should not happen!')

    def loadParmDict(self,computeSU=False):
        '''Load the GSAS-II refinable parameters from the tree into a dict (self.parmDict). Update
        refined values to those from the last cycle and set the uncertainties for the
        refined parameters in another dict (self.sigDict).

        Expands the parm & sig dicts to include values derived from constraints.

        This could be made faster for sequential fits as info for each histogram is loaded
        later when iterating over histograms.
        '''
        from . import GSASIIdataGUI as G2gd
        self.G2frame.CheckNotebook()
        self.parmDict = {}
        self.sigDict = {} # dict with s.u. values, currently used only for CIF & Bracket exports
        self.RBsuDict = {} # dict with s.u. values for atoms in a rigid body
        self.constList = [] # constraint entries from data tree
        rigidbodyDict = {}
        covDict = {}
        consDict = {}
        Histograms,Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        if self.G2frame.GPXtree.IsEmpty(): return # nothing to do
        item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.GPXtree.GetItemText(item)
            if name == 'Rigid bodies':
                 rigidbodyDict = self.G2frame.GPXtree.GetItemPyData(item)
            elif name == 'Covariance':
                 covDict = self.G2frame.GPXtree.GetItemPyData(item)
            elif name == 'Constraints':
                 consDict = self.G2frame.GPXtree.GetItemPyData(item)
            item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
        rbVary,rbDict =  G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        self.parmDict.update(rbDict)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave =  \
            G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)
        self.parmDict.update(phaseDict)
        hapVary,hapDict,controlDict =  G2stIO.GetHistogramPhaseData(Phases,Histograms,Print=False,resetRefList=False)
        self.parmDict.update(hapDict)
        histVary,histDict,controlDict =  G2stIO.GetHistogramData(Histograms,Print=False)
        self.parmDict.update(histDict)
        self.parmDict.update(zip(
            covDict.get('varyList',[]),
            covDict.get('variables',[])))
        self.sigDict = dict(zip(
            covDict.get('varyList',[]),
            covDict.get('sig',[])))
        # expand to include constraints: first compile a list of constraints
        self.constList = []
        for item in consDict:
            if item.startswith('_'): continue
            self.constList += consDict[item]
        # now process the constraints
        G2mv.InitVars()
        constrDict,fixedList,ignored = G2mv.ProcessConstraints(self.constList)
        varyList = covDict.get('varyListStart',[])
        if varyList is None and len(constrDict) == 0:
            # no constraints can use varyList
            varyList = covDict.get('varyList')
        elif varyList is None:
            varyList = []
            # # old GPX file from before pre-constraint varyList is saved
            # print (' *** Old refinement: Please use Calculate/Refine to redo  ***')
            # raise Exception(' *** Export aborted ***')
        else:
            varyList = list(varyList)
        # add symmetry-generated constraints
        rigidbodyDict = self.G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,'Rigid bodies'))
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)  # done twice, needed?
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave = \
            G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False) # generates atom symmetry constraints
        msg = G2mv.EvaluateMultipliers(constrDict,phaseDict)
        if msg:
            print('Unable to interpret multiplier(s): '+msg)
            raise Exception(' *** CIF creation aborted ***')
        errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,self.parmDict)
        if errmsg:
            # this really should not happen
            print (' *** ERROR - constraints are internally inconsistent ***')
            print ('Errors: ',errmsg)
            if warnmsg: print ('Warnings'+warnmsg)
            raise Exception(' *** CIF creation aborted ***')
        G2mv.Map2Dict(self.parmDict,varyList)   # changes varyList
        G2mv.Dict2Map(self.parmDict)   # add the constrained values to the parameter dictionary
        # and add their uncertainties into the esd dictionary (sigDict)
        if covDict.get('covMatrix') is not None and computeSU:
            self.sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],covDict['varyList'],noSym=True))
            if 'depSigDict' in self.OverallParms['Covariance']:
                self.sigDict.update(
                    {i:v[1] for i,v in self.OverallParms['Covariance']['depSigDict'].items()})
            # compute the s.u.'s on rigid bodies
            from . import GSASIIstrMath as G2stMth
            self.RBsuDict = G2stMth.computeRBsu(self.parmDict,Phases,rigidbodyDict,
                            covDict['covMatrix'],covDict['varyList'],covDict['sig'])

    def loadTree(self):
        '''Load the contents of the data tree into a set of dicts
        (self.OverallParms, self.Phases and self.Histogram as well as self.powderDict
        & self.xtalDict)

        * The childrenless data tree items are overall parameters/controls for the
          entire project and are placed in self.OverallParms
        * Phase items are placed in self.Phases
        * Data items are placed in self.Histogram. The key for these data items
          begin with a keyword, such as PWDR, IMG, HKLF,... that identifies the data type.
        '''
        from . import GSASIIdataGUI as G2gd
        self.OverallParms = {}
        self.powderDict = {}
        self.sasdDict = {}
        self.refdDict = {}
        self.xtalDict = {}
        self.Phases = {}
        self.Histograms = {}
        self.SeqRefdata = None
        self.SeqRefhist = None
        self.DelayOpen = False
        if self.G2frame.GPXtree.IsEmpty(): return # nothing to do
        histType = None
        if self.currentExportType == 'phase':
            # if exporting phases load them here
            sub = G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,'Phases')
            if not sub:
                print ('no phases found')
                return True
            item, cookie = self.G2frame.GPXtree.GetFirstChild(sub)
            while item:
                phaseName = self.G2frame.GPXtree.GetItemText(item)
                self.Phases[phaseName] =  self.G2frame.GPXtree.GetItemPyData(item)
                item, cookie = self.G2frame.GPXtree.GetNextChild(sub, cookie)
            # Get rigid body info into self.OverallParms
            for key in ('Rigid bodies','Covariance'):
                item = G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,key)
                if item:
                    self.OverallParms[key] = self.G2frame.GPXtree.GetItemPyData(item)
                item, cookie = self.G2frame.GPXtree.GetNextChild(sub, cookie)
            return
        elif self.currentExportType == 'single':
            histType = 'HKLF'
        elif self.currentExportType == 'powder':
            histType = 'PWDR'
        elif self.currentExportType == 'image':
            histType = 'IMG'
        elif self.currentExportType == 'sasd':
            histType = 'SASD'
        elif self.currentExportType == 'refd':
            histType = 'REFD'

        if histType: # Loading just one kind of tree entry
            item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
            while item:
                name = self.G2frame.GPXtree.GetItemText(item)
                if name.startswith(histType):
                    if self.Histograms.get(name): # there is already an item with this name
                        print('Histogram name '+str(name)+' is repeated. Renaming')
                        if name[-1] == '9':
                            name = name[:-1] + '10'
                        elif name[-1] in '012345678':
                            name = name[:-1] + str(int(name[-1])+1)
                        else:
                            name += '-1'
                    self.Histograms[name] = {}
                    # the main info goes into Data, but the 0th
                    # element contains refinement results, carry
                    # that over too now.
                    self.Histograms[name]['Data'] = self.G2frame.GPXtree.GetItemPyData(item)[1]
                    self.Histograms[name][0] = self.G2frame.GPXtree.GetItemPyData(item)[0]
                    item2, cookie2 = self.G2frame.GPXtree.GetFirstChild(item)
                    while item2:
                        child = self.G2frame.GPXtree.GetItemText(item2)
                        self.Histograms[name][child] = self.G2frame.GPXtree.GetItemPyData(item2)
                        item2, cookie2 = self.G2frame.GPXtree.GetNextChild(item, cookie2)
                item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
            # index powder and single crystal histograms by number
            for hist in self.Histograms:
                if hist.startswith("PWDR"):
                    d = self.powderDict
                elif hist.startswith("HKLF"):
                    d = self.xtalDict
                elif hist.startswith("SASD"):
                    d = self.sasdDict
                elif hist.startswith("REFD"):
                    d = self.refdDict
                else:
                    return
                i = self.Histograms[hist].get('hId')
                if i is None and not d.keys():
                    i = 0
                elif i is None or i in d.keys():
                    i = max(d.keys())+1
                d[i] = hist
            return
        # else standard load: using all interlinked phases and histograms
        self.Histograms,self.Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.GPXtree.GetItemText(item)
            if name == 'Restraints':
                self.OverallParms[name] = self.G2frame.GPXtree.GetItemPyData(item)
            else:
                item2, cookie2 = self.G2frame.GPXtree.GetFirstChild(item)
                if not item2:
                    self.OverallParms[name] = self.G2frame.GPXtree.GetItemPyData(item)
            item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
        # index powder and single crystal histograms
        for hist in self.Histograms:
            i = self.Histograms[hist]['hId']
            if hist.startswith("PWDR"):
                self.powderDict[i] = hist
            elif hist.startswith("HKLF"):
                self.xtalDict[i] = hist
            elif hist.startswith("SASD"):
                self.sasdDict[i] = hist
            elif hist.startswith("REFD"):
                self.refdDict[i] = hist

    def dumpTree(self,mode='type'):
        '''Print out information on the data tree dicts loaded in loadTree.
        Used for testing only.
        '''
        if self.SeqRefdata and self.SeqRefhist:
            print('Note that dumpTree does not show sequential results')
        print ('\nOverall')
        if mode == 'type':
            def Show(arg): return type(arg)
        else:
            def Show(arg): return arg
        for key in self.OverallParms:
            print ('  '+key+Show(self.OverallParms[key]))
        print ('Phases')
        for key1 in self.Phases:
            print ('    '+key1+Show(self.Phases[key1]))
        print ('Histogram')
        for key1 in self.Histograms:
            print ('    '+key1+Show(self.Histograms[key1]))
            for key2 in self.Histograms[key1]:
                print ('      '+key2+Show(self.Histograms[key1][key2]))

    def defaultSaveFile(self):
        return os.path.abspath(
            os.path.splitext(self.G2frame.GSASprojectfile
                             )[0]+self.extension)

    def askSaveFile(self):
        '''Ask the user to supply a file name

        :returns: a file name (str) or None if Cancel is pressed

        '''
        from . import GSASIIctrlGUI as G2G
        #pth = G2G.GetExportPath(self.G2frame)
        if self.G2frame.GSASprojectfile:
            defnam = os.path.splitext(
                os.path.split(self.G2frame.GSASprojectfile)[1]
                )[0]+self.extension
        else:
            defnam = 'default' + self.extension
        return G2G.askSaveFile(self.G2frame,defnam,self.extension,self.longFormatName)

    def askSaveDirectory(self):
        '''Ask the user to supply a directory name. Path name is used as the
        starting point for the next export path search.

        :returns: a directory name (str) or None if Cancel is pressed

        TODO: Can this be replaced with G2G.askSaveDirectory?
        '''
        import wx
        from . import GSASIIctrlGUI as G2G
        pth = G2G.GetExportPath(self.G2frame)
        dlg = wx.DirDialog(
            self.G2frame, 'Input directory where file(s) will be written', pth,
            wx.DD_DEFAULT_STYLE)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                self.G2frame.LastExportDir = filename
            else:
                filename = None
        finally:
            dlg.Destroy()
        return filename

    # Tools for file writing.
    def OpenFile(self,fil=None,mode='w',delayOpen=False):
        '''Open the output file

        :param str fil: The name of the file to open. If None (default)
          the name defaults to self.dirname + self.filename.
          If an extension is supplied, it is not overridded,
          but if not, the default extension is used.
        :param str mode:  The mode can 'w' to write a file, or 'a' to append to it. If
          the mode is 'd' (for debug), output is displayed on the console.
        :returns: the file object opened by the routine which is also
          saved as self.fp
        '''
        if mode == 'd': # debug mode
            self.fullpath = '(stdout)'
            self.fp = sys.stdout
            return
        if not fil:
            if not os.path.splitext(self.filename)[1]:
                self.filename += self.extension
            fil = os.path.join(self.dirname,self.filename)
        self.fullpath = os.path.abspath(fil)
        self.DelayOpen = False
        if delayOpen:
            self.DelayOpen = True
            self.fp = None
            return
        self.fp = open(self.fullpath,mode)
        return self.fp

    def openDelayed(self,mode='w'):
        self.DelayOpen = False
        self.fp = open(self.fullpath,mode)
        return self.fp

    def Write(self,line):
        '''write a line of output, attaching a line-end character

        :param str line: the text to be written.
        '''
        if self.fp is None:
            raise Exception('Attempt to Write without use of OpenFile')
        self.fp.write(line+'\n')

    def CloseFile(self,fp=None):
        '''Close a file opened in OpenFile

        :param file fp: the file object to be closed. If None (default)
          file object self.fp is closed.
        '''
        if self.fp is None and self.DelayOpen:
            if GSASIIpath.GetConfigValue('debug'):
                print('Delayed open: close before uncreated file')
            return
        if self.fp is None:
            if GSASIIpath.GetConfigValue('debug'):
                raise Exception('Attempt to CloseFile without use of OpenFile')
            else:
                print('Attempt to CloseFile without use of OpenFile')
                return
        if self.fp == sys.stdout: return # debug mode
        if fp is None:
            fp = self.fp
            self.fp = None
        if fp is not None: fp.close()

    def SetSeqRef(self,data,hist):
        '''Set the exporter to retrieve results from a sequential refinement
        rather than the main tree
        '''
        self.SeqRefdata = data
        self.SeqRefhist = hist
        data_name = data[hist]
        for i,val in zip(data_name['varyList'],data_name['sig']):
            self.sigDict[i] = val
            self.sigDict[striphist(i)] = val
        for i in data_name['parmDict']:
            self.parmDict[striphist(i)] = data_name['parmDict'][i]
            self.parmDict[i] = data_name['parmDict'][i]
            # zero out the dA[xyz] terms, they would only bring confusion
            key = i.split(':')
            if len(key) < 3: continue
            if key[2].startswith('dA'):
                self.parmDict[i] = 0.0
        for i,(val,sig) in data_name.get('depParmDict',{}).items():
            self.parmDict[i] = val
            self.sigDict[i] = sig
        #GSASIIpath.IPyBreak()

    # Tools to pull information out of the data arrays
    def GetCell(self,phasenam,unique=False):
        """Gets the unit cell parameters and their s.u.'s for a selected phase

        :param str phasenam: the name for the selected phase
        :param bool unique: when True, only directly refined parameters
          (a in cubic, a & alpha in rhombohedral cells) are assigned
          positive s.u. values. Used as True for CIF generation.
        :returns: `cellList,cellSig` where each is a 7 element list corresponding
          to a, b, c, alpha, beta, gamma, volume where `cellList` has the
          cell values and `cellSig` has their uncertainties.
        """
        if self.SeqRefdata and self.SeqRefhist:
            return self.GetSeqCell(phasenam,self.SeqRefdata[self.SeqRefhist])
        phasedict = self.Phases[phasenam] # pointer to current phase info
        try:
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            cellSig = G2lat.getCellEsd(pfx,phasedict['General']['SGData'],A,
                self.OverallParms['Covariance'],unique=unique)  # returns 7 vals, includes sigVol
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            return cellList,cellSig
        except KeyError:
            cell = phasedict['General']['Cell'][1:]
            return cell,7*[0]

    def GetSeqCell(self,phasenam,data_name):
        """Gets the unit cell parameters and their s.u.'s for a selected phase
        and histogram in a sequential fit

        :param str phasenam: the name for the selected phase
        :param dict data_name: the sequential refinement parameters for the selected histogram
        :returns: `cellList,cellSig` where each is a 7 element list corresponding
          to a, b, c, alpha, beta, gamma, volume where `cellList` has the
          cell values and `cellSig` has their uncertainties.
        """
        phasedict = self.Phases[phasenam]
        SGdata = phasedict['General']['SGData']
        pId = phasedict['pId']
        RecpCellTerms = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        ESDlookup = {}
        Dlookup = {}
        varied = [striphist(i) for i in data_name['varyList']]
        for item,val in data_name['newCellDict'].items():
            if item in varied:
                ESDlookup[val[0]] = item
                Dlookup[item] = val[0]
        A = RecpCellTerms[:]
        for i in range(6):
            var = str(pId)+'::A'+str(i)
            if var in ESDlookup:
                A[i] = data_name['newCellDict'][ESDlookup[var]][1] # override with refined value
        cellDict = dict(zip([str(pId)+'::A'+str(i) for i in range(6)],A))
        zeroDict = {i:0.0 for i in cellDict}
        A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata,cellDict,zeroDict)
        covData = {
            'varyList': [Dlookup.get(striphist(v),v) for v in data_name['varyList']],
            'covMatrix': data_name['covMatrix']
            }
        return list(G2lat.A2cell(A)) + [G2lat.calc_V(A)], G2lat.getCellEsd(str(pId)+'::',SGdata,A,covData)

    def GetAtoms(self,phasenam):
        """Gets the atoms associated with a phase. Can be used with standard
        or macromolecular phases

        :param str phasenam: the name for the selected phase
        :returns: a list of items for eac atom where each item is a list containing:
          label, typ, mult, xyz, and td, where

          * label and typ are the atom label and the scattering factor type (str)
          * mult is the site multiplicity (int)
          * xyz is contains a list with four pairs of numbers:
            x, y, z and fractional occupancy and
            their standard uncertainty (or a negative value)
          * td is contains a list with either one or six pairs of numbers:
            if one number it is U\\ :sub:`iso` and with six numbers it is
            U\\ :sub:`11`, U\\ :sub:`22`, U\\ :sub:`33`, U\\ :sub:`12`, U\\ :sub:`13` & U\\ :sub:`23`
            paired with their standard uncertainty (or a negative value)
        """
        phasedict = self.Phases[phasenam] # pointer to current phase info
        cx,ct,cs,cia = phasedict['General']['AtomPtrs']
        cfrac = cx+3
        fpfx = str(phasedict['pId'])+'::Afrac:'
        atomslist = []
        for i,at in enumerate(phasedict['Atoms']):
            if phasedict['General']['Type'] == 'macromolecular':
                label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
            else:
                label = at[ct-1]
            fval = self.parmDict.get(fpfx+str(i),at[cfrac])
            fsig = self.sigDict.get(fpfx+str(i),-0.009)
            mult = at[cs+1]
            typ = at[ct]
            xyz = []
            for j,v in enumerate(('x','y','z')):
                val = at[cx+j]
                pfx = str(phasedict['pId']) + '::A' + v + ':' + str(i)
                val = self.parmDict.get(pfx, val)
                dpfx = str(phasedict['pId'])+'::dA'+v+':'+str(i)
                sig = self.sigDict.get(dpfx,-0.000009)
                xyz.append((val,sig))
            xyz.append((fval,fsig))
            td = []
            if at[cia] == 'I':
                pfx = str(phasedict['pId'])+'::AUiso:'+str(i)
                val = self.parmDict.get(pfx,at[cia+1])
                sig = self.sigDict.get(pfx,-0.0009)
                td.append((val,sig))
            else:
                for i,var in enumerate(('AU11','AU22','AU33','AU12','AU13','AU23')):
                    pfx = str(phasedict['pId'])+'::'+var+':'+str(i)
                    val = self.parmDict.get(pfx,at[cia+2+i])
                    sig = self.sigDict.get(pfx,-0.0009)
                    td.append((val,sig))
            atomslist.append((label,typ,mult,xyz,td))
        return atomslist

if __name__ == '__main__':
    for i in (1.23456789e-129,1.23456789e129,1.23456789e-99,1.23456789e99,-1.23456789e-99,-1.23456789e99):
        print (FormatSigFigs(i),i)
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000):
        print (FormatSigFigs(1.23456789e-9*i),1.23456789e-9*i)
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000):
        print (FormatSigFigs(1.23456789e9/i),1.23456789e9/i)

    print (FormatSigFigs(200,10,3))
