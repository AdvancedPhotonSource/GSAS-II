# -*- coding: utf-8 -*-
'''
This module should not contain any references to wxPython so that it
can be imported for scriptable use or potentially on clients where
wx is not installed.

Future refactoring: Module :mod:`GSASIIIO` needs some work to
move non-wx routines to here and wx routines to a GSASII*GUI.py file. 
It will likely make sense to rename the GSASIIIO module after that is done.
'''
from __future__ import division, print_function
import platform
import os
import sys
import glob
import inspect
import numpy as np

import GSASIIpath

if not sys.platform.startswith('win'):
    try:
        from dmp import dump2tmp,undumptmp
    except:
        print('Note: Import of dmp skipped')

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

def LoadImportRoutines(prefix, errprefix=None, traceback=False):
    '''Routine to locate GSASII importers matching a prefix string.

    Warns if more than one file with the same name is in the path 
    or if a file is found that is not in the main directory tree. 
    '''
    if errprefix is None:
        errprefix = prefix

    readerlist = []
    import_files = {}
    if '.' not in sys.path: sys.path.append('.')
    for path in sys.path:
        for filename in glob.iglob(os.path.join(path, 'G2'+prefix+'*.py')):
            pkg = os.path.splitext(os.path.split(filename)[1])[0]
            if pkg in import_files:
                G2Print('Warning: importer {} overrides {}'.format(import_files[pkg],os.path.abspath(filename)))
            elif not filename.startswith(GSASIIpath.path2GSAS2):
                G2Print('Note: found importer in non-standard location:'+
                            f'\n\t{os.path.abspath(filename)}')
                import_files[pkg] = filename
            else:
                import_files[pkg] = filename

    for pkg in sorted(import_files.keys()):
        try:
            exec('import '+pkg)
            #print(eval(pkg+'.__file__'))
            for name, value in inspect.getmembers(eval(pkg)):
                if name.startswith('_'):
                    continue
                if inspect.isclass(value):
                    for method in 'Reader', 'ExtensionValidator', 'ContentsValidator':
                        if not hasattr(value, method):
                            break
                        if not callable(getattr(value, method)):
                            break
                    else:
                        reader = value()
                        if reader.UseReader:
                            readerlist.append(reader)
        except AttributeError:
            G2Print ('Import_' + errprefix + ': Attribute Error ' + import_files[pkg])
            if traceback:
                traceback.print_exc(file=sys.stdout)
        except Exception as exc:
            G2Print ('\nImport_' + errprefix + ': Error importing file ' + import_files[pkg])
            G2Print (u'Error message: {}\n'.format(exc))
            if traceback:
                traceback.print_exc(file=sys.stdout)
    return readerlist

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
    if pkg: condaRequestList.update(pkg)

def LoadExportRoutines(parent, traceback=False):
    '''Routine to locate GSASII exporters. Warns if more than one file
    with the same name is in the path or if a file is found that is not 
    in the main directory tree. 
    '''
    exporterlist = []
    export_files = {}
    if '.' not in sys.path: sys.path.append('.')
    for path in sys.path:
        for filename in glob.iglob(os.path.join(path,"G2export*.py")):
            pkg = os.path.splitext(os.path.split(filename)[1])[0]
            if pkg in export_files:
                G2Print('Warning: exporter {} overrides {}'.format(export_files[pkg],os.path.abspath(filename)))
            elif not filename.startswith(GSASIIpath.path2GSAS2):
                G2Print('Note, found non-standard exporter: {}'.format(os.path.abspath(filename)))
                export_files[pkg] = filename
            else:
                export_files[pkg] = filename
    # go through the routines and import them, saving objects that
    # have export routines (method Exporter)
    for pkg in sorted(export_files.keys()):
        try:
            exec('import '+pkg)
            for clss in inspect.getmembers(eval(pkg)): # find classes defined in package
                if clss[0].startswith('_'): continue
                if not inspect.isclass(clss[1]): continue
                # check if we have the required methods
                if not hasattr(clss[1],'Exporter'): continue
                if not callable(getattr(clss[1],'Exporter')): continue
                if parent is None:
                    if not hasattr(clss[1],'Writer'): continue
                else:
                    if not hasattr(clss[1],'loadParmDict'): continue
                    if not callable(getattr(clss[1],'loadParmDict')): continue
                try:
                    exporter = clss[1](parent) # create an export instance
                except AttributeError:
                    pass
                except Exception as exc:
                    G2Print ('\nExport init: Error substantiating class ' + clss[0])
                    G2Print (u'Error message: {}\n'.format(exc))
                    if traceback:
                        traceback.print_exc(file=sys.stdout)
                    continue
                exporterlist.append(exporter)
        except AttributeError:
            G2Print ('Export Attribute Error ' + export_files[pkg])
            if traceback:
                traceback.print_exc(file=sys.stdout)
        except Exception as exc:
            G2Print ('\nExport init: Error importing file ' + export_files[pkg])
            G2Print (u'Error message: {}\n'.format(exc))
            if traceback:
                traceback.print_exc(file=sys.stdout)
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

def RereadImageData(ImageReaderlist,imagefile,ImageTag=None,FormatName=''):
    '''Read a single image with an image importer. This is called to 
    reread an image after it has already been imported, so it is not 
    necessary to reload metadata.

    Based on :func:`GetImageData.GetImageData` which this can replace
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
        if key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds']:
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
        import GSASIImath as G2mth
        import GSASIIlattice as G2lat       
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
      the GSASII.py file in the same directory as this file will
      be used. 
    :param str pythonapp: the Python interpreter to be used. 
      Defaults to sys.executable which is usually what is wanted.
    :param str terminal: a name for a preferred terminal emulator
    '''
    import subprocess
    if g2script is None:
        g2script = os.path.join(os.path.dirname(__file__),'GSASII.py')
    
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
        
if __name__ == '__main__':
    for i in (1.23456789e-129,1.23456789e129,1.23456789e-99,1.23456789e99,-1.23456789e-99,-1.23456789e99):
        print (FormatSigFigs(i),i)
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000):
        print (FormatSigFigs(1.23456789e-9*i),1.23456789e-9*i)
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000):
        print (FormatSigFigs(1.23456789e9/i),1.23456789e9/i)

    print (FormatSigFigs(200,10,3))
