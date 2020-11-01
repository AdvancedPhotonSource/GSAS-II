# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIfiles: data (non-GUI) I/O routines*
==========================================

Module with miscellaneous routines for input and output from files.

This module should not contain any references to wxPython so that it
can be imported for scriptable use or potentially on clients where
wx is not installed.

Future refactoring: This module and GSASIIIO.py needs some work to
move non-wx routines here. It may will likely make sense to rename the module(s)
at that point.
'''
from __future__ import division, print_function
import platform
import os
import sys
import glob
import inspect
import numpy as np
import imp

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")

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
'''This defines the level of output from calls to :func:`G2Print`, 
which should  be used in place of print() within this module. 
Settings for this are 'all', 'warn', 'error' or 'none'. Also see:
:func:`G2Print` and :func:`G2SetPrintLevel`.
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
    if G2printLevel is 'none': return
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

def ReadPowderInstprm(instLines, bank, databanks, rd):
    '''Read lines from a GSAS-II (new) instrument parameter file
    similar to G2pwdGUI.OnLoad
    If instprm file has multiple banks each with header #Bank n: ..., this
    finds matching bank no. to load - problem with nonmatches?
    
    Note that this routine performs a similar role to :meth:`GSASIIdataGUI.GSASII.ReadPowderInstprm`,
    but that will call a GUI routine for selection when needed. This routine will raise exceptions
    on errors and will select the first bank when a choice might be appropriate.
    TODO: refactor to combine the two routines. 
    
    :param list instLines: strings from GSAS-II parameter file; can be concatenated with ';'
    :param int  bank: bank number to check when instprm file has '#BANK n:...' strings
         when bank = n then use parameters; otherwise skip that set. Ignored if BANK n:
         not present. NB: this kind of instprm file made by a Save all profile command in Instrument Par     ameters
    :return dict: Inst  instrument parameter dict if OK, or
             str: Error message if failed
    
    (transliterated from GSASIIdataGUI.py:1235 (rev 3008), function of the same name)
     ''' 
    if 'GSAS-II' not in instLines[0]:
        raise ValueError("Not a valid GSAS-II instprm file")

    newItems = []
    newVals = []
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
        rd.powderentry[2] = bank
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
                newItems.append(item)
                try:
                    newVals.append(float(val))
                except ValueError:
                    newVals.append(val)
            il += 1
            continue
        # read multiline values, delimited by ''' or """
        item, val = S.strip().split(':', 1)
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
    if 'Lam1' in newItems:
        rd.Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
            'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
    else:
        rd.Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],'DisplaceY':[0.,False]})
    return [makeInstDict(newItems, newVals, len(newVals)*[False]), {}]

def LoadImportRoutines(prefix, errprefix=None, traceback=False):
    '''Routine to locate GSASII importers matching a prefix string
    '''
    if errprefix is None:
        errprefix = prefix

    readerlist = []
    pathlist = sys.path[:]
    if '.' not in pathlist:
        pathlist.append('.')

    potential_files = []
    for path in pathlist:
        for filename in glob.iglob(os.path.join(path, 'G2'+prefix+'*.py')):
            potential_files.append(filename)

    potential_files = sorted(list(set(potential_files)))
    for filename in potential_files:
        path, rootname = os.path.split(filename)
        pkg = os.path.splitext(rootname)[0]

        try:
            fp = None
            fp, fppath, desc = imp.find_module(pkg, [path])
            pkg = imp.load_module(pkg, fp, fppath, desc)
            for name, value in inspect.getmembers(pkg):
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
            G2Print ('Import_' + errprefix + ': Attribute Error ' + filename)
            if traceback:
                traceback.print_exc(file=sys.stdout)
        except Exception as exc:
            G2Print ('\nImport_' + errprefix + ': Error importing file ' + filename)
            G2Print (u'Error message: {}\n'.format(exc))
            if traceback:
                traceback.print_exc(file=sys.stdout)
        finally:
            if fp:
                fp.close()

    return readerlist

def LoadExportRoutines(parent, traceback=False):
    '''Routine to locate GSASII exporters 
    '''
    exporterlist = []
    pathlist = sys.path
    filelist = []
    for path in pathlist:
        for filename in glob.iglob(os.path.join(path,"G2export*.py")):
                    filelist.append(filename)    
    filelist = sorted(list(set(filelist))) # remove duplicates
    # go through the routines and import them, saving objects that
    # have export routines (method Exporter)
    for filename in filelist:
        path,rootname = os.path.split(filename)
        pkg = os.path.splitext(rootname)[0]
        try:
            fp = None
            fp, fppath,desc = imp.find_module(pkg,[path,])
            pkg = imp.load_module(pkg,fp,fppath,desc)
            for clss in inspect.getmembers(pkg): # find classes defined in package
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
            G2Print ('Export Attribute Error ' + filename)
            if traceback:
                traceback.print_exc(file=sys.stdout)
        except Exception as exc:
            G2Print ('\nExport init: Error importing file ' + filename)
            G2Print (u'Error message: {}\n'.format(exc))
            if traceback:
                traceback.print_exc(file=sys.stdout)
        finally:
            if fp:
                fp.close()
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
        grfile.write('#L r($\AA$)  G($\AA^{-2}$)\n')            
        for r,gr in grnew:
            grfile.write("%15.2F %15.6F\n" % (r,gr))
        grfile.close()
        G2Print (' G(R) saved to: '+grfilename)
        
    if PDFsaves[5]: #RMCProfile files for F(Q) & g(r) overwrites any above
        
        fqfilename = fileroot+'.fq'
        fqdata = PDFControls['F(Q)'][1]
        fqfxn = scintp.interp1d(fqdata[0],fqdata[1],kind='linear')
        fqfile = open(fqfilename,'w')
        qnew = np.arange(fqdata[0][0],fqdata[0][-1],0.005)
        nq = qnew.shape[0]
        fqfile.write('%20d\n'%nq-1)
        fqfile.write(fqfilename+'\n')
        fqnew = zip(qnew,fqfxn(qnew))
        for q,fq in fqnew[1:]:
            fqfile.write("%15.6g %15.6g\n" % (q,fq))
        fqfile.close()
        G2Print (' F(Q) saved to: '+fqfilename)
        
        grfilename = fileroot+'.gr'
        grdata = PDFControls['g(r)'][1]
        grfxn = scintp.interp1d(grdata[0],grdata[1],kind='linear')
        grfile = open(grfilename,'w')
        rnew = np.arange(grdata[0][0],grdata[0][-1],0.010)
        nr = rnew.shape[0]
        grfile.write('%20d\n'%nr-1)
        grfile.write(grfilename+'\n')
        grnew = zip(rnew,grfxn(rnew))
        for r,gr in grnew[1:]:
            grfile.write("%15.6g %15.6g\n" % (r,gr))
        grfile.close()
        G2Print (' G(R) saved to: '+grfilename)

