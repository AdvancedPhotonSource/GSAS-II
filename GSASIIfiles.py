# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*GSASIIfile: data (non-GUI) I/O routines*
=========================================

Module with miscellaneous routines for input and output from files.

This module should not contain any references to wxPython so that it
can be imported for scriptable use or potentially on clients where
wx is not installed.

Future refactoring: This module and GSASIIIO.py needs some work to
move non-wx routines here. It may will likely make sense to rename the module(s)
at that point.
'''

import os
import sys
import glob
import imp
import inspect
import traceback
import platform

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2957 $")

# N.B. This is duplicated in G2IO
def sfloat(S):
    'Convert a string to float. An empty field or a unconvertable value is treated as zero'
    if S.strip():
        try:
            return float(S)
        except ValueError:
            pass
    return 0.0

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
            names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth']
            v = (v[0],v[2],v[4])
            codes = [0,0,0,0]
            rd.Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],'DisplaceY':[0.,False]})
        else:
            names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','SH/L','Azimuth']
            codes = [0,0,0,0,0,0]
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
                data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],azm)])  #get LX, LY, S+H/L & azimuth
            else:
                data.extend([0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
        else:
            v1 = Iparm['INS  1PRCF1 '].split()
            v = Iparm['INS  1PRCF11'].split()
            data.extend([float(v[0]),float(v[1]),float(v[2])])                  #get GU, GV & GW - always here
            azm = float(Iparm.get('INS  1DETAZM','0.0'))
            v = Iparm['INS  1PRCF12'].split()
            if v1[0] == 3:
                data.extend([float(v[0]),float(v[1]),float(v[2])+float(v[3],azm)])  #get LX, LY, S+H/L & azimuth
            else:
                data.extend([0.0,0.0,0.002,azm])                                      #OK defaults if fxn #3 not 1st in iprm file
        codes.extend([0,0,0,0,0,0,0])
        Iparm1 = makeInstDict(names,data,codes)
        Iparm1['Source'] = [Irads[irad],Irads[irad]]
        Iparm1['Bank'] = [Bank,Bank,0]
        return [Iparm1,{}]
    elif 'T' in DataType:
        names = ['Type','fltPath','2-theta','difC','difA', 'difB','Zero','alpha','beta-0','beta-1',
            'beta-q','sig-0','sig-1','sig-2','sig-q', 'X','Y','Azimuth',]
        codes = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,]
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
                data.extend([0.0,0.0,sfloat(s[1]),sfloat(s[2]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
            elif abs(pfType) in [3,4,5]:
                data.extend([sfloat(s[0]),sfloat(s[1]),sfloat(s[2])]) #alpha, beta-0, beta-1
                if abs(pfType) == 4:
                    data.extend([0.0,0.0,sfloat(s[3]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                else:
                    s = Iparm['INS  1PRCF 2'].split()
                    data.extend([0.0,0.0,sfloat(s[0]),sfloat(s[1]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
            elif abs(pfType) == 2:
                data.extend([sfloat(s[1]),0.0,1./sfloat(s[3])]) #alpha, beta-0, beta-1
                data.extend([0.0,0.0,sfloat(s[1]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
        else:
            s = Iparm['INS  1PRCF1 '].split()
            pfType = int(s[0])
            s = Iparm['INS  1PRCF11'].split()
            if abs(pfType) == 1:
                data.extend([sfloat(s[1]),sfloat(s[2]),sfloat(s[3])]) #alpha, beta-0, beta-1
                s = Iparm['INS  1PRCF12'].split()
                data.extend([0.0,0.0,sfloat(s[1]),sfloat(s[2]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
            elif abs(pfType) in [3,4,5]:
                data.extend([sfloat(s[0]),sfloat(s[1]),sfloat(s[2])]) #alpha, beta-0, beta-1
                if abs(pfType) == 4:
                    data.extend([0.0,0.0,sfloat(s[3]),0.0,0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
                else:
                    s = Iparm['INS  1PRCF12'].split()
                    data.extend([0.0,0.0,sfloat(s[0]),sfloat(s[1]),0.0,0.0,0.0,azm])    #beta-q, sig-0, sig-1, sig-2, sig-q, X, Y
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
            print 'Import_' + errprefix + ': Attribute Error ' + filename
            if traceback:
                traceback.print_exc(file=sys.stdout)
        except Exception as exc:
            print '\nImport_' + errprefix + ': Error importing file ' + filename
            print u'Error message: {}\n'.format(exc)
            if traceback:
                traceback.print_exc(file=sys.stdout)
        finally:
            if fp:
                fp.close()

    return readerlist

def LoadExportRoutines():
    '''Placeholder that will someday retrieve the exporters
    '''
    pass
