#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
Classes and routines defined in :mod:`~GSASII.GSASIIscriptable` follow.
A script will create one or more :class:`G2Project` objects by reading
a GSAS-II project (.gpx) file or creating a new one and will then
perform actions such as adding a histogram (method :meth:`G2Project.add_powder_histogram`),
adding a phase (method :meth:`G2Project.add_phase`),
or setting parameters and performing a refinement
(method :meth:`G2Project.do_refinements`).

To change settings within histograms, images and phases, one usually needs to use
methods inside :class:`G2PwdrData`, :class:`G2Image` or :class:`G2Phase`.
"""
# Note that documentation for GSASIIscriptable.py has been moved
# to file docs/source/GSASIIscriptable.rst

#============================================================================
# Notes for adding a new object type
# 1) add a new object class (e.g. G2PDF)
# 2) add the wrapper into G2Project (e.g. _pdfs, pdf, pdfs)
# 3) add a new method to add the object into a project (G2Project.add_PDF)
# Document: (in ../docs/source/GSASIIscriptable.rst)
# 4) add to documentation in section :class:`G2Project`
# 5) add a new documentation section for the new class
#============================================================================

from __future__ import division, print_function
import argparse
import os.path as ospath
import sys
import platform
import pickle
import copy
import os
import random as ran

import numpy as np
import numpy.ma as ma

from . import GSASIIpath
GSASIIpath.SetBinaryPath(True)  # for now, this is needed before some of these modules can be imported
from . import GSASIIobj as G2obj
from . import GSASIIpwd as G2pwd
from . import GSASIIstrMain as G2strMain
from . import GSASIIstrIO as G2stIO
from . import GSASIIspc as G2spc
from . import GSASIIElem as G2elem
from . import GSASIIfiles as G2fil
from . import GSASIIimage as G2img
from . import GSASIIlattice as G2lat
from . import GSASIImapvars as G2mv
from . import GSASIImath as G2mth

# Delay imports loading to not slow down small scripts that don't need them
Readers = {'Pwdr': [], 'Phase': [], 'Image': [], 'importErrpkgs': []}
'''Readers by reader type'''
exportersByExtension = {}
'''Specifies the list of extensions that are supported for Powder data export'''
npsind = lambda x: np.sin(x*np.pi/180.)

def SetPrintLevel(level):
    '''Set the level of output from calls to :func:`GSASIIfiles.G2Print`,
    which should be used in place of print() where possible. This is a
    wrapper for :func:`GSASIIfiles.G2SetPrintLevel` so that this routine is
    documented here.

    :param str level: a string used to set the print level, which may be
      'all', 'warn', 'error' or 'none'.
      Note that capitalization and extra letters in level are ignored, so
      'Warn', 'warnings', etc. will all set the mode to 'warn'
    '''
    G2fil.G2SetPrintLevel(level)
    global printLevel
    for mode in  'all', 'warn', 'error', 'none':
        if mode in level.lower():
            printLevel = mode
            return
def SetDebugMode(mode):
    '''Set the debug configuration mode on (mode=True) or off (mode=False).
    This will provide some additional output that may help with 
    tracking down problems in the code.
    '''
    GSASIIpath.SetConfigValue({'debug':bool(mode)})

def installScriptingShortcut():
    '''Creates a file named G2script in the current Python site-packages directory.
    This is equivalent to the "Install GSASIIscriptable shortcut" command in the GUI's
    File menu. Once this is done, a shortcut for calling GSASIIscriptable is created,
    where the command:

    >>> import G2script as G2sc

    will provide access to GSASIIscriptable without changing the sys.path; also see
    :ref:`ScriptingShortcut`.

    Note that this only affects the current Python installation. If more than one
    Python installation will be used with GSAS-II (for example because different
    conda environments are used), this command should be called from within each
    Python environment.

    If more than one GSAS-II installation will be used with a Python installation,
    this shortcut can only be used with one of them.
    '''
    f = GSASIIpath.makeScriptShortcut()
    if f:
        G2fil.G2Print(f'success creating {f}')
    else:
        raise G2ScriptException('error creating G2script')

def ShowVersions():
    '''Show the versions all of required Python packages, etc.
    '''
    out = ''
    pkgList = [('Python',None), ('numpy',np)]
    try:
        import scipy
        pkgList.append(('scipy',scipy))
    except:
        pass
    try:
        import IPython
        pkgList.append(('IPython',IPython))
    except:
        pass
    for s,m in pkgList:
        msg = ''
        if s == 'Python':
            pkgver = platform.python_version()
            #prefix = ''
            msg = f"from {format(sys.executable)}"
        else:
            pkgver = m.__version__
            #prefix = 'Package '
        out += f"  {s:12s}{pkgver}:  {msg}\n"
    out += GSASIIpath.getG2VersionInfo()
    out += "\n\n"
    try:
        out += f"GSAS-II location: {GSASIIpath.path2GSAS2}\n"
    except:
        out += "GSAS-II location: not set\n"
    try:
        if GSASIIpath.binaryPath is None:
            from . import pyspg
            loc = os.path.dirname(pyspg.__file__)
        else:
            loc = GSASIIpath.binaryPath
        try:
            f = os.path.join(loc,'GSASIIversion.txt')
            with open(f,'r') as fp:
                version = fp.readline().strip()
                vnum = fp.readline().strip()
            dated = f'{vnum}, {version}'
        except:
            dated = 'undated'
        out += f"Binary location:  {loc} ({dated})\n"
    except:
        out += "Binary location:  not found\n"
    return out

def LoadG2fil():
    '''Setup GSAS-II importers.
    Delay importing this module when possible, it is slow.
    Multiple calls are not. Only the first does anything.
    '''
    if len(Readers['Pwdr']) > 0: return
    # initialize imports
    before = list(G2fil.condaRequestList.keys())
    Readers['Pwdr'] = G2fil.LoadImportRoutines("pwd", "Powder_Data")
    Readers['Phase'] = G2fil.LoadImportRoutines("phase", "Phase")
    Readers['Image'] = G2fil.LoadImportRoutines("img", "Image")
    Readers['HKLF'] = G2fil.LoadImportRoutines('sfact','Struct_Factor')
    # save list of importers that could not be loaded
    Readers['importErrpkgs'] = [i for i in G2fil.condaRequestList.keys()
                                if i not in before]

    # initialize exports
    for obj in G2fil.LoadExportRoutines(None):
        try:
            obj.Writer
        except AttributeError:
            continue
        for typ in obj.exporttype:
            if typ not in exportersByExtension:
                exportersByExtension[typ] = {obj.extension:obj}
            elif obj.extension in exportersByExtension[typ]:
                if type(exportersByExtension[typ][obj.extension]) is list:
                    exportersByExtension[typ][obj.extension].append(obj)
                else:
                    exportersByExtension[typ][obj.extension] = [
                        exportersByExtension[typ][obj.extension],
                        obj]
            else:
                exportersByExtension[typ][obj.extension] = obj

def LoadDictFromProjFile(ProjFile):
    '''Read a GSAS-II project file and load items to dictionary

    :param str ProjFile: GSAS-II project (name.gpx) full file name
    :returns: Project,nameList, where

      * Project (dict) is a representation of gpx file following the GSAS-II tree structure
        for each item: key = tree name (e.g. 'Controls','Restraints',etc.), data is dict
        data dict = {'data':item data whch may be list, dict or None,'subitems':subdata (if any)}
      * nameList (list) has names of main tree entries & subentries used to reconstruct project file

    Example for fap.gpx::

      Project = {                 #NB:dict order is not tree order
        'Phases':{'data':None,'fap':{phase dict}},
        'PWDR FAP.XRA Bank 1':{'data':[histogram data list],'Comments':comments,'Limits':limits, etc},
        'Rigid bodies':{'data': {rigid body dict}},
        'Covariance':{'data':{covariance data dict}},
        'Controls':{'data':{controls data dict}},
        'Notebook':{'data':[notebook list]},
        'Restraints':{'data':{restraint data dict}},
        'Constraints':{'data':{constraint data dict}}]
        }
      nameList = [                #NB: reproduces tree order
        ['Notebook',],
        ['Controls',],
        ['Covariance',],
        ['Constraints',],
        ['Restraints',],
        ['Rigid bodies',],
        ['PWDR FAP.XRA Bank 1',
             'Comments',
             'Limits',
             'Background',
             'Instrument Parameters',
             'Sample Parameters',
             'Peak List',
             'Index Peak List',
             'Unit Cells List',
             'Reflection Lists'],
        ['Phases', 'fap']
        ]
    '''
    # Let IOError be thrown if file does not exist
    if not ospath.exists(ProjFile):
        G2fil.G2Print ('\n*** Error attempt to open project file that does not exist: \n    {}'.
                   format(ProjFile))
        raise IOError('GPX file {} does not exist'.format(ProjFile))
    errmsg = ''
    try:
        Project, nameList = G2stIO.GetFullGPX(ProjFile)
    except Exception as msg:
        errmsg = msg
    if errmsg: raise IOError(errmsg)
    return Project,nameList

def SaveDictToProjFile(Project,nameList,ProjFile):
    '''Save a GSAS-II project file from dictionary/nameList created by LoadDictFromProjFile

    :param dict Project: representation of gpx file following the GSAS-II
        tree structure as described for LoadDictFromProjFile
    :param list nameList: names of main tree entries & subentries used to reconstruct project file
    :param str ProjFile: full file name for output project.gpx file (including extension)
    '''
    file = open(ProjFile,'wb')
    try:
        for name in nameList:
            data = []
            item = Project[name[0]]
            data.append([name[0],item['data']])
            for item2 in name[1:]:
                data.append([item2,item[item2]])
            pickle.dump(data,file,1)
    finally:
        file.close()
    G2fil.G2Print('gpx file saved as %s'%ProjFile)

def PreSetup(data):
    '''Create part of an initial (empty) phase dictionary

    from GSASIIphsGUI.py, near end of UpdatePhaseData

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    '''
    if 'RBModels' not in data:
        data['RBModels'] = {}
    if 'MCSA' not in data:
        data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
    if 'dict' in str(type(data['MCSA']['Results'])):
        data['MCSA']['Results'] = []
    if 'Modulated' not in data['General']:
        data['General']['Modulated'] = False
#    if 'modulated' in data['General']['Type']:
#        data['General']['Modulated'] = True
#        data['General']['Type'] = 'nuclear'


def SetupGeneral(data, dirname):
    '''Initialize phase data.
    '''
    errmsg = ''
    try:
        G2elem.SetupGeneral(data,dirname)
    except ValueError as msg:
        errmsg = msg
    if errmsg: raise G2ScriptException(errmsg)

def make_empty_project(author=None, filename=None):
    """Creates an dictionary in the style of GSASIIscriptable, for an empty
    project.

    If no author name or filename is supplied, 'no name' and
    <current dir>/test_output.gpx are used , respectively.

    Returns: project dictionary, name list

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    if not filename:
        filename = 'test_output.gpx'
    filename = os.path.abspath(filename)
    LoadG2fil()

    controls_data = dict(G2obj.DefaultControls)
    if author:
        controls_data['Author'] = author

    output = {'Constraints': {'data': {'Hist': [], 'HAP': [], 'Phase': [],
                                       'Global': [],
                                       '_seqmode':'auto-wildcard',
                                       '_seqhist':0}},
              'Controls': {'data': controls_data},
              'Covariance': {'data': {}},
              'Notebook': {'data': ['']},
              'Restraints': {'data': {}},
              'Rigid bodies': {'data':
                    {'Vector':{'AtInfo':{}},'Residue':{'AtInfo':{}},
                         "Spin": {},
                         'RBIds':{'Vector':[], 'Residue':[],'Spin':[]}}}
             }
    names = [['Notebook'], ['Controls'], ['Covariance'],
             ['Constraints'], ['Restraints'], ['Rigid bodies']]

    return output, names

def GenerateReflections(spcGrp,cell,Qmax=None,dmin=None,TTmax=None,wave=None):
    """Generates the crystallographically unique powder diffraction reflections
    for a lattice and space group (see :func:`GSASIIlattice.GenHLaue`).

    :param str spcGrp: A GSAS-II formatted space group (with spaces between
       axial fields, e.g. 'P 21 21 21' or 'P 42/m m c'). Note that non-standard
       space groups, such as 'P 21/n' or 'F -1' are allowed (see
       :func:`GSASIIspc.SpcGroup`).
    :param list cell: A list/tuple with six unit cell constants,
      (a, b, c, alpha, beta, gamma) with values in Angstroms/degrees.
      Note that the cell constants are not checked for consistency
      with the space group.
    :param float Qmax: Reflections up to this Q value are computed
       (do not use with dmin or TTmax)
    :param float dmin: Reflections with d-space above this value are computed
       (do not use with Qmax or TTmax)
    :param float TTmax: Reflections up to this 2-theta value are computed
       (do not use with dmin or Qmax, use of wave is required.)
    :param float wave: wavelength in Angstroms for use with TTmax (ignored
       otherwise.)
    :returns: a list of reflections, where each reflection contains four items:
       h, k, l, d, where d is the d-space (Angstroms)

    Example:

    >>> refs = G2sc.GenerateReflections('P 1',
    ...                     (5.,6.,7.,90.,90.,90),
    ...                     TTmax=20,wave=1)
    >>> for r in refs: print(r)
    ...
    [0, 0, 1, 7.0]
    [0, 1, 0, 6.0]
    [1, 0, 0, 5.0]
    [0, 1, 1, 4.55553961419178]
    [0, 1, -1, 4.55553961419178]
    [1, 0, 1, 4.068667356033675]
    [1, 0, -1, 4.068667356033674]
    [1, 1, 0, 3.8411063979868794]
    [1, -1, 0, 3.8411063979868794]
    """

    if len(cell) != 6:
        raise G2ScriptException("GenerateReflections: Invalid unit cell:" + str(cell))
    opts = (Qmax is not None) + (dmin is not None) + (TTmax is not None)
    if Qmax:
        dmin = 2 * np.pi / Qmax
        #print('Q,d',Qmax,dmin)
    elif TTmax and wave is None:
        raise G2ScriptException("GenerateReflections: specify a wavelength with TTmax")
    elif TTmax:
        dmin = wave / (2.0 * np.sin(np.pi*TTmax/360.))
        #print('2theta,d',TTmax,dmin)
    if opts != 1:
        raise G2ScriptException("GenerateReflections: specify one Qmax, dmin or TTmax")
    err,SGData = G2spc.SpcGroup(spcGrp)
    if err != 0:
        print('GenerateReflections space group error:',G2spc.SGErrors(err))
        raise G2ScriptException("GenerateReflections: Invalid space group: " + str(spcGrp))
    A = G2lat.cell2A(cell)
    return G2lat.GenHLaue(dmin,SGData,A)

class G2ImportException(Exception):
    pass

class G2ScriptException(Exception):
    pass

def downloadFile(URL,download_loc=None):
    '''Download the URL 
    '''
        
    import requests
    fname = os.path.split(URL)[1]
    if download_loc is None:
       import tempfile
       download_loc = tempfile.gettempdir()
    elif os.path.isdir(download_loc) and os.path.exists(download_loc):
        pass
    elif os.path.exists(os.path.dirname(download_loc)):
        download_loc,fname = os.path.split(download_loc)
        pass
    else:
        raise G2ScriptException(f"Import error: Cannot download to {download_loc}")
    G2fil.G2Print(f'Preparing to download {URL}')
    response = requests.get(URL)
    filename = os.path.join(download_loc,fname)
    with open(filename,'wb') as fp:
        fp.write(response.content)
    G2fil.G2Print(f'File {filename} written')
    return filename

def import_generic(filename, readerlist, fmthint=None, bank=None,
                       URL=False, download_loc=None):
    """Attempt to import a filename, using a list of reader objects.

    Returns the first reader object which worked."""
    if URL is True:
        filename = downloadFile(filename,download_loc)
    # Translated from OnImportGeneric method in GSASII.py
    primaryReaders, secondaryReaders = [], []
    hintcount = 0
    for reader in readerlist:
        if fmthint is not None and fmthint not in reader.formatName: continue
        hintcount += 1
        flag = reader.ExtensionValidator(filename)
        if flag is None:
            secondaryReaders.append(reader)
        elif flag:
            primaryReaders.append(reader)
    if not secondaryReaders and not primaryReaders:
        # common reason for read error -- package needed?
        l = []
        for i in Readers['importErrpkgs']:
            for j in G2fil.condaRequestList[i]:
                if j not in l: l.append(j)
        print(f'\nReading of file {filename!r} failed.')
        if fmthint is not None and hintcount == 0:
            print(f'\nNo readers matched hint {fmthint!r}\n')
        if Readers['importErrpkgs']: 
            print('Not all importers are available due to uninstalled Python packages.')
            print('The following importer(s) are not available:\n'+
                  f'\t{", ".join(Readers["importErrpkgs"])}')
            print('because the following optional Python package(s) are not installed:\n'+
                  f'\t{", ".join(l)}\n')
        print('Available importers:')
        for reader in readerlist:
            print(f'\t{reader.longFormatName}')
        raise G2ImportException(f"Could not read file: {filename}")

    with open(filename, 'r'):
        rd_list = []

        for rd in primaryReaders + secondaryReaders:
            # Initialize reader
            rd.selections = []
            if bank is None:
                rd.selections = []
            else:
                try:
                    rd.selections = [i-1 for i in bank]
                except TypeError:
                    rd.selections = [bank-1]
            rd.dnames = []
            rd.ReInitialize()
            # Rewind file
            rd.errors = ""
            if not rd.ContentsValidator(filename):
                # Report error
                G2fil.G2Print("Warning: File {} has a validation error, continuing".format(filename))
            #if len(rd.selections) > 1:
            #    raise G2ImportException("File {} has {} banks. Specify which bank to read with databank param."
            #                    .format(filename,len(rd.selections)))

            block = 0
            rdbuffer = {}
            repeat = True
            while repeat:
                repeat = False
                block += 1
                rd.objname = os.path.basename(filename)
                try:
                    flag = rd.Reader(filename,buffer=rdbuffer, blocknum=block)
                except Exception as msg:
                    if GSASIIpath.GetConfigValue('debug'):
                        print('Reader exception',msg)
                    flag = False
                if flag:
                    # Omitting image loading special cases
                    rd.readfilename = filename
                    rd_list.append(copy.deepcopy(rd))
                    repeat = rd.repeat
                else:
                    G2fil.G2Print("Warning: {} Reader failed to read {}".format(rd.formatName,filename))
            if rd_list:
                if rd.warnings:
                    G2fil.G2Print("Read warning by", rd.formatName, "reader:",
                          rd.warnings)
                elif bank is None:
                    G2fil.G2Print("{} read by Reader {}"
                              .format(filename,rd.formatName))
                else:
                    G2fil.G2Print("{} block # {} read by Reader {}"
                              .format(filename,bank,rd.formatName))
                return rd_list
    raise G2ImportException(f"No reader could read file: {filename}")


def load_iprms(instfile, reader, bank=None):
    """Loads instrument parameters from a file, and edits the
    given reader.

    Returns a 2-tuple of (Iparm1, Iparm2) parameters
    """
    LoadG2fil()
    ext = os.path.splitext(instfile)[1]

    if ext.lower() == '.instprm':
        # New GSAS File, load appropriate bank
        with open(instfile) as f:
            lines = f.readlines()
        if bank is None:
            bank = reader.powderentry[2]
        nbank,iparms = G2fil.ReadInstprm(lines, bank, reader.Sample)

        reader.instfile = instfile
        reader.instmsg = '{} (G2 fmt) bank {}'.format(instfile,bank)
        return iparms
    elif ext.lower() not in ('.prm', '.inst', '.ins'):
        raise ValueError('Expected .prm file, found: ', instfile)

    # It's an old GSAS file, load appropriately
    Iparm = {}
    with open(instfile, 'r') as fp:
        for line in fp:
            if '#' in line:
                continue
            Iparm[line[:12]] = line[12:-1]
    ibanks = int(Iparm.get('INS   BANK  ', '1').strip())
    if bank is not None:
        # pull out requested bank # bank from the data, and change the bank to 1
        Iparm,IparmC = {},Iparm
        for key in IparmC:
            if 'INS' not in key[:3]: continue   #skip around rubbish lines in some old iparm
            if key[4:6] == "  ":
                Iparm[key] = IparmC[key]
            elif int(key[4:6].strip()) == bank:
                Iparm[key[:4]+' 1'+key[6:]] = IparmC[key]
        reader.instbank = bank
    elif ibanks == 1:
        reader.instbank = 1
    else:
        raise G2ImportException(f"Instrument parameter file has {ibanks} banks, select one with instbank param.")
    reader.powderentry[2] = 1
    reader.instfile = instfile
    reader.instmsg = '{} bank {}'.format(instfile,reader.instbank)
    return G2fil.SetPowderInstParms(Iparm, reader)

def load_pwd_from_reader(reader, instprm, existingnames=[],bank=None):
    """Loads powder data from a reader object, and assembles it into a GSASII data tree.

    :returns: (name, tree) - 2-tuple of the histogram name (str), and data

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    HistName = 'PWDR ' + G2obj.StripUnicode(reader.idstring, '_')
    HistName = G2obj.MakeUniqueLabel(HistName, existingnames)

    try:
        Iparm1, Iparm2 = instprm
    except ValueError:
        Iparm1, Iparm2 = load_iprms(instprm, reader, bank=bank)
        G2fil.G2Print('Instrument parameters read:',reader.instmsg)
    except TypeError:  # instprm is None, get iparms from reader
        Iparm1, Iparm2 = reader.pwdparms['Instrument Parameters']

    if 'T' in Iparm1['Type'][0]:
        if not reader.clockWd and reader.GSAS:
            reader.powderdata[0] *= 100.        #put back the CW centideg correction
        cw = np.diff(reader.powderdata[0])
        reader.powderdata[0] = reader.powderdata[0][:-1]+cw/2.
        if reader.GSAS:     #NB: old GSAS wanted intensities*CW even if normalized!
            npts = min(len(reader.powderdata[0]),len(reader.powderdata[1]),len(cw))
            reader.powderdata[1] = reader.powderdata[1][:npts]/cw[:npts]
            reader.powderdata[2] = reader.powderdata[2][:npts]*cw[:npts]**2  #1/var=w at this point
        else:       #NB: from topas/fullprof type files
            reader.powderdata[1] = reader.powderdata[1][:-1]
            reader.powderdata[2] = reader.powderdata[2][:-1]
        if 'Itype' in Iparm2:
            Ibeg = np.searchsorted(reader.powderdata[0],Iparm2['Tminmax'][0])
            Ifin = np.searchsorted(reader.powderdata[0],Iparm2['Tminmax'][1])
            reader.powderdata[0] = reader.powderdata[0][Ibeg:Ifin]
            YI,WYI = G2pwd.calcIncident(Iparm2,reader.powderdata[0])
            reader.powderdata[1] = reader.powderdata[1][Ibeg:Ifin]/YI
            var = 1./reader.powderdata[2][Ibeg:Ifin]
            var += WYI*reader.powderdata[1]**2
            var /= YI**2
            reader.powderdata[2] = 1./var
        reader.powderdata[1] = np.where(np.isinf(reader.powderdata[1]),0.,reader.powderdata[1])
        reader.powderdata[3] = np.zeros_like(reader.powderdata[0])
        reader.powderdata[4] = np.zeros_like(reader.powderdata[0])
        reader.powderdata[5] = np.zeros_like(reader.powderdata[0])

    Ymin = np.min(reader.powderdata[1])
    Ymax = np.max(reader.powderdata[1])
    valuesdict = {'wtFactor': 1.0,
                  'Dummy': False,
                  'ranId': ran.randint(0, sys.maxsize),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'Yminmax': [Ymin, Ymax]}
    # apply user-supplied corrections to powder data
    if 'CorrectionCode' in Iparm1:
        G2fil.G2Print('Applying corrections from instprm file')
        corr = Iparm1['CorrectionCode'][0]
        try:
            exec(corr)
            G2fil.G2Print('done')
        except Exception as err:
            print('error: {}'.format(err))
            print('with commands -------------------')
            print(corr)
            print('---------------------------------')
        finally:
            del Iparm1['CorrectionCode']
    reader.Sample['ranId'] = valuesdict['ranId']

    # Ending keys:
    # ['Reflection Lists',
    #  'Limits',
    #  'data',
    #  'Index Peak List',
    #  'Comments',
    #  'Unit Cells List',
    #  'Sample Parameters',
    #  'Peak List',
    #  'Background',
    #  'Instrument Parameters']
    Tmin = np.min(reader.powderdata[0])
    Tmax = np.max(reader.powderdata[0])
    Tmin1 = Tmin
    if 'NT' in Iparm1['Type'][0] and G2lat.Pos2dsp(Iparm1,Tmin) < 0.4:
        Tmin1 = G2lat.Dsp2pos(Iparm1,0.4)

    default_background = [['chebyschev-1', False, 3, 1.0, 0.0, 0.0],
        {'nDebye': 0, 'debyeTerms': [], 'nPeaks': 0,
        'peaksList': [],'background PWDR':['',1.0,False]}]

    output_dict = {'Reflection Lists': {},
                   'Limits': reader.pwdparms.get('Limits', [(Tmin, Tmax), [Tmin1, Tmax]]),
                   'data': [valuesdict, reader.powderdata, HistName],
                   'Index Peak List': [[], []],
                   'Comments': reader.comments,
                   'Unit Cells List': [],
                   'Sample Parameters': reader.Sample,
                   'Peak List': {'peaks': [], 'sigDict': {}},
                   'Background': reader.pwdparms.get('Background', default_background),
                   'Instrument Parameters': [Iparm1, Iparm2],
                   }

    names = ['Comments',
             'Limits',
             'Background',
             'Instrument Parameters',
             'Sample Parameters',
             'Peak List',
             'Index Peak List',
             'Unit Cells List',
             'Reflection Lists']

    # TODO controls?? GSASII.py:1664-7

    return HistName, [HistName] + names, output_dict

def _deep_copy_into(from_, into):
    """Helper function for reloading .gpx file. See G2Project.reload()

    :author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    if isinstance(from_, dict) and isinstance(into, dict):
        combined_keys = set(from_.keys()).union(into.keys())
        for key in combined_keys:
            if key in from_ and key in into:
                both_dicts = (isinstance(from_[key], dict)
                              and isinstance(into[key], dict))
                both_lists = (isinstance(from_[key], list)
                              and isinstance(into[key], list))
                if both_dicts or both_lists:
                    _deep_copy_into(from_[key], into[key])
                else:
                    into[key] = from_[key]
            elif key in from_:
                into[key] = from_[key]
            else:  # key in into
                del into[key]
    elif isinstance(from_, list) and isinstance(into, list):
        if len(from_) == len(into):
            for i in range(len(from_)):
                both_dicts = (isinstance(from_[i], dict)
                              and isinstance(into[i], dict))
                both_lists = (isinstance(from_[i], list)
                              and isinstance(into[i], list))
                if both_dicts or both_lists:
                    _deep_copy_into(from_[i], into[i])
                else:
                    into[i] = from_[i]
        else:
            into[:] = from_

def _getCorrImage(ImageReaderlist,proj,imageRef):
    '''Gets image & applies dark, background & flat background corrections.
    based on :func:`GSASIIimgGUI.GetImageZ`. Expected to be for internal
    use only.

    :param list ImageReaderlist: list of Reader objects for images
    :param object proj: references a :class:`G2Project` project
    :param imageRef: A reference to the desired image in the project.
      Either the Image tree name (str), the image's index (int) or
      a image object (:class:`G2Image`)

    :return: array sumImg: corrected image for background/dark/flat back
    '''
    ImgObj = proj.image(imageRef)
    Controls = ImgObj.data['Image Controls']
    formatName = Controls.get('formatName','')
    imagefile = ImgObj.data['data'][1]
    if isinstance(imagefile, tuple) or isinstance(imagefile, list):
        imagefile, ImageTag =  imagefile # fix for multiimage files
    else:
        ImageTag = None # single-image file
    sumImg = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
    if sumImg is None:
        return []
    sumImg = np.array(sumImg,dtype=np.int32)
    darkImg = False
    if 'dark image' in Controls:
        darkImg,darkScale = Controls['dark image']
        if darkImg:
            dImgObj = proj.image(darkImg)
            formatName = dImgObj.data['Image Controls'].get('formatName','')
            imagefile = dImgObj.data['data'][1]
            if type(imagefile) is tuple:
                imagefile,ImageTag  = imagefile
            darkImage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if darkImg is None:
                raise Exception('Error reading dark image {}'.format(imagefile))
            sumImg += np.array(darkImage*darkScale,dtype=np.int32)
    if 'background image' in Controls:
        backImg,backScale = Controls['background image']
        if backImg:     #ignores any transmission effect in the background image
            bImgObj = proj.image(backImg)
            formatName = bImgObj.data['Image Controls'].get('formatName','')
            imagefile = bImgObj.data['data'][1]
            ImageTag = None # fix this for multiimage files
            backImage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if backImage is None:
                raise Exception('Error reading background image {}'.format(imagefile))
            if darkImg:
                backImage += np.array(darkImage*darkScale/backScale,dtype=np.int32)
            else:
                sumImg += np.array(backImage*backScale,dtype=np.int32)
    if 'Gain map' in Controls:
        gainMap = Controls['Gain map']
        if gainMap:
            gImgObj = proj.image(gainMap)
            formatName = gImgObj.data['Image Controls'].get('formatName','')
            imagefile = gImgObj.data['data'][1]
            ImageTag = None # fix this for multiimage files
            GMimage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if GMimage is None:
                raise Exception('Error reading Gain map image {}'.format(imagefile))
            sumImg = sumImg*GMimage/1000
    sumImg -= int(Controls.get('Flat Bkg',0))
    Imax = np.max(sumImg)
    Controls['range'] = [(0,Imax),[0,Imax]]
    return np.asarray(sumImg,dtype=np.int32)

def _constr_type(var):
    '''returns the constraint type based on phase/histogram use
    in a variable
    '''
    if var.histogram and var.phase:
        return 'HAP'
    elif var.phase:
        return 'Phase'
    elif var.histogram:
        return 'Hist'
    else:
        return 'Global'

class G2ObjectWrapper(object):
    """Base class for all GSAS-II object wrappers.

    The underlying GSAS-II format can be accessed as `wrapper.data`. A number
    of overrides are implemented so that the wrapper behaves like a dictionary.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    def __init__(self, datadict):
        self.data = datadict

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, key):
        return key in self.data

    def get(self, k, d=None):
        return self.data.get(k, d)

    def keys(self):
        return self.data.keys()

    def values(self):
        return self.data.values()

    def items(self):
        return self.data.items()


class G2Project(G2ObjectWrapper):
    """Represents an entire GSAS-II project. The object contains these
    class variables:

     * G2Project.filename: contains the .gpx filename
     * G2Project.names: contains the contents of the project "tree" as a list
       of lists. Each top-level entry in the tree is an item in the list. The
       name of the top-level item is the first item in the inner list. Children
       of that item, if any, are subsequent entries in that list.
     * G2Project.data: contains the entire project as a dict. The keys
       for the dict are the top-level names in the project tree (initial items
       in the G2Project.names inner lists) and each top-level item is stored
       as a dict.

         * The contents of Top-level entries will be found in the item
           named 'data', as an example, ``G2Project.data['Notebook']['data']``

         * The contents of child entries will be found in the item
           using the names of the parent and child, for example
           ``G2Project.data['Phases']['NaCl']``

    :param str gpxfile: Existing .gpx file to be loaded. If nonexistent,
            creates an empty project.
    :param str author: Author's name (not yet implemented)
    :param str newgpx: The filename the project should be saved to in
            the future. If both newgpx and gpxfile are present, the project is
            loaded from the file named by gpxfile and then when saved will
            be written to the file named by newgpx.
    :param str filename: To be deprecated. Serves the same function as newgpx,
            which has a somewhat more clear name.
            (Do not specify both newgpx and filename).

    There are two ways to initialize this object:

    >>> # Load an existing project file
    >>> proj = G2Project('filename.gpx')

    >>> # Create a new project
    >>> proj = G2Project(newgpx='new_file.gpx')

    Histograms can be accessed easily.

    >>> # By name
    >>> hist = proj.histogram('PWDR my-histogram-name')

    >>> # Or by index
    >>> hist = proj.histogram(0)
    >>> assert hist.id == 0

    >>> # Or by random id
    >>> assert hist == proj.histogram(hist.ranId)

    Phases can be accessed the same way.

    >>> phase = proj.phase('name of phase')

    New data can also be loaded via :meth:`~G2Project.add_phase` and
    :meth:`~G2Project.add_powder_histogram`.

    >>> hist = proj.add_powder_histogram('some_data_file.chi',
                                         'instrument_parameters.prm')
    >>> phase = proj.add_phase('my_phase.cif', histograms=[hist])

    Parameters for Rietveld refinement can be turned on and off at the project level
    as well as described in
    :meth:`~G2Project.set_refinement`, :meth:`~G2Project.iter_refinements` and
    :meth:`~G2Project.do_refinements`.
    """
    def __init__(self, gpxfile=None, author=None, filename=None, newgpx=None):
        if filename is not None and newgpx is not None:
            raise G2ScriptException('Do not use filename and newgpx together')
        elif filename is not None:
            G2fil.G2Print("Warning - recommending use of parameter newgpx rather than filename\n\twhen creating a G2Project")
        elif newgpx:
            filename = newgpx
        if not gpxfile:
            filename = os.path.abspath(os.path.expanduser(filename))
            self.filename = filename
            self.data, self.names = make_empty_project(author=author, filename=filename)
        elif os.path.exists(os.path.expanduser(gpxfile)):
            # TODO set author
            self.data, self.names = LoadDictFromProjFile(gpxfile)
            self.update_ids()
            if filename:
                filename = os.path.abspath(os.path.expanduser(filename))
                dr = os.path.split(filename)[0]
                if not os.path.exists(dr):
                    raise Exception("Directory {} for filename/newgpx does not exist".format(dr))
                self.filename = filename
            else:
                self.filename = os.path.abspath(os.path.expanduser(gpxfile))
        else:
            raise ValueError("Not sure what to do with gpxfile {}. Does not exist?".format(gpxfile))

    @classmethod
    def from_dict_and_names(cls, gpxdict, names, filename=None):
        """Creates a :class:`G2Project` directly from
        a dictionary and a list of names. If in doubt, do not use this.

        :returns: a :class:`G2Project`
        """
        out = cls()
        if filename:
            filename = os.path.abspath(os.path.expanduser(filename))
            out.filename = filename
            gpxdict['Controls']['data']['LastSavedAs'] = filename
        else:
            try:
                out.filename = gpxdict['Controls']['data']['LastSavedAs']
            except KeyError:
                out.filename = None
        out.data = gpxdict
        out.names = names

    def save(self, filename=None):
        """Saves the project, either to the current filename, or to a new file.

        Updates self.filename if a new filename provided"""
        controls_data = self.data['Controls']['data']
        # place info in .gpx about GSAS-II
        #  report versions of modules; ignore what has not been used
        modList = [np]
        try:
            modList += [mpl]
        except NameError:
            pass
        try:
            modList += [sp]
        except NameError:
            pass
        controls_data['PythonVersions'] = G2fil.get_python_versions(modList)
        #    G2 version info
        controls_data['LastSavedUsing'] = str(GSASIIpath.GetVersionNumber())
        try:
            if GSASIIpath.HowIsG2Installed().startswith('git'):
                g2repo = GSASIIpath.openGitRepo(GSASIIpath.path2GSAS2)
                commit = g2repo.head.commit
                controls_data['LastSavedUsing'] += f" git {commit.hexsha[:8]} script"
            else:
                gv = getSavedVersionInfo()
                if gv is not None:
                    Controls['LastSavedUsing'] += f" static {gv.git_version[:8]}"
        except:
            pass
        #    .gpx name
        if filename:
            filename = os.path.abspath(os.path.expanduser(filename))
            controls_data['LastSavedAs'] = filename
            self.filename = filename
        elif not self.filename:
            raise AttributeError("No file name to save to")
        SaveDictToProjFile(self.data, self.names, self.filename)

    def add_powder_histogram(self, datafile, iparams=None, phases=[],
                    fmthint=None, databank=None, instbank=None,
                    multiple=False, URL=False):
        """Loads a powder data histogram or multiple powder histograms
        into the project.

        Note that the data type (x-ray/CW neutron/TOF) for the histogram
        will be set from the instrument parameter file. The instrument
        geometry is assumed to be Debye-Scherrer except for
        dual-wavelength x-ray, where Bragg-Brentano is assumed.

        :param str datafile: A filename with the powder data file to read.
          Note that in unix fashion, "~" can be used to indicate the
          home directory (e.g. ~/G2data/data.fxye).
        :param str iparams: A filename for an instrument parameters file,
            or a pair of instrument parameter dicts from :func:`load_iprms`.
            This may be omitted for readers that provide the instrument
            parameters in the file. (Only a few importers do this.)
        :param list phases: A list of phases to link to the new histogram,
           phases can be references by object, name, rId or number.
           Alternately, use 'all' to link to all phases in the project.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu).
        :param int databank: Specifies a dataset number to read, if file contains
          more than set of data. This should be 1 to read the first bank in
          the file (etc.) regardless of the number on the Bank line, etc.
          Default is None which means the first dataset in the file is read.
          When multiple is True, optionally a list of dataset numbers can
          be supplied here.
        :param int instbank: Specifies an instrument parameter set to read, if
          the instrument parameter file contains more than set of parameters.
          This will match the INS # in an GSAS type file so it will typically
          be 1 to read the first parameter set in the file (etc.)
          Default is None which means there should only be one parameter set
          in the file.
        :param bool multiple: If False (default) only one dataset is read, but if
          specified as True, all selected banks of data (see databank)
          are read in.
        :param bool URL: if True, the contents of datafile is a URL and
          if not a dict, the contents of iparams is also a URL.
          both files will be downloaded to a temporary location 
          and read. The downloaded files will not be saved. 
          If URL is specified and the Python requests package is 
          not installed, a `ModuleNotFoundError` Exception will occur. 
          will occur.
        :returns: A :class:`G2PwdrData` object representing
            the histogram, or if multiple is True, a list of :class:`G2PwdrData`
            objects is returned.
        """
        LoadG2fil()
        if not URL:
            datafile = os.path.abspath(os.path.expanduser(datafile))
        pwdrreaders = import_generic(datafile, Readers['Pwdr'],fmthint=fmthint,
                                         bank=databank, URL=URL)
        if not multiple: pwdrreaders = pwdrreaders[0:1]
        histlist = []
        if URL:
            iparmfile = downloadFile(iparams)
        else:
            try:
                iparmfile = os.path.abspath(os.path.expanduser(iparams))
            except:
                pass
        for r in pwdrreaders:
            histname, new_names, pwdrdata = load_pwd_from_reader(r, iparmfile,
                                          [h.name for h in self.histograms()],bank=instbank)
            if histname in self.data:
                G2fil.G2Print("Warning - redefining histogram", histname)
            elif self.names[-1][0] == 'Phases':
                self.names.insert(-1, new_names)
            else:
                self.names.append(new_names)
            self.data[histname] = pwdrdata
            self.update_ids()

            if phases == 'all':
                phases = self.phases()
            for phase in phases:
                phase = self.phase(phase)
                self.link_histogram_phase(histname, phase)
            histlist.append(self.histogram(histname))

        if multiple:
            return histlist
        else:
            return histlist[0]

    def clone_powder_histogram(self, histref, newname, Y, Yerr=None):
        '''Creates a copy of a powder diffraction histogram with new Y values.
        The X values are not changed. The number of Y values must match the
        number of X values.

        :param histref: The histogram object, the name of the histogram (str), or ranId
           or histogram index.
        :param str newname: The name to be assigned to the new histogram
        :param list Y: A set of intensity values
        :param list Yerr: A set of uncertainties for the intensity values (may be None,
           sets all weights to unity)
        :returns: the new histogram object (type G2PwdrData)
        '''
        hist = self.histogram(histref)
        for i in self.names:
            if i[0] == hist.name:
                subkeys = i[1:]
                break
        else:
            raise Exception("error in self.names, hist not found")
        orighist = hist.name
        newhist = 'PWDR '+newname
        if len(Y) != len(self[orighist]['data'][1][0]):
            raise Exception("clone error: length of Y does not match number of X values ({})"
                                .format(len(self[orighist]['data'][1][0])))
        if Yerr is not None and len(Yerr) != len(self[orighist]['data'][1][0]):
            raise Exception("clone error: length of Yerr does not match number of X values ({})"
                                .format(len(self[orighist]['data'][1][0])))

        self[newhist] = copy.deepcopy(self[orighist])
        # intensities
        yo = self[newhist]['data'][1][1] = ma.MaskedArray(Y,mask=self[orighist]['data'][1][1].mask)

        Ymin,Ymax = yo.min(),yo.max()
        # set to zero: weights, calc, bkg, obs-calc
        for i in [2,3,4,5]:
            self[newhist]['data'][1][i] *= 0
        # weights
        if Yerr is not None:
            self[newhist]['data'][1][2] += 1./np.array(Yerr)**2
        else:
            self[newhist]['data'][1][2] += 1            # set all weights to 1
        self[newhist]['data'][0] = {'wtFactor': 1.0, 'Dummy': False,
                  'ranId': ran.randint(0, sys.maxsize),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'Yminmax': [Ymin, Ymax]}
        self[newhist]['Comments'].insert(0,'Cloned from '+orighist)
        self[newhist]['Reflection Lists'] = {}
        self[newhist]['Index Peak List'] = [[], []]
        self[newhist]['Unit Cells List'] = []
        self[newhist]['Peak List'] = {'peaks': [], 'sigDict': {}}
        self.names.append([newhist]+subkeys)
        self.update_ids()
        return self.histogram(newhist)

    def add_simulated_powder_histogram(self, histname, iparams, Tmin, Tmax, Tstep=None,
                                       wavelength=None, scale=None, phases=[], ibank=None,
                                           Npoints=None):
        """Create a simulated powder data histogram for the project.

        Requires an instrument parameter file.
        Note that in unix fashion, "~" can be used to indicate the
        home directory (e.g. ~/G2data/data.prm). The instrument parameter file
        will determine if the histogram is x-ray, CW neutron, TOF, etc. as well
        as the instrument type.

        :param str histname: A name for the histogram to be created.
        :param str iparams: The instrument parameters file, a filename.
        :param float Tmin: Minimum 2theta or TOF (millisec) for dataset to be simulated
        :param float Tmax: Maximum 2theta or TOF (millisec) for dataset to be simulated
        :param float Tstep: Step size in 2theta or deltaT/T (TOF) for simulated dataset.
           Default is to compute this from Npoints.
        :param float wavelength: Wavelength for CW instruments, overriding the value
           in the instrument parameters file if specified. For single-wavelength histograms,
           this should be a single float value, for K alpha 1,2 histograms, this should
           be a list or tuple with two values.
        :param float scale: Histogram scale factor which multiplies the pattern. Note that
           simulated noise is added to the pattern, so that if the maximum intensity is
           small, the noise will mask the computed pattern. The scale needs to be a large
           number for neutrons.
           The default, None, provides a scale of 1 for x-rays, 10,000 for CW neutrons
           and 100,000 for TOF.
        :param list phases: Phases to link to the new histogram. Use proj.phases() to link to
           all defined phases.
        :param int ibank: provides a bank number for the instrument parameter file. The
           default is None, corresponding to load the first bank.
        :param int Νpoints: the number of data points to be used for computing the
            diffraction pattern. Defaults as None, which sets this to 2500. Do not specify
            both Npoints and Tstep. Due to roundoff the actual number of points used may differ
            by +-1 from Npoints. Must be below 25,000.

        :returns: A :class:`G2PwdrData` object representing the histogram
        """
        LoadG2fil()
        iparams = os.path.abspath(os.path.expanduser(iparams))
        if not os.path.exists(iparams):
            raise G2ScriptException("File does not exist:"+iparams)
        rd = G2obj.ImportPowderData( # Initialize a base class reader
            extensionlist=tuple(),
            strictExtension=False,
            formatName = 'Simulate dataset',
            longFormatName = 'Compute a simulated pattern')
        rd.powderentry[0] = '' # no filename
        rd.powderentry[2] = 1 # only one bank
        rd.comments.append('This is a dummy dataset for powder pattern simulation')
        rd.idstring = histname
        #Iparm1, Iparm2 = load_iprms(iparams, rd)
        if Tmax < Tmin:
            Tmin,Tmax = Tmax,Tmin
        if Tstep is not None and Npoints is not None:
            raise G2ScriptException("Error: Tstep and Npoints both specified")
        elif Tstep is not None:
            Tstep = abs(Tstep)
        elif Npoints is None:
            Npoints = 2500
        Iparm1, Iparm2 = load_iprms(iparams, rd, bank=ibank)
        #G2fil.G2Print('Instrument parameters read:',reader.instmsg)
        if 'T' in Iparm1['Type'][0]:
            # patch -- anticipate TOF values in microsec from buggy version
            if Tmax > 200.:
                print('Error: Tmax is too large. Note that input for TOF Tmin & Tmax has changed.')
                print('       Tmin & Tmax are now in milliseconds not microsec. Step is now deltaT/T.')
                raise G2ScriptException("Error: Tmax is too large")
            if Npoints:
                N = Npoints
                Tstep = (np.log(Tmax)-np.log(Tmin))/N
            else:
                N = (np.log(Tmax)-np.log(Tmin))/Tstep
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.exp((np.arange(0,N))*Tstep+np.log(Tmin*1000.))
            N = len(x)
            unit = 'millisec'
        else:
            if Npoints:
                N = Npoints
            else:
                N = int((Tmax-Tmin)/Tstep)+1
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.linspace(Tmin,Tmax,N,True)
            N = len(x)
            unit = 'degrees 2theta'
        if N < 3:
            raise G2ScriptException("Error: Range is too small or step is too large, <3 points")
        G2fil.G2Print('Simulating {} points from {} to {} {}'.format(N,Tmin,Tmax,unit))
        rd.powderdata = [
            np.array(x), # x-axis values
            np.zeros_like(x), # powder pattern intensities
            np.ones_like(x), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(x), # calc. intensities (zero)
            np.zeros_like(x), # calc. background (zero)
            np.zeros_like(x), # obs-calc profiles
            ]
        Tmin = rd.powderdata[0][0]
        Tmax = rd.powderdata[0][-1]
        histname, new_names, pwdrdata = load_pwd_from_reader(rd, iparams,
                                            [h.name for h in self.histograms()],ibank)
        if histname in self.data:
            G2fil.G2Print("Warning - redefining histogram", histname)
        elif self.names[-1][0] == 'Phases':
            self.names.insert(-1, new_names)
        else:
            self.names.append(new_names)
        if scale is not None:
            pwdrdata['Sample Parameters']['Scale'][0] = scale
        elif pwdrdata['Instrument Parameters'][0]['Type'][0].startswith('PNC'):
            pwdrdata['Sample Parameters']['Scale'][0] = 10000.
        elif pwdrdata['Instrument Parameters'][0]['Type'][0].startswith('PNT'):
            pwdrdata['Sample Parameters']['Scale'][0] = 100000.
        if wavelength is not None:
            if 'Lam1' in pwdrdata['Instrument Parameters'][0]:
                # have alpha 1,2 here
                errmsg = ''
                try:
                    pwdrdata['Instrument Parameters'][0]['Lam1'][0] = float(wavelength[0])
                    pwdrdata['Instrument Parameters'][0]['Lam1'][1] = float(wavelength[0])
                    pwdrdata['Instrument Parameters'][0]['Lam2'][0] = float(wavelength[1])
                    pwdrdata['Instrument Parameters'][0]['Lam2'][1] = float(wavelength[1])
                except:
                    errmsg = "add_simulated_powder_histogram Error: only one wavelength with alpha 1+2 histogram?"
                if errmsg: raise G2ScriptException(errmsg)
            elif 'Lam' in pwdrdata['Instrument Parameters'][0]:
                errmsg = ''
                try:
                    pwdrdata['Instrument Parameters'][0]['Lam'][0] = float(wavelength)
                    pwdrdata['Instrument Parameters'][0]['Lam'][1] = float(wavelength)
                except:
                    errmsg = "add_simulated_powder_histogram Error: invalid wavelength?"
                if errmsg: raise G2ScriptException(errmsg)
            else:
                raise G2ScriptException("add_simulated_powder_histogram Error: can't set a wavelength for a non-CW dataset")
        self.data[histname] = pwdrdata
        self.update_ids()

        for phase in phases:
            phase = self.phase(phase)
            self.link_histogram_phase(histname, phase)

        self.set_Controls('cycles', 0)
        self.data[histname]['data'][0]['Dummy'] = True
        return self.histogram(histname)

    def add_phase(self, phasefile=None, phasename=None, histograms=[],
                      fmthint=None, mag=False,
                      spacegroup='P 1',cell=None, URL=False):
        """Loads a phase into the project, usually from a .cif file

        :param str phasefile: The CIF file (or other file type, see fmthint)
          that the phase will be read from.
          May be left as None (the default) if the phase will be constructed
          a step at a time.
        :param str phasename: The name of the new phase, or None for the
          default. A phasename must be specified when a phasefile is not.
        :param list histograms: The names of the histograms to associate with
            this phase. Use proj.histograms() to add to all histograms.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu). Specifying this
          is optional but is strongly encouraged.
        :param bool mag: Set to True to read a magCIF
        :param str spacegroup: The space group name as a string. The
          space group must follow the naming rules used in
          :func:`GSASIIspc.SpcGroup`. Defaults to 'P 1'. Note that
          this is only used when phasefile is None.
        :param list cell: a list with six unit cell constants
            (a, b, c, alpha, beta and gamma in Angstrom/degrees).
        :param bool URL: if True, the contents of phasefile is a URL 
          and the file will be downloaded to a temporary location 
          and read. The downloaded file will not be saved. 
          If URL is specified and the Python requests package is 
          not installed, a `ModuleNotFoundError` Exception will occur. 
          will occur. 

        :returns: A :class:`G2Phase` object representing the
            new phase.
        """
        LoadG2fil()
        histograms = [self.histogram(h).name for h in histograms]
        if phasefile is None:  # create a phase, rather than reading it
            if phasename is None:
                raise Exception('add_phase: phasefile and phasename cannot both be None')
            phaseNameList = [p.name for p in self.phases()]
            phasename = G2obj.MakeUniqueLabel(phasename, phaseNameList)
            self.data['Phases'] = self.data.get('Phases', {'data': None})
            err,SGData=G2spc.SpcGroup(spacegroup)
            if err != 0:
                print('Space group error:', G2spc.SGErrors(err))
                raise Exception('Space group error')
            self.data['Phases'][phasename] = G2obj.SetNewPhase(
                Name=phasename,SGData=SGData)
            self.data['Phases'][phasename]['General']['Name'] = phasename
            for hist in histograms:
                self.link_histogram_phase(hist, phasename)

            for obj in self.names:
                if obj[0] == 'Phases':
                    phasenames = obj
                    break
            else:
                phasenames = ['Phases']
                self.names.append(phasenames)
            phasenames.append(phasename)

            data = self.data['Phases'][phasename]
            if cell:
                if len(cell) != 6:
                    raise Exception('Error: Unit cell must have 6 entries')
                data['General']['Cell'][1:7] = cell
                data['General']['Cell'][7] = G2lat.calc_V(G2lat.cell2A(cell))
            SetupGeneral(data, None)
            self.index_ids()

            self.update_ids()
            return self.phase(phasename)

        if not URL: 
            phasefile = os.path.abspath(os.path.expanduser(phasefile))
            if not os.path.exists(phasefile):
                raise G2ImportException(f'File {phasefile} does not exist')
        else:
            print(f'reading phase from URL {phasefile}')
        # TODO: handle multiple phases in a file
        phasereaders = import_generic(phasefile, Readers['Phase'],
                                          fmthint=fmthint, URL=URL)
        phasereader = phasereaders[0]

        if phasereader.MPhase and mag:
            phObj = phasereader.MPhase
        else:
            phObj = phasereader.Phase

        phasename = phasename or phObj['General']['Name']
        phaseNameList = [p.name for p in self.phases()]
        phasename = G2obj.MakeUniqueLabel(phasename, phaseNameList)
        phObj['General']['Name'] = phasename

        if 'Phases' not in self.data:
            self.data['Phases'] = { 'data': None }
        assert phasename not in self.data['Phases'], "phase names should be unique"
        self.data['Phases'][phasename] = phObj

        # process constraints, currently generated only from ISODISTORT CIFs
        if phasereader.Constraints:
            Constraints = self.data['Constraints']
            for i in phasereader.Constraints:
                if isinstance(i, dict):
                    if '_Explain' not in Constraints:
                        Constraints['_Explain'] = {}
                    Constraints['_Explain'].update(i)
                else:
                    Constraints['Phase'].append(i)

        data = self.data['Phases'][phasename]
#        generalData = data['General']
#        SGData = generalData['SGData']
#        NShkl = len(G2spc.MustrainNames(SGData))
#        NDij = len(G2spc.HStrainNames(SGData))
#        Super = generalData.get('Super', 0)
#        if Super:
#            SuperVec = np.array(generalData['SuperVec'][0])
#        else:
#            SuperVec = []
#        UseList = data['Histograms']

        for hist in histograms:
            self.link_histogram_phase(hist, phasename)

        for obj in self.names:
            if obj[0] == 'Phases':
                phasenames = obj
                break
        else:
            phasenames = ['Phases']
            self.names.append(phasenames)
        phasenames.append(phasename)

        # TODO should it be self.filename, not phasefile?
        SetupGeneral(data, os.path.dirname(phasefile))
        self.index_ids()

        self.update_ids()
        return self.phase(phasename)

    def add_single_histogram(self, datafile, phase=None, fmthint=None):
        """Loads a powder data histogram or multiple powder histograms
        into the project.

        :param str datafile: A filename with the single crystal data file
          to read. Note that in unix fashion, "~" can be used to indicate the
          home directory (e.g. ~/G2data/data.hkl).
        :param phases: A phase to link to the new histogram. A
           phase can be referenced by object, name, rId or number.
           If not specified, no phase will be linked.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified,
          an error will be generated, as the file format will not distinguish
          well between different data types.
        :returns: A :class:`G2Single` object representing
            the histogram
        """
        LoadG2fil()
        datafile = os.path.abspath(os.path.expanduser(datafile))
        if fmthint is None:
            choices = ', '.join(['"'+i.formatName+'"' for i in Readers['HKLF']])
            raise ValueError("add_single_histogram: A value is required for the fmthint parameter"+
                                 f'\n\nChoices are: {choices}')
        readers = import_generic(datafile, Readers['HKLF'],fmthint=fmthint)
        #histlist = []
        for r in readers:  # only expect 1
            histname = 'HKLF ' + G2obj.StripUnicode(
                os.path.split(r.readfilename)[1],'_')
            #HistName = G2obj.MakeUniqueLabel(HistName, existingnames)
            newnames = [histname,'Instrument Parameters','Reflection List']
            if histname in self.data:
                G2fil.G2Print("Warning - redefining histogram", histname)
            elif self.names[-1][0] == 'Phases':
                self.names.insert(-1, newnames)  # add histogram into tree before phases
            else:
                self.names.append(newnames)
            data = [{
                'wtFactor': 1.0,
                'Dummy': False,
                'ranId': ran.randint(0, sys.maxsize),
                'histTitle': ''},r.RefDict]
            self.data[histname] = {
                'data':data,
                'Instrument Parameters':r.Parameters,
                'Reflection List': {}
                }
            self.update_ids()
            if phase is None: return
            self.link_histogram_phase(histname, phase)

    def link_histogram_phase(self, histogram, phase):
        """Associates a given histogram and phase.

        .. seealso::

            :meth:`~G2Project.histogram`
            :meth:`~G2Project.phase`"""
        hist = self.histogram(histogram)
        phase = self.phase(phase)

        generalData = phase['General']

        if hist.name.startswith('HKLF '):
            G2mth.UpdateHKLFvals(hist.name, phase, hist.data['data'][1])
        elif hist.name.startswith('PWDR '):
            hist['Reflection Lists'][generalData['Name']] = {}
            UseList = phase['Histograms']
            SGData = generalData['SGData']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            UseList[hist.name] = G2mth.SetDefaultDData(
                'PWDR', hist.name, NShkl=NShkl, NDij=NDij)
            UseList[hist.name]['hId'] = hist.id
            for key, val in [('Use', True), ('LeBail', False),
                             ('Babinet', {'BabA': [0.0, False],
                                          'BabU': [0.0, False]})]:
                if key not in UseList[hist.name]:
                    UseList[hist.name][key] = val
        else:
            raise RuntimeError("Unexpected histogram" + hist.name)

    def reload(self):
        """Reload self from self.filename"""
        data, names = LoadDictFromProjFile(self.filename)
        self.names = names
        # Need to deep copy the new data file data into the current tree,
        # so that any existing G2Phase, or G2PwdrData objects will still be
        # valid
        _deep_copy_into(from_=data, into=self.data)

    def refine(self, newfile=None, printFile=None, makeBack=False):
        '''Invoke a refinement for the project. The project is written to
        the currently selected gpx file and then either a single or sequential refinement
        is performed depending on the setting of 'Seq Data' in Controls
        (set in :meth:`get_Controls`).
        '''
        seqSetting = self.data['Controls']['data'].get('Seq Data',[])
        if not seqSetting:
            self.index_ids()    # index_ids will automatically save the project
            # TODO: migrate to RefineCore G2strMain does not properly use printFile
            # G2strMain.RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
            #      calcControls,pawleyLookup,ifPrint,printFile,dlg)

            # check that constraints are OK
            errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.filename)
            if errmsg:
                G2fil.G2Print('Constraint error',errmsg)
                raise Exception('Constraint error')
            if warnmsg:
                G2fil.G2Print('\nNote these constraint warning(s):\n'+warnmsg)
                G2fil.G2Print('Generated constraints\n'+G2mv.VarRemapShow([],True))
            G2strMain.Refine(self.filename, makeBack=makeBack)
        else:
            self._seqrefine()
        self.reload() # get file from GPX

    def _seqrefine(self):
        '''Perform a sequential refinement.
        '''

        self.data['Controls']['data']['ShowCell'] = True
        # add to tree item to project, if not present
        if 'Sequential results' not in self.data:
            self.data['Sequential results'] = {'data':{}}
            self.names.append(['Sequential results'])
        self.index_ids()    # index_ids will automatically save the project
        #GSASIIpath.IPyBreak_base()

        # check that constraints are OK
        errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.filename)
        if errmsg:
            G2fil.G2Print('Constraint error',errmsg)
            raise Exception('Constraint error')
        if warnmsg:
            G2fil.G2Print('\nNote these constraint warning(s):\n'+warnmsg)
            G2fil.G2Print('Generated constraints\n'+G2mv.VarRemapShow([],True))
        OK,Msg = G2strMain.SeqRefine(self.filename,None)

    def histogram(self, histname):
        """Returns the histogram object associated with histname, or None
        if it does not exist.

        :param histname: The name of the histogram (str), or ranId or
          (for powder) the histogram index.
        :returns: A :class:`G2PwdrData` object, or :class:`G2Single` object, or
           None if the histogram does not exist

        .. seealso::
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
        """
        if isinstance(histname, G2PwdrData) or isinstance(histname, G2Single):
            if histname.proj == self:
                return histname
            else:
                raise Exception(f'Histogram object {histname} (type {type(histname)}) is project {histname.proj}')
        if histname in self.data:
            if histname.startswith('HKLF '):
                return G2Single(self.data[histname], self, histname)
            elif histname.startswith('PWDR '):
                return G2PwdrData(self.data[histname], self, histname)
        try:
            # see if histname is an id or ranId
            histname = int(histname)
        except ValueError:
            return
        for histogram in self.histograms():
            if histogram.name.startswith('HKLF '):
                if histogram.data['data'][0]['ranId'] == histname:
                    return histogram
            elif histogram.name.startswith('PWDR '):
                if histogram.id == histname or histogram.ranId == histname:
                    return histogram

    def histType(self, histname):
        """Returns the type for histogram object associated with histname, or
        None if it does not exist.

        :param histname: The name of the histogram (str), or ranId or
          (for powder) the histogram index.
        :returns: 'PWDR' for a Powder histogram,
           'HKLF' for a single crystal histogram, or
           None if the histogram does not exist

        .. seealso::
            :meth:`~G2Project.histogram`
        """
        if isinstance(histname, G2PwdrData):
            if histname.proj == self:
                return 'PWDR'
            else:
                raise Exception(f'Histogram object {histname} (type {type(histname)}) is project {histname.proj}')
        if isinstance(histname, G2Single):
            if histname.proj == self:
                return 'HKLF'
            else:
                raise Exception(f'Histogram object {histname} (type {type(histname)}) is project {histname.proj}')
        if histname in self.data:
            if histname.startswith('HKLF '):
                return 'HKLF'
            elif histname.startswith('PWDR '):
                return 'PWDR'
        try:
            # see if histname is an id or ranId
            histname = int(histname)
        except ValueError:
            return
        for histogram in self.histograms():
            if histogram.name.startswith('HKLF '):
                if histogram.data['data'][0]['ranId'] == histname:
                    return 'HKLF'
            elif histogram.name.startswith('PWDR '):
                if histogram.id == histname or histogram.ranId == histname:
                    return 'PWDR'

    def histograms(self, typ=None):
        """Return a list of all histograms, as :class:`G2PwdrData` objects

        For now this only finds Powder/Single Xtal histograms, since that is all that is
        currently implemented in this module.

        :param ste typ: The prefix (type) the histogram such as 'PWDR ' for
          powder or 'HKLF ' for single crystal. If None
          (the default) all known histograms types are found.
        :returns: a list of objects

        .. seealso::
            :meth:`~G2Project.histogram`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
        """
        output = []
        # loop through each tree entry. If it is more than one level (more than one item in the
        # list of names). then it must be a histogram, unless labeled Phases or Restraints
        if typ is None:
            for obj in self.names:
                if obj[0].startswith('PWDR ') or obj[0].startswith('HKLF '):
                    output.append(self.histogram(obj[0]))
        else:
            for obj in self.names:
                if len(obj) > 1 and obj[0].startswith(typ):
                    output.append(self.histogram(obj[0]))
        return output

    def phase(self, phasename):
        """
        Gives an object representing the specified phase in this project.

        :param str phasename: A reference to the desired phase. Either the phase
            name (str), the phase's ranId, the phase's index (both int) or
            a phase object (:class:`G2Phase`)
        :returns: A :class:`G2Phase` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
            """
        if isinstance(phasename, G2Phase):
            if phasename.proj == self:
                return phasename
        phases = self.data['Phases']
        if phasename in phases:
            return G2Phase(phases[phasename], phasename, self)

        try:
            # phasename should be phase index or ranId
            phasename = int(phasename)
        except ValueError:
            return

        for phase in self.phases():
            if phase.id == phasename or phase.ranId == phasename:
                return phase

    def phases(self):
        """
        Returns a list of all the phases in the project.

        :returns: A list of :class:`G2Phase` objects

        .. seealso::
            :meth:`~G2Project.histogram`
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            """
        for obj in self.names:
            if obj[0] == 'Phases':
                return [self.phase(p) for p in obj[1:]]
        return []

    def _images(self):
        """Returns a list of all the images in the project.
        """
        return [i[0] for i in self.names if i[0].startswith('IMG ')]

    def image(self, imageRef):
        """
        Gives an object representing the specified image in this project.

        :param str imageRef: A reference to the desired image. Either the Image
          tree name (str), the image's index (int) or
          a image object (:class:`G2Image`)
        :returns: A :class:`G2Image` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.images`
        """
        if isinstance(imageRef, G2Image):
            if imageRef.proj == self:
                return imageRef
            else:
                raise Exception("Image {} not in current selected project".format(imageRef.name))
        if imageRef in self._images():
            return G2Image(self.data[imageRef], imageRef, self)

        errmsg = ''
        try:
            # imageRef should be an index
            num = int(imageRef)
            imageRef = self._images()[num]
            return G2Image(self.data[imageRef], imageRef, self)
        except ValueError:
            errmsg = "imageRef {imageRef} not an object, name or image index in current selected project"
        except IndexError:
            errmsg = "imageRef {imageRef} out of range (max={len(self._images())-1)}) in current selected project"
        if errmsg: raise G2ScriptException(errmsg)

    def images(self):
        """
        Returns a list of all the images in the project.

        :returns: A list of :class:`G2Image` objects
        """
        return [G2Image(self.data[i],i,self) for i in self._images()]

    def _pdfs(self):
        """Returns a list of all the PDF entries in the project.
        """
        return [i[0] for i in self.names if i[0].startswith('PDF ')]

    def pdf(self, pdfRef):
        """
        Gives an object representing the specified PDF entry in this project.

        :param pdfRef: A reference to the desired image. Either the PDF
          tree name (str), the pdf's index (int) or
          a PDF object (:class:`G2PDF`)
        :returns: A :class:`G2PDF` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.pdfs`
            :class:`~G2PDF`
        """
        if isinstance(pdfRef, G2PDF):
            if pdfRef.proj == self:
                return pdfRef
            else:
                raise Exception(f"PDF {pdfRef.name} not in current selected project")
        if pdfRef in self._pdfs():
            return G2PDF(self.data[pdfRef], pdfRef, self)

        errmsg = ''
        try:
            # pdfRef should be an index
            num = int(pdfRef)
            pdfRef = self._pdfs()[num]
            return G2PDF(self.data[pdfRef], pdfRef, self)
        except ValueError:
            errmsg = f"pdfRef {pdfRef} not an object, name or PDF index in current selected project"
        except IndexError:
            errmsg = f"pdfRef {pdfRef} out of range (max={len(G2SmallAngle)-1}) in current selected project"
        if errmsg: raise Exception(errmsg)
    def pdfs(self):
        """
        Returns a list of all the PDFs in the project.

        :returns: A list of :class:`G2PDF` objects
        """
        return [G2PDF(self.data[i],i,self) for i in self._pdfs()]

    def copy_PDF(self, PDFobj, histogram):
        '''Creates a PDF entry that can be used to compute a PDF
        as a copy of settings in an existing PDF (:class:`G2PDF`)
        object.
        This places an entry in the project but :meth:`G2PDF.calculate`
        must be used to actually perform the PDF computation.

        :param PDFobj: A :class:`G2PDF` object which may be
          in a separate project or the dict associated with the
          PDF object (G2PDF.data).
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :returns: A :class:`G2PDF` object for the PDF entry
        '''
        LoadG2fil()
        PDFname = 'PDF ' + self.histogram(histogram).name[4:]
        PDFdict = {'data':None}
        for i in 'PDF Controls', 'PDF Peaks':
            PDFdict[i] = copy.deepcopy(PDFobj[i])
        self.names.append([PDFname]+['PDF Controls', 'PDF Peaks'])
        self.data[PDFname] = PDFdict
        for i in 'I(Q)','S(Q)','F(Q)','G(R)','g(r)':
            self.data[PDFname]['PDF Controls'][i] = []
        G2fil.G2Print('Adding "{}" to project'.format(PDFname))
        return G2PDF(self.data[PDFname], PDFname, self)

    def add_PDF(self, prmfile, histogram):
        '''Creates a PDF entry that can be used to compute a PDF.
        Note that this command places an entry in the project,
        but :meth:`G2PDF.calculate` must be used to actually perform
        the computation.

        :param str datafile: The powder data file to read, a filename.
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :returns: A :class:`G2PDF` object for the PDF entry
        '''

        LoadG2fil()
        PDFname = 'PDF ' + self.histogram(histogram).name[4:]
        peaks = {'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]}
        Controls = {
            'Sample':{'Name':self.histogram(histogram).name,'Mult':1.0},
            'Sample Bkg.':{'Name':'','Mult':-1.0,'Refine':False},
            'Container':{'Name':'','Mult':-1.0,'Refine':False},
            'Container Bkg.':{'Name':'','Mult':-1.0},
            'ElList':{},
            'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':0.0,'Flat Bkg':0,
            'DetType':'Area detector','ObliqCoeff':0.2,'Ruland':0.025,'QScaleLim':[20,25],
            'Lorch':False,'BackRatio':0.0,'Rmax':100.,'noRing':False,'IofQmin':1.0,'Rmin':1.0,
            'I(Q)':[],'S(Q)':[],'F(Q)':[],'G(R)':[],'g(r)':[]}

        fo = open(prmfile,'r')
        S = fo.readline()
        while S:
            if '#' in S:
                S = fo.readline()
                continue
            key,val = S.split(':',1)
            try:
                Controls[key] = eval(val)
            except:
                Controls[key] = val.strip()
            S = fo.readline()
        fo.close()
        Controls['Sample']['Name'] = self.histogram(histogram).name
        for i in 'Sample Bkg.','Container','Container Bkg.':
            Controls[i]['Name'] = ''
        PDFdict = {'data':None,'PDF Controls':Controls, 'PDF Peaks':peaks}
        self.names.append([PDFname]+['PDF Controls', 'PDF Peaks'])
        self.data[PDFname] = PDFdict
        G2fil.G2Print('Adding "{}" to project'.format(PDFname))
        return G2PDF(self.data[PDFname], PDFname, self)

    def seqref(self):
        """
        Returns a sequential refinement results object, if present

        :returns: A :class:`G2SeqRefRes` object or None if not present
        """
        if 'Sequential results' not in self.data: return
        return G2SeqRefRes(self.data['Sequential results']['data'], self)

    def update_ids(self):
        """Makes sure all phases and histograms have proper hId and pId"""
        # Translated from GetUsedHistogramsAndPhasesfromTree,
        #   GSASIIdataGUI.py:4107
        for i, h in enumerate(self.histograms()):
            h.id = i
        for i, p in enumerate(self.phases()):
            p.id = i

    def do_refinements(self, refinements=[{}], histogram='all', phase='all',
                       outputnames=None, makeBack=False):
        """Conducts one or a series of refinements according to the
           input provided in parameter refinements. This is a wrapper
           around :meth:`iter_refinements`

        :param list refinements: A list of dictionaries specifiying changes to be made to
            parameters before refinements are conducted.
            See the :ref:`Refinement_recipe` section for how this is defined.
            If not specified, the default value is ``[{}]``, which performs a single
            refinement step is performed with the current refinement settings.
        :param str histogram: Name of histogram for refinements to be applied
            to, or 'all'; note that this can be overridden for each refinement
            step via a "histograms" entry in the dict.
        :param str phase: Name of phase for refinements to be applied to, or
            'all'; note that this can be overridden for each refinement
            step via a "phases" entry in the dict.
        :param list outputnames: Provides a list of project (.gpx) file names
            to use for each refinement step (specifying None skips the save step).
            See :meth:`save`.
            Note that this can be overridden using an "output" entry in the dict.
        :param bool makeBack: determines if a backup ).bckX.gpx) file is made
            before a refinement is performed. The default is False.

        To perform a single refinement without changing any parameters, use this
        call:

        .. code-block::  python

            my_project.do_refinements([])
        """

        for proj in self.iter_refinements(refinements, histogram, phase,
                                          outputnames, makeBack):
            pass
        return self

    def iter_refinements(self, refinements, histogram='all', phase='all',
                         outputnames=None, makeBack=False):
        """Conducts a series of refinements, iteratively. Stops after every
        refinement and yields this project, to allow error checking or
        logging of intermediate results. Parameter use is the same as for
        :meth:`do_refinements` (which calls this method).

        >>> def checked_refinements(proj):
        ...     for p in proj.iter_refinements(refs):
        ...         # Track intermediate results
        ...         log(p.histogram('0').residuals)
        ...         log(p.phase('0').get_cell())
        ...         # Check if parameter diverged, nonsense answer, or whatever
        ...         if is_something_wrong(p):
        ...             raise Exception("I need a human!")


        """
        if outputnames:
            if len(refinements) != len(outputnames):
                raise ValueError("Should have same number of outputs to"
                                 "refinements")
        else:
            outputnames = [None for r in refinements]

        for output, refinedict in zip(outputnames, refinements):
            if 'histograms' in refinedict:
                hist = refinedict['histograms']
            else:
                hist = histogram
            if 'phases' in refinedict:
                ph = refinedict['phases']
            else:
                ph = phase
            if 'output' in refinedict:
                output = refinedict['output']
            self.set_refinement(refinedict, hist, ph)
            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinedict:
                temp = {'set': refinedict['once']}
                self.set_refinement(temp, hist, ph)

            if output:
                self.save(output)

            if 'skip' not in refinedict:
                self.refine(makeBack=makeBack)
            yield self

            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinedict:
                temp = {'clear': refinedict['once']}
                self.set_refinement(temp, hist, ph)
            if 'call' in refinedict:
                fxn = refinedict['call']
                if callable(fxn):
                    fxn(*refinedict.get('callargs',[self]))
                elif callable(eval(fxn)):
                    eval(fxn)(*refinedict.get('callargs',[self]))
                else:
                    raise G2ScriptException("Error: call value {} is not callable".format(fxn))

    def set_refinement(self, refinement, histogram='all', phase='all'):
        """Set refinment flags at the project level to specified histogram(s)
        or phase(s).

        :param dict refinement: The refinements to be conducted
        :param histogram: Specifies either 'all' (default), a single histogram or
          a list of histograms. Histograms may be specified as histogram objects
          (see :class:`G2PwdrData`), the histogram name (str) or the index number (int)
          of the histogram in the project, numbered starting from 0.
          Omitting the parameter or the string 'all' indicates that parameters in
          all histograms should be set.
        :param phase: Specifies either 'all' (default), a single phase or
          a list of phases. Phases may be specified as phase objects
          (see :class:`G2Phase`), the phase name (str) or the index number (int)
          of the phase in the project, numbered starting from 0.
          Omitting the parameter or the string 'all' indicates that parameters in
          all phases should be set.

        Note that refinement parameters are categorized as one of three types:

        1. Histogram parameters
        2. Phase parameters
        3. Histogram-and-Phase (HAP) parameters

        .. seealso::
            :meth:`G2PwdrData.set_refinements`
            :meth:`G2PwdrData.clear_refinements`
            :meth:`G2Phase.set_refinements`
            :meth:`G2Phase.clear_refinements`
            :meth:`G2Phase.set_HAP_refinements`
            :meth:`G2Phase.clear_HAP_refinements`
            :meth:`G2Single.set_refinements`
        """

        if histogram == 'all':
            hists = self.histograms()
        elif isinstance(histogram, list) or isinstance(histogram, tuple):
            hists = []
            for h in histogram:
                if isinstance(h, str) or isinstance(h, int):
                    hists.append(self.histogram(h))
                else:
                    hists.append(h)
        elif isinstance(histogram, str) or isinstance(histogram, int):
            hists = [self.histogram(histogram)]
        else:
            hists = [histogram]

        if phase == 'all':
            phases = self.phases()
        elif isinstance(phase, list) or isinstance(phase, tuple):
            phases = []
            for ph in phase:
                if isinstance(ph, str) or isinstance(ph, int):
                    phases.append(self.phase(ph))
                else:
                    phases.append(ph)
        elif isinstance(phase, str) or isinstance(phase, int):
            phases = [self.phase(phase)]
        else:
            phases = [phase]

        pwdr_set = {}
        phase_set = {}
        hap_set = {}
        xtal_set = {}
        # divide the refinement set commands across the various objects
        for key, val in refinement.get('set', {}).items():
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_set[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_set[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                if key == 'PhaseFraction': key = 'Scale'
                hap_set[key] = val
            elif G2Single.is_valid_refinement_key(key):
                xtal_set[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            if isinstance(hist, G2PwdrData):
                hist.set_refinements(pwdr_set)
            elif isinstance(hist, G2Single):
                hist.set_refinements(xtal_set)
            else:
                raise KeyError(f"set_refinement error: unknown hist type {hist}")
        for phase in phases:
            phase.set_refinements(phase_set)
        for phase in phases:
            phase.set_HAP_refinements(hap_set, hists)

        pwdr_clear = {}
        phase_clear = {}
        hap_clear = {}
        xtal_clear = {}
        for key, val in refinement.get('clear', {}).items():
            # Clear refinement options
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_clear[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_clear[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                if key == 'PhaseFraction': key = 'Scale'
                hap_clear[key] = val   # was _set, seems wrong
            elif G2Single.is_valid_refinement_key(key):
                xtal_clear[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            if isinstance(hist, G2PwdrData):
                hist.clear_refinements(pwdr_clear)
            else:
                hist.clear_refinements(xtal_clear)
        for phase in phases:
            phase.clear_refinements(phase_clear)
        for phase in phases:
            phase.clear_HAP_refinements(hap_clear, hists)

    def index_ids(self):
        self.save()
        return G2stIO.GetUsedHistogramsAndPhases(self.filename)

    def _sasd(self):
        """Returns a list of all the SASD entries in the project.
        """
        return [i[0] for i in self.names if i[0].startswith('SASD ')]

    def SAS(self, sasRef):
        """
        Gives an object representing the specified SAS entry in this project.

        :param sasRef: A reference to the desired SASD entry. Either the SASD
          tree name (str), the SASD's index (int) or
          a SASD object (:class:`G2SmallAngle`)
        :returns: A :class:`G2SmallAngle` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.SASs`
            :class:`~G2PDF`
        """
        if isinstance(sasRef, G2SmallAngle):
            if sasRef.proj == self:
                return sasRef
            else:
                raise Exception(f"SASD {sasRef.name} not in current selected project")
        if sasRef in self._sasd():
            return G2SmallAngle(self.data[sasRef], sasRef, self)

        errmsg = ''
        try:
            # sasRef should be an index
            num = int(sasRef)
            sasRef = self._sasd()[num]
            return G2PDF(self.data[sasRef], sasRef, self)
        except ValueError:
            errmsg = f"sasRef {sasRef} not an object, name or SAS index in current selected project"
        except IndexError:
            errmsg = "sasRef {sasRef} out of range (max={len(self._sasd())-1}) in current selected project"
        if errmsg: raise Exception(errmsg)

    def SASs(self):
        """
        Returns a list of all the Small Angle histograms in the project.

        :returns: A list of :class:`G2SmallAngle` objects
        """
        return [G2SmallAngle(self.data[i],i,self) for i in self._sasd()]

    def add_SmallAngle(self, datafile):
        '''Placeholder for an eventual routine that will read a small angle dataset
        from a file.

        :param str datafile: The SASD data file to read, a filename.
        :returns: A :class:`G2SmallAngle` object for the SASD entry
        '''
        raise G2ScriptException("Error: add_SmallAngle not yet implemented.")
        #return G2SmallAngle(self.data[PDFname], PDFname, self)

    def get_Constraints(self,ctype):
        '''Returns a list of constraints of the type selected.

        :param str ctype: one of the following keywords: 'Hist', 'HAP', 'Phase', 'Global'
        :returns: a list of constraints, see the
          :ref:`constraint definition descriptions <Constraint_definitions_table>`. Note that
          if this list is changed (for example by deleting elements or by changing them)
          the constraints in the project are changed.
        '''
        if ctype in ('Hist', 'HAP', 'Phase', 'Global'):
            return self.data['Constraints']['data'][ctype]
        else:
            raise Exception(("get_Constraints error: value of ctype ({})"
                    +" must be 'Hist', 'HAP', 'Phase', or 'Global'.")
                                .format(ctype))

    def add_HoldConstr(self,varlist,reloadIdx=True,override=False):
        '''Set a hold constraint on a list of variables.

        Note that this will cause the project to be saved if not
        already done so. It will always save the .gpx file before
        creating constraint(s) if reloadIdx is True.

        :param list varlist: A list of variables to hold.
          Each value in the list may be one of the following three items:
          (A) a :class:`GSASIIobj.G2VarObj` object,
          (B) a variable name (str), or
          (C) a list/tuple of arguments for :meth:`make_var_obj`.
        :param bool reloadIdx: If True (default) the .gpx file will be
          saved and indexed prior to use. This is essential if atoms, phases
          or histograms have been added to the project.
        :param bool override: This routine looks up variables using
          :func:`GSASIIobj.getDescr` (which is not comprehensive). If
          not found, the routine will throw an exception, unless
          override=True is specified.

        Example::

            gpx.add_HoldConstr(('0::A4','0:1:D12',':0:Lam'))

        '''
        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()
        vardefwarn = False
        for var in varlist:
            # make var object
            if isinstance(var, str):
                var = self.make_var_obj(var,reloadIdx=False)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var,reloadIdx=False)
            if G2obj.getDescr(var.name) is None:
                vardefwarn = True
                print(f'Constraint var warning: No definition for variable {var.name}. Name correct?')
            # make constraint
            self.add_constraint_raw(_constr_type(var), [[1.0, var], None, None, 'h'])
        if vardefwarn and not override:
            raise Exception('Constraint var error: undefined variables.\n\tIf correct, use override=True in call.\n\tAlso, please let GSAS-II developers know so definition can be added.')

    def add_EquivConstr(self,varlist,multlist=[],reloadIdx=True,override=False):
        '''Set a equivalence on a list of variables.

        Note that this will cause the project to be saved if not
        already done so. It will always save the .gpx file before
        creating a constraint if reloadIdx is True.

        :param list varlist: A list of variables to make equivalent to the
          first item in the list.
          Each value in the list may be one of the following three items:
          (A) a :class:`GSASIIobj.G2VarObj` object,
          (B) a variable name (str), or
          (C) a list/tuple of arguments for :meth:`make_var_obj`.
        :param list multlist: a list of multipliers for each variable in
          varlist. If there are fewer values than supplied for varlist
          then missing values will be set to 1. The default is [] which
          means that all multipliers are 1.
        :param bool reloadIdx: If True (default) the .gpx file will be
          saved and indexed prior to use. This is essential if atoms, phases
          or histograms have been added to the project.
        :param bool override: This routine looks up variables using
          :func:`GSASIIobj.getDescr` (which is not comprehensive). If
          not found, the routine will throw an exception, unless
          override=True is specified.

        Examples::

            gpx.add_EquivConstr(('0::AUiso:0','0::AUiso:1','0::AUiso:2'))
            gpx.add_EquivConstr(('0::dAx:0','0::dAx:1'),[1,-1])

        '''
        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()
        if len(varlist) < 2:
            raise Exception('add_EquivConstr Error: varlist must have at least 2 variables')
        constr = []
        typ_prev = None
        vardefwarn = False
        for i,var in enumerate(varlist):
            m = 1.
            try:
                m = float(multlist[i])
            except IndexError:
                pass
            # make var object
            if isinstance(var, str):
                var = self.make_var_obj(var,reloadIdx=False)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var,reloadIdx=False)
            if G2obj.getDescr(var.name) is None:
                vardefwarn = True
                print(f'Constraint var warning: No definition for variable {var.name}. Name correct?')
            # make constraint
            constr.append([m,var])
            typ = _constr_type(var)
            if typ_prev is None:
                typ_prev = typ
                var_prev = var
            if typ_prev != typ:
                msg = 'Type ({}) for var {} is different from {} ({})'.format(typ,var,var_prev,typ_prev)
                raise Exception('add_EquivConstr Error: '+msg)
            typ_prev = typ
            var_prev = var
        if vardefwarn and not override:
            raise Exception('Constraint var error: undefined variables.\n\tIf correct, use override=True in call.\n\tAlso, please let GSAS-II developers know so definition can be added.')
        constr += [None, None, 'e']
        self.add_constraint_raw(typ, constr)

    def add_EqnConstr(self,total,varlist,multlist=[],reloadIdx=True,override=False):
        '''Set a constraint equation on a list of variables.

        Note that this will cause the project to be saved if not
        already done so. It will always save the .gpx file before
        creating a constraint if reloadIdx is True.

        :param float total: A value that the constraint must equal
        :param list varlist: A list of variables to use in the equation.
          Each value in the list may be one of the following three items:
          (A) a :class:`GSASIIobj.G2VarObj` object,
          (B) a variable name (str), or
          (C) a list/tuple of arguments for :meth:`make_var_obj`.
        :param list multlist: a list of multipliers for each variable in
          varlist. If there are fewer values than supplied for varlist
          then missing values will be set to 1. The default is [] which
          means that all multipliers are 1.
        :param bool reloadIdx: If True (default) the .gpx file will be
          saved and indexed prior to use. This is essential if atoms, phases
          or histograms have been added to the project.
        :param bool override: This routine looks up variables using
          :func:`GSASIIobj.getDescr` (which is not comprehensive). If
          not found, the routine will throw an exception, unless
          override=True is specified.

        Example::

            gpx.add_EqnConstr(1.0,('0::Ax:0','0::Ax:1'),[1,1])

        '''
        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()
        if len(varlist) < 2:
            raise Exception('Constraint var error: varlist must have at least 2 variables')
        try:
            float(total)
        except:
            raise Exception('Constraint var error: total be a valid float')
        constr = []
        typ_prev = None
        vardefwarn = False
        for i,var in enumerate(varlist):
            m = 1.
            try:
                m = float(multlist[i])
            except IndexError:
                pass
            # make var object
            if isinstance(var, str):
                var = self.make_var_obj(var,reloadIdx=False)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var,reloadIdx=False)
            if G2obj.getDescr(var.name) is None:
                vardefwarn = True
                print(f'Constraint var warning: No definition for variable {var.name}. Name correct?')
            # make constraint
            constr.append([m,var])
            typ = _constr_type(var)
            if typ_prev is None:
                typ_prev = typ
                var_prev = var
            if typ_prev != typ:
                msg = 'Type ({}) for var {} is different from {} ({})'.format(typ,var,var_prev,typ_prev)
                raise Exception('add_EquivConstr Error: '+msg)
            typ_prev = typ
            var_prev = var
        if vardefwarn and not override:
            raise Exception('Constraint var error: undefined variables.\n\tIf correct, use override=True in call.\n\tAlso, please let GSAS-II developers know so definition can be added.')
        constr += [float(total), None, 'c']
        self.add_constraint_raw(typ, constr)

    def add_NewVarConstr(self,varlist,multlist=[],name=None,vary=False,
                             reloadIdx=True,override=False):
        '''Set a new-variable constraint from a list of variables to
        create a new parameter from two or more predefined parameters.

        Note that this will cause the project to be saved, if not
        already done so. It will always save the .gpx file before
        creating a constraint if reloadIdx is True.

        :param list varlist: A list of variables to use in the expression.
          Each value in the list may be one of the following three items:
          (A) a :class:`GSASIIobj.G2VarObj` object,
          (B) a variable name (str), or
          (C) a list/tuple of arguments for :meth:`make_var_obj`.
        :param list multlist: a list of multipliers for each variable in
          varlist. If there are fewer values than supplied for varlist
          then missing values will be set to 1. The default is [] which
          means that all multipliers are 1.
        :param str name: An optional string to be supplied as a name for this
          new parameter.
        :param bool vary: Determines if the new variable should be flagged to
          be refined.
        :param bool reloadIdx: If True (default) the .gpx file will be
          saved and indexed prior to use. This is essential if atoms, phases
          or histograms have been added to the project.
        :param bool override: This routine looks up variables using
          :func:`GSASIIobj.getDescr` (which is not comprehensive). If
          not found, the routine will throw an exception, unless
          override=True is specified.

        Examples::

            gpx.add_NewVarConstr(('0::AFrac:0','0::AFrac:1'),[0.5,0.5],'avg',True)
            gpx.add_NewVarConstr(('0::AFrac:0','0::AFrac:1'),[1,-1],'diff',False,False)

        The example above is a way to treat two variables that are closely correlated.
        The first variable, labeled as avg, allows the two variables to refine in tandem
        while the second variable (diff) tracks their difference. In the initial stages of
        refinement only avg would be refined, but in the final stages, it might be possible
        to refine diff. The second False value in the second example prevents the
        .gpx file from being saved.
        '''
        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()
        if len(varlist) < 2:
            raise Exception('add_EquivEquation Error: varlist must have at least 2 variables')
        constr = []
        if name is not None:
            name = str(name)
        typ_prev = None
        vardefwarn = False
        for i,var in enumerate(varlist):
            m = 1.
            try:
                m = float(multlist[i])
            except IndexError:
                pass
            # make var object
            if isinstance(var, str):
                var = self.make_var_obj(var,reloadIdx=False)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var,reloadIdx=False)
            if G2obj.getDescr(var.name) is None:
                vardefwarn = True
                print(f'Constraint var warning: No definition for variable {var.name}. Name correct?')
            # make constraint
            constr.append([m,var])
            typ = _constr_type(var)
            if typ_prev is None:
                typ_prev = typ
                var_prev = var
            if typ_prev != typ:
                msg = 'Type ({}) for var {} is different from {} ({})'.format(typ,var,var_prev,typ_prev)
                raise Exception('add_EquivConstr Error: '+msg)
            typ_prev = typ
            var_prev = var
        if vardefwarn and not override:
            raise Exception('Constraint var error: undefined variables.\n\tIf correct, use override=True in call.\n\tAlso, please let GSAS-II developers know so definition can be added.')
        constr += [name, bool(vary), 'f']
        self.add_constraint_raw(typ, constr)

    def add_constraint_raw(self, cons_scope, constr):
        """Adds a constraint to the project.

        :param str cons_scope: should be one of "Hist", "Phase", "HAP", or "Global".
        :param list constr: a constraint coded with  :class:`GSASIIobj.G2VarObj`
          objects as described in the
          :ref:`constraint definition descriptions <Constraint_definitions_table>`.

        WARNING this function does not check the constraint is well-constructed.
        Please use :meth:`G2Project.add_HoldConstr` or
        :meth:`G2Project.add_EquivConstr` (etc.) instead, unless you are really
        certain you know what you are doing.
        """
        constrs = self.data['Constraints']['data']
        if 'Global' not in constrs:
            constrs['Global'] = []
        constrs[cons_scope].append(constr)

    def hold_many(self, vars, ctype):
        """Apply holds for all the variables in vars, for constraint of a given type.
        This routine has been superceeded by :meth:`add_Hold`

        :param list vars: A list of variables to hold. Each may be a
          :class:`GSASIIobj.G2VarObj` object, a variable name (str), or a
          list/tuple of arguments for :meth:`make_var_obj`.
        :param str ctype: A string constraint type specifier, passed directly to
          :meth:`add_constraint_raw` as consType. Should be one of "Hist", "Phase",
          or "HAP" ("Global" not implemented).
        """
        G2fil.G2Print('G2Phase.hold_many Warning: replace calls to hold_many() with add_Hold()')
        for var in vars:
            if isinstance(var, str):
                var = self.make_var_obj(var)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var)
            self.add_constraint_raw(ctype, [[1.0, var], None, None, 'h'])

    def make_var_obj(self, phase=None, hist=None, varname=None, atomId=None,
                     reloadIdx=True):
        """Wrapper to create a G2VarObj. Takes either a string representation ("p:h:name:a")
        or individual names of phase, histogram, varname, and atomId.

        Automatically converts string phase, hist, or atom names into the ID required
        by G2VarObj.

        Note that this will cause the project to be saved if not
        already done so.
        """

        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()

        # If string representation, short circuit
        if hist is None and varname is None and atomId is None:
            if isinstance(phase, str) and ':' in phase:
                return G2obj.G2VarObj(phase)

        # Get phase index
        phaseObj = None
        if isinstance(phase, G2Phase):
            phaseObj = phase
            phase = G2obj.PhaseRanIdLookup[phase.ranId]
        elif phase in self.data['Phases']:
            phaseObj = self.phase(phase)
            phaseRanId = phaseObj.ranId
            phase = G2obj.PhaseRanIdLookup[phaseRanId]

        # Get histogram index
        if isinstance(hist, G2PwdrData):
            hist = G2obj.HistRanIdLookup[hist.ranId]
        elif hist in self.data:
            histRanId = self.histogram(hist).ranId
            hist = G2obj.HistRanIdLookup[histRanId]

        # Get atom index (if any)
        if isinstance(atomId, G2AtomRecord):
            atomId = G2obj.AtomRanIdLookup[phase][atomId.ranId]
        elif phaseObj:
            atomObj = phaseObj.atom(atomId)
            if atomObj:
                atomRanId = atomObj.ranId
                atomId = G2obj.AtomRanIdLookup[phase][atomRanId]

        return G2obj.G2VarObj(phase, hist, varname, atomId)

    def add_image(self, imagefile, fmthint=None, defaultImage=None,
                      indexList=None, cacheImage=False,
                      URL=False, download_loc=None):
        """Load an image into a project

        :param str imagefile: The image file to read, a filename.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu).
        :param str defaultImage: The name of an image to use as a default for
          setting parameters for the image file to read.
        :param list indexList: specifies the image numbers (counting from zero)
          to be used from the file when a file has multiple images. A value of
          ``[0,2,3]`` will cause the only first, third and fourth images in the file
          to be included in the project.
        :param bool cacheImage: When True, the image is cached to save
          in rereading it later. Default is False (no caching).
        :param bool URL: if True, the contents of imagefile is a URL 
          and the file will be downloaded and saved. The file will be 
          written in the specified directory (see `download_loc`)
          or a temporary location, if not specified. Note that 
          if a temporary location, if the proiject (.gpx) file is 
          saved, the image may not be accessible if the .gpx file
          is later reopened. Default is False.
          If URL is specified and the Python requests package is 
          not installed, a `ModuleNotFoundError` Exception will occur. 
          will occur. 
        :param str download_loc: a location or file name where the 
          image will be saved. Note that for almost all image types, 
          the image cannot be read if the file extension does not
          match what is expected for the format. (This can be determined
          by looking at the importer code; if `strictExtension=True`, 
          the extension must be in the `extensionlist` list.)
          If only a directory is specified, the file name will be taken 
          from the URL, which will likely cause problems if it does
          not match the needed extension. 
          If URL is specified and the default download_loc
          value is used (None), the image will be saved in a temporary 
          location that will persist until the OS removes it.
        :returns: a list of :class:`G2Image` object(s) for the added image(s)
        """
        LoadG2fil()
        if not URL: imagefile = os.path.abspath(os.path.expanduser(imagefile))
        readers = import_generic(imagefile, Readers['Image'],
                    fmthint=fmthint, URL=URL, download_loc=download_loc)
        objlist = []
        for i,rd in enumerate(readers):
            if indexList is not None and i not in indexList:
                G2fil.G2Print("Image {} skipped".format(i))
                continue
            if rd.SciPy:        #was default read by scipy; needs 1 time fixes
                G2fil.G2Print('Warning: Image {} read by generic SciPy import. Image parameters likely wrong'.format(imagefile))
                #see this: G2IO.EditImageParms(self,rd.Data,rd.Comments,rd.Image,imagefile)
                rd.SciPy = False
            rd.readfilename = imagefile
            if rd.repeatcount == 1 and not rd.repeat: # skip image number if only one in set
                rd.Data['ImageTag'] = None
            else:
                rd.Data['ImageTag'] = rd.repeatcount
            rd.Data['formatName'] = rd.formatName
            if rd.sumfile:
                rd.readfilename = rd.sumfile
            # Load generic metadata, as configured
            G2fil.GetColumnMetadata(rd)
            # Code from G2IO.LoadImage2Tree(rd.readfilename,self,rd.Comments,rd.Data,rd.Npix,rd.Image)
            Imax = np.amax(rd.Image)
            ImgNames = [i[0] for i in self.names if i[0].startswith('IMG ')]
            TreeLbl = 'IMG '+os.path.basename(imagefile)
            ImageTag = rd.Data.get('ImageTag')
            if ImageTag:
                TreeLbl += ' #'+'%04d'%(ImageTag)
                imageInfo = (imagefile,ImageTag)
            else:
                imageInfo = imagefile
            TreeName = G2obj.MakeUniqueLabel(TreeLbl,ImgNames)
            # MT dict to contain image info
            ImgDict = {}
            ImgDict['data'] = [rd.Npix,imageInfo]
            ImgDict['Comments'] = rd.Comments
            if defaultImage:
                if isinstance(defaultImage, G2Image):
                    if defaultImage.proj == self:
                        defaultImage = self.data[defaultImage.name]['data']
                    else:
                        raise Exception("Image {} not in current selected project".format(defaultImage.name))
                elif defaultImage in self._images():
                    defaultImage = self.data[defaultImage]['data']
                else:
                    defaultImage = None
            Data = rd.Data
            if defaultImage:
                Data = copy.copy(defaultImage)
                Data['showLines'] = True
                Data['ring'] = []
                Data['rings'] = []
                Data['cutoff'] = 10.
                Data['pixLimit'] = 20
                Data['edgemin'] = 100000000
                Data['calibdmin'] = 0.5
                Data['calibskip'] = 0
                Data['ellipses'] = []
                Data['calibrant'] = ''
                Data['GonioAngles'] = [0.,0.,0.]
                Data['DetDepthRef'] = False
            else:
                Data['type'] = 'PWDR'
                Data['color'] = GSASIIpath.GetConfigValue('Contour_color','Paired')
                if 'tilt' not in Data:          #defaults if not preset in e.g. Bruker importer
                    Data['tilt'] = 0.0
                    Data['rotation'] = 0.0
                    Data['pixLimit'] = 20
                    Data['calibdmin'] = 0.5
                    Data['cutoff'] = 10.
                Data['showLines'] = False
                Data['calibskip'] = 0
                Data['ring'] = []
                Data['rings'] = []
                Data['edgemin'] = 100000000
                Data['ellipses'] = []
                Data['GonioAngles'] = [0.,0.,0.]
                Data['DetDepth'] = 0.
                Data['DetDepthRef'] = False
                Data['calibrant'] = ''
                Data['IOtth'] = [5.0,50.0]
                Data['LRazimuth'] = [0.,180.]
                Data['azmthOff'] = 0.0
                Data['outChannels'] = 2500
                Data['outAzimuths'] = 1
                Data['centerAzm'] = False
                Data['fullIntegrate'] = GSASIIpath.GetConfigValue('fullIntegrate',True)
                Data['setRings'] = False
                Data['background image'] = ['',-1.0]
                Data['dark image'] = ['',-1.0]
                Data['Flat Bkg'] = 0.0
                Data['Oblique'] = [0.5,False]
            Data['varyList'] = {'dist':True,'det-X':True,'det-Y':True,'tilt':True,'phi':True,'dep':False,'wave':False}
            Data['setDefault'] = False
            Data['range'] = [(0,Imax),[0,Imax]]
            ImgDict['Image Controls'] = Data
            ImgDict['Masks'] = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],
                'Frames':[],'Thresholds':[(0,Imax),[0,Imax]],'SpotMask':{'esdMul':3.,'spotMask':None},}
            ImgDict['Stress/Strain']  = {'Type':'True','d-zero':[],'Sample phi':0.0,
                'Sample z':0.0,'Sample load':0.0}
            self.names.append([TreeName]+['Comments','Image Controls','Masks','Stress/Strain'])
            self.data[TreeName] = ImgDict
            if cacheImage:
                objlist.append(G2Image(self.data[TreeName], TreeName, self, image=rd.Image))
            else:
                objlist.append(G2Image(self.data[TreeName], TreeName, self))
                del rd.Image
        return objlist

    def imageMultiDistCalib(self,imageList=None,verbose=False):
        '''Invokes a global calibration fit (same as Image Controls/Calibration/Multi-distance Recalibrate
        menu command) with images as multiple distance settings.
        Note that for this to work properly, the initial calibration parameters
        (center, wavelength, distance & tilts) must be close enough to converge.
        This may produce a better result if run more than once.

        See :ref:`MultiDist_Example` for example code.

        :param str imageList: the images to include in the fit, if not specified
          all images in the project will be included.

        :returns: parmDict,covData where parmDict has the refined parameters
          and their values and covData is a dict containing the covariance matrix ('covMatrix'),
          the number of ring picks ('obs') the reduced Chi-squared ('chisq'),
          the names of the variables ('varyList') and their values ('variables')
        '''
        if imageList is None:
            imageList = self.images()

        # code based on GSASIIimgGUI..OnDistRecalib
        obsArr = np.array([]).reshape(0,4)
        parmDict = {}
        varList = []
        HKL = {}

        for img in imageList:
            name = img.name
            G2fil.G2Print ('getting rings for',name)
            Data = img.data['Image Controls']
            key = str(int(Data['setdist']))
            # create a parameter dict for combined fit
            if 'wavelength' not in parmDict:
                parmDict['wavelength'] = Data['wavelength']
                if Data['varyList']['wave']:
                    varList += ['wavelength']
                if Data['varyList']['dist']:
                    raise Exception(
                                'You cannot vary individual detector positions and the global wavelength.\n\nChange flags for 1st image.',
                                'Conflicting vars')
                parmDict['dep'] = Data['DetDepth']
                if Data['varyList']['dep']:
                    varList += ['dep']
                # distance flag determines if individual values are refined
                if not Data['varyList']['dist']:
                    # starts as zero, single variable, always refined
                    parmDict['deltaDist'] = 0.
                    varList += ['deltaDist']
                parmDict['phi'] = Data['rotation']
                if Data['varyList']['phi']:
                    varList += ['phi']
                parmDict['tilt'] = Data['tilt']
                if Data['varyList']['tilt']:
                    varList += ['tilt']

            ImageZ = _getCorrImage(Readers['Image'],self,img)
            Data['setRings'] = True
            Masks = img.data['Masks']
            result = G2img.ImageRecalibrate(None,ImageZ,Data,Masks,getRingsOnly=True)
            if not len(result):
                raise Exception('calibrant missing from local image calibrants files')
            rings,HKL[key] = result
            # add detector set dist into data array, create a single really large array
            distarr = np.zeros_like(rings[:,2:3])
            if 'setdist' not in Data:
                raise Exception('Distance (setdist) not in image metadata')
            distarr += Data['setdist']
            obsArr = np.concatenate((
                        obsArr,
                        np.concatenate((rings[:,0:2],distarr,rings[:,2:3]),axis=1)),axis=0)
            if 'deltaDist' not in parmDict:
                # starts as zero, variable refined for each image
                parmDict['delta'+key] = 0
                varList += ['delta'+key]
            for i,z in enumerate(['X','Y']):
                v = 'det-'+z
                if v+key in parmDict:
                    raise Exception('Error: two images with setdist ~=',key)
                parmDict[v+key] = Data['center'][i]
                if Data['varyList'][v]:
                    varList += [v+key]
        #GSASIIpath.IPyBreak()
        G2fil.G2Print('\nFitting',obsArr.shape[0],'ring picks and',len(varList),'variables...')
        result = G2img.FitMultiDist(obsArr,varList,parmDict,covar=True,Print=verbose)

        for img in imageList: # update GPX info with fit results
            name = img.name
            #print ('updating',name)
            Data = img.data['Image Controls']
            Data['wavelength'] = parmDict['wavelength']
            key = str(int(Data['setdist']))
            Data['center'] = [parmDict['det-X'+key],parmDict['det-Y'+key]]
            if 'deltaDist' in parmDict:
                Data['distance'] = Data['setdist'] - parmDict['deltaDist']
            else:
                Data['distance'] = Data['setdist'] - parmDict['delta'+key]
            Data['rotation'] = np.mod(parmDict['phi'],360.0)
            Data['tilt'] = parmDict['tilt']
            Data['DetDepth'] = parmDict['dep']
            N = len(Data['ellipses'])
            Data['ellipses'] = []           #clear away individual ellipse fits
            for H in HKL[key][:N]:
                ellipse = G2img.GetEllipse(H[3],Data)
                Data['ellipses'].append(copy.deepcopy(ellipse+('b',)))

        covData = {'title':'Multi-distance recalibrate','covMatrix':result[3],
                       'obs':obsArr.shape[0],'chisq':result[0],
                       'varyList':varList,'variables':result[1]}
        return parmDict,covData

    def set_Frozen(self, variable=None, histogram=None, mode='remove'):
        '''Removes one or more Frozen variables (or adds one)
        (See :ref:`Parameter Limits<ParameterLimits>` description.)
        Note that use of this
        will cause the project to be saved if not already done so.

        :param str variable: a variable name as a str or
          (as a :class:`GSASIIobj.G2VarObj` object). Should
          not contain wildcards.
          If None (default), all frozen variables are deleted
          from the project, unless a sequential fit and
          a histogram is specified.
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
          Used for sequential fits only.
        :param str mode: The default mode is to remove variables
          from the appropriate Frozen list, but if the mode
          is specified as 'add', the variable is added to the
          list.
        :returns: True if the variable was added or removed, False
          otherwise. Exceptions are generated with invalid requests.
        '''
        Controls = self.data['Controls']['data']
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        if G2obj.TestIndexAll(): self.index_ids()
        G2obj.patchControls(Controls)

        if mode == 'remove':
            if variable is None and histogram is None:
                Controls['parmFrozen'] = {}
                return True
        elif mode == 'add':
            if variable is None:
                raise Exception('set_Frozen error: variable must be specified')
        else:
            raise Exception('Undefined mode ({}) in set_Frozen'.format(mode))

        if histogram is None:
            h = 'FrozenList'
        else:
            hist = self.histogram(histogram)
            if hist:
                h = hist.name
                if h not in Controls['parmFrozen']: # no such list, not found
                    return False
            else:
                raise Exception('set_Frozen error: histogram {} not found'.format(histogram))
        if mode == 'remove':
            if variable is None:
                Controls['parmFrozen'][h] = []
                return True
            if h not in Controls['parmFrozen']:
                return True
            delList = []
            for i,v in enumerate(Controls['parmFrozen'][h]):
                if v == variable: delList.append(i)
            if delList:
                for i in reversed(delList):
                    del Controls['parmFrozen'][h][i]
                return True
            return False
        elif mode == 'add':
            if type(variable) is str:
                variable = G2obj.G2VarObj(variable)
            elif type(v) is not G2obj.G2VarObj:
                raise Exception(
                    'set_Frozen error: variable {} wrong type ({})'
                    .format(variable,type(variable)))
            if h not in Controls['parmFrozen']: Controls['parmFrozen'][h] = []
            Controls['parmFrozen'][h].append(variable)
            return True

    def get_Frozen(self, histogram=None):
        '''Gets a list of Frozen variables.
        (See :ref:`Parameter Limits<ParameterLimits>` description.)
        Note that use of this
        will cause the project to be saved if not already done so.

        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number. Used
          for sequential fits only. If left as the default (None)
          for a sequential fit, all Frozen variables in all
          histograms are returned.
        :returns: a list containing variable names, as str values
        '''
        Controls = self.data['Controls']['data']
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        if G2obj.TestIndexAll(): self.index_ids()
        G2obj.patchControls(Controls)
        if histogram:
            hist = self.histogram(histogram)
            if hist:
                h = hist.name
                return [str(i) for i in Controls['parmFrozen'][h]]
            elif histogram.lower() == 'all':
                names = set()
                for h in self.histograms():
                    if h.name not in Controls['parmFrozen']: continue
                    names = names.union([str(i) for i in Controls['parmFrozen'][h.name]])
                return list(names)
            else:
                raise Exception('histogram {} not recognized'.format(histogram))
        elif 'FrozenList' in Controls['parmFrozen']:
            return [str(i) for i in Controls['parmFrozen']['FrozenList']]
        else:
            return []

    def get_Controls(self, control, variable=None):
        '''Return project controls settings

        :param str control: the item to be returned. See below for allowed values.
        :param str variable: a variable name as a str or
          (as a :class:`GSASIIobj.G2VarObj` object).
          Used only with control set to "parmMin" or "parmMax".
        :returns: The value for the control.

        Allowed values for parameter control:

            * cycles: the maximum number of cycles (returns int)
            * sequential: the histograms used for a sequential refinement as a list
              of histogram names or an empty list when in non-sequential mode.
            * Reverse Seq: returns True or False. True indicates that fitting of the
              sequence of histograms proceeds in reversed order.
            * seqCopy: returns True or False. True indicates that results from
              each sequential fit are used as the starting point for the next
              histogram.
            * parmMin & parmMax: retrieves a maximum or minimum value for
              a refined parameter. Note that variable will be a GSAS-II
              variable name, optionally with * specified for a histogram
              or atom number. Return value will be a float.
              (See :ref:`Parameter Limits<ParameterLimits>` description.)
            * Anything else returns the value in the Controls dict, if present. An
              exception is raised if the control value is not present.

        .. seealso::
            :meth:`set_Controls`
        '''
        if control.startswith('parmM'):
            if not variable:
                raise Exception('set_Controls requires a value for variable for control=parmMin or parmMax')
            for key in ('parmMinDict','parmMaxDict','parmFrozen'):
                if key not in self.data['Controls']['data']: self.data['Controls']['data'][key] = {}
            if G2obj.TestIndexAll(): self.index_ids()
            G2obj.patchControls(self.data['Controls']['data'])
        if control == 'cycles':
            return self.data['Controls']['data']['max cyc']
        elif control == 'sequential':
            return self.data['Controls']['data']['Seq Data']
        elif control == 'Reverse Seq':
            return self.data['Controls']['data']['Reverse Seq']
        elif control == 'seqCopy':
            return self.data['Controls']['data']['Copy2Next']
        elif control == 'parmMin' or control == 'parmMax':
            key = G2obj.G2VarObj(variable)
            return G2obj.prmLookup(
                variable,
                self.data['Controls']['data'][control+'Dict'])[1]
        elif control in self.data['Controls']['data']:
            return self.data['Controls']['data'][control]
        else:
            G2fil.G2Print('Defined Controls:',self.data['Controls']['data'].keys())
            raise Exception('{} is an invalid control value'.format(control))

    def set_Controls(self, control, value, variable=None):
        '''Set project controls.

        Note that use of this with control set to parmMin or parmMax
        will cause the project to be saved if not already done so.

        :param str control: the item to be set. See below for allowed values.
        :param value: the value to be set.
        :param str variable: used only with control set to "parmMin" or "parmMax"

        Allowed values for *control* parameter:

        * ``'cycles'``: sets the maximum number of cycles (value must be int)
        * ``'sequential'``: sets the histograms to be used for a sequential refinement.
          Use an empty list to turn off sequential fitting.
          The values in the list may be the name of the histogram (a str), or
          a ranId or index (int values), see :meth:`histogram`.
        * ``'seqCopy'``: when True, the results from each sequential fit are used as
          the starting point for the next. After each fit is is set to False.
          Ignored for non-sequential fits.
        * ``'Reverse Seq'``: when True, sequential refinement is performed on the
          reversed list of histograms.
        * ``'parmMin'`` & ``'parmMax'``: set a maximum or minimum value for a refined
          parameter. Note that variable will be a GSAS-II variable name,
          optionally with * specified for a histogram or atom number and
          value must be a float.
          (See :ref:`Parameter Limits<ParameterLimits>` description.)

        .. seealso::
            :meth:`get_Controls`
        '''
        if control.startswith('parmM'):
            if not variable:
                raise Exception('set_Controls requires a value for variable for control=parmMin or parmMax')
            for key in ('parmMinDict','parmMaxDict','parmFrozen'):
                if key not in self.data['Controls']['data']: self.data['Controls']['data'][key] = {}
            if G2obj.TestIndexAll(): self.index_ids()
            G2obj.patchControls(self.data['Controls']['data'])
        if control == 'cycles':
            self.data['Controls']['data']['max cyc'] = int(value)
        elif control == 'seqCopy':
            self.data['Controls']['data']['Copy2Next'] = bool(value)
        elif control == 'Reverse Seq':
            self.data['Controls']['data']['Reverse Seq'] = bool(value)
        elif control == 'sequential':
            histlist = []
            for i,j in enumerate(value):
                h = self.histogram(j)
                if h:
                    histlist.append(h.name)
                else:
                    raise Exception('item #{} ({}) is an invalid histogram value'
                                        .format(i,j))
            self.data['Controls']['data']['Seq Data'] = histlist
        elif control == 'parmMin' or control == 'parmMax':
            key = G2obj.G2VarObj(variable)
            self.data['Controls']['data'][control+'Dict'][key] = float(value)
        else:
            raise Exception('{} is an invalid control value'.format(control))

    def copyHistParms(self,sourcehist,targethistlist='all',modelist='all'):
        '''Copy histogram information from one histogram to others

        :param sourcehist: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`),
            a histogram name or the index number of the histogram.
            if the string 'all' (default value), then all histograms in
            the project are used.
        :param list modelist: May be a list of sections to copy, which may
           include 'Background', 'Instrument Parameters', 'Limits' and
           'Sample Parameters' (items may be shortened to uniqueness and
           capitalization is ignored, so ['b','i','L','s'] will work.)
           The default value, 'all' causes the listed sections to

        '''
        sections = ('Background','Instrument Parameters','Limits',
                        'Sample Parameters')
        hist_in = self.histogram(sourcehist)
        if not hist_in:
            raise Exception('{} is not a valid histogram'.format(sourcehist))
        if targethistlist == "all":
            targethistlist = self.histograms()
        if 'all' in modelist:
            copysections = sections
        else:
            copysections = set()
            for s in sections:
                for m in modelist:
                    if s.lower().startswith(m.lower()):
                        copysections.add(s)
        for h in targethistlist:
            hist_out = self.histogram(h)
            if not hist_out:
                raise Exception('{} is not a valid histogram'.format(h))
            for key in copysections:
                hist_out[key] = copy.deepcopy(hist_in[key])

    def get_VaryList(self):
        '''Returns a list of the refined variables in the
        last refinement cycle

        :returns: a list of variables or None if no refinement has been
          performed.
        '''
        try:
            return self['Covariance']['data']['varyList']
        except:
            return

    def get_ParmList(self):
        '''Returns a list of all the parameters defined in the
        last refinement cycle

        :returns: a list of parameters or None if no refinement has been
          performed.
        '''
        if 'parmDict' not in self['Covariance']['data']:
            raise G2ScriptException('No parameters found in project, has a refinement been run?')
        try:
            return list(self['Covariance']['data']['parmDict'].keys())
        except:
            return

    def get_Variable(self,var):
        '''Returns the value and standard uncertainty (esd) for a variable
        parameters, as defined in the last refinement cycle

        :param str var: a variable name of form '<p>:<h>:<name>', such as
          ':0:Scale'
        :returns: (value,esd) if the parameter is refined or
          (value, None) if the variable is in a constraint or is not
          refined or None if the parameter is not found.
        '''
        if 'parmDict' not in self['Covariance']['data']:
            raise G2ScriptException('No parameters found in project, has a refinement been run?')
        if var not in self['Covariance']['data']['parmDict']:
            return None
        val = self['Covariance']['data']['parmDict'][var]
        try:
            pos = self['Covariance']['data']['varyList'].index(var)
            esd = np.sqrt(self['Covariance']['data']['covMatrix'][pos,pos])
            return (val,esd)
        except ValueError:
            return (val,None)

    def get_Covariance(self,varList):
        '''Returns the values and covariance matrix for a series of variable
        parameters. as defined in the last refinement cycle

        :param tuple varList: a list of variable names of form '<p>:<h>:<name>'
        :returns: (valueList,CovMatrix) where valueList contains the (n) values
          in the same order as varList (also length n) and CovMatrix is a
          (n x n) matrix. If any variable name is not found in the varyList
          then None is returned.

        Use this code, where sig provides standard uncertainties for
        parameters and where covArray provides the correlation between
        off-diagonal terms::

            sig = np.sqrt(np.diag(covMatrix))
            xvar = np.outer(sig,np.ones_like(sig))
            covArray = np.divide(np.divide(covMatrix,xvar),xvar.T)

        '''
        missing = [i for i in varList if i not in self['Covariance']['data']['varyList']]
        if missing:
            G2fil.G2Print('Warning: Variable(s) {} were not found in the varyList'.format(missing))
            return None
        if 'parmDict' not in self['Covariance']['data']:
            raise G2ScriptException('No parameters found in project, has a refinement been run?')
        vals = [self['Covariance']['data']['parmDict'][i] for i in varList]
        cov = G2mth.getVCov(varList,
                            self['Covariance']['data']['varyList'],
                            self['Covariance']['data']['covMatrix'])
        return (vals,cov)

    def ComputeWorstFit(self):
        '''Computes the worst-fit parameters in a model.

        :returns: (keys, derivCalcs, varyList) where:

          * keys is a list of parameter names
            where the names are ordered such that first entry in the list
            will produce the largest change in the fit if refined and the last
            entry will have the smallest change;

          * derivCalcs is a dict where the key is a variable name and the
            value is a list with three partial derivative values for
            d(Chi**2)/d(var) where the derivatives are computed
            for values v-d to v; v-d to v+d; v to v+d where v is
            the current value for the variable and d is a small delta
            value chosen for that variable type;

          * varyList is a list of the parameters that are currently set to
            be varied.
        '''
        self.save()
        derivCalcs,varyList = G2strMain.Refine(self.filename,None,allDerivs=True)
        keys = sorted(derivCalcs,key=lambda x:abs(derivCalcs[x][1]),reverse=True)
        return (keys,derivCalcs,varyList)

class G2AtomRecord(G2ObjectWrapper):
    """Wrapper for an atom record. Allows many atom properties to be access
    and changed. See the :ref:`Atom Records description <Atoms_table>`
    for the details on what information is contained in an atom record.

    Scripts should not try to create a :class:`G2AtomRecord` object directly as
    these objects are created via access from a :class:`G2Phase` object.

    Example showing some uses of :class:`G2AtomRecord` methods:

    >>> atom = some_phase.atom("O3")
    >>> # We can access the underlying data structure (a list):
    >>> atom.data
    ['O3', 'O-2', '', ... ]
    >>> # We can also use wrapper accessors to get or change atom info:
    >>> atom.coordinates
    (0.33, 0.15, 0.5)
    >>> atom.coordinates = [1/3, .1, 1/2]
    >>> atom.coordinates
    (0.3333333333333333, 0.1, 0.5)
    >>> atom.refinement_flags
    'FX'
    >>> atom.ranId
    4615973324315876477
    >>> atom.occupancy
    1.0
    """
    def __init__(self, data, indices, proj):
        self.data = data
        self.cx, self.ct, self.cs, self.cia = indices
        self.proj = proj

#    @property
#    def X(self):
#        '''Get or set the associated atom's X.
#        Use as ``x = atom.X`` to obtain the value and
#        ``atom.X = x`` to set the value.
#        '''
#        pass
#    @X.setter
#    def X(self, val):
#        pass

    @property
    def label(self):
        '''Get the associated atom's label.
        Use as ``x = atom.label`` to obtain the value and
        ``atom.label = x`` to set the value.
        '''
        return self.data[self.ct-1]
    @label.setter
    def label(self, val):
        self.data[self.ct-1] = str(val)

    @property
    def type(self):
        '''Get or set the associated atom's type. Call as a variable
        (``x = atom.type``) to obtain the value or use
        ``atom.type = x`` to change the type. It is the user's
        responsibility to make sure that the atom type is valid;
        no checking is done here.

        .. seealso::
            :meth:`element`
        '''
        return self.data[self.ct]
    @type.setter
    def type(self, val):
        # TODO: should check if atom type is defined
        self.data[self.ct] = str(val)

    @property
    def element(self):
        '''Parses element symbol from the atom type symbol for the atom
        associated with the current object.

        .. seealso::
            :meth:`type`
        '''
        import re
        try:
            return re.match('^([A-Z][a-z]?)',self.data[self.ct]).group(1)
        except:
            raise Exception("element parse error with type {}".
                                format(self.data[self.ct]))

    @property
    def refinement_flags(self):
        '''Get or set refinement flags for the associated atom.
        Use as ``x = atom.refinement_flags`` to obtain the flags and
        ``atom.refinement_flags = "XU"`` (etc) to set the value.
        '''
        return self.data[self.ct+1]
    @refinement_flags.setter
    def refinement_flags(self, other):
        # Automatically check it is a valid refinement
        for c in other:
            if c not in ' FXU':
                raise ValueError("Invalid atom refinement: ", other)
        self.data[self.ct+1] = other

    @property
    def coordinates(self):
        '''Get or set the associated atom's coordinates.
        Use as ``x = atom.coordinates`` to obtain a tuple with
        the three (x,y,z) values and ``atom.coordinates = (x,y,z)``
        to set the values.

        Changes needed to adapt for changes in site symmetry have not yet been
        implemented:
        '''
        return tuple(self.data[self.cx:self.cx+3])
    @coordinates.setter
    def coordinates(self, val):
        if len(val) != 3:
            raise ValueError(f"coordinates are of wrong length {val}")
        try:
            self.data[self.cx:self.cx+3] = [float(i) for i in val]
        except:
            raise ValueError(f"conversion error with coordinates {val}")
        # TODO: should recompute the atom site symmetries here

    @property
    def occupancy(self):
        '''Get or set the associated atom's site fraction.
        Use as ``x = atom.occupancy`` to obtain the value and
        ``atom.occupancy = x`` to set the value.
        '''
        return self.data[self.cx+3]
    @occupancy.setter
    def occupancy(self, val):
        self.data[self.cx+3] = float(val)

    @property
    def mult(self):
        '''Get the associated atom's multiplicity value. Should not be
        changed by user.
        '''
        return self.data[self.cs+1]

    @property
    def ranId(self):
        '''Get the associated atom's Random Id number. Don't change this.
        '''
        return self.data[self.cia+8]

    @property
    def adp_flag(self):
        '''Get the associated atom's iso/aniso setting. The value
        will be 'I' or 'A'. No API provision is offered to change
        this.
        '''
        # Either 'I' or 'A'
        return self.data[self.cia]

    @property
    def ADP(self):
        '''Get or set the associated atom's Uiso or Uaniso value(s).
        Use as ``x = atom.ADP`` to obtain the value(s) and
        ``atom.ADP = x`` to set the value(s). For isotropic atoms
        a single float value is returned (or used to set). For
        anisotropic atoms a list of six values is used.

        .. seealso::
            :meth:`adp_flag`
            :meth:`uiso`
        '''
        if self.adp_flag == 'I':
            return self.data[self.cia+1]
        else:
            return self.data[self.cia+2:self.cia+8]
    @ADP.setter
    def ADP(self, value):
        if self.adp_flag == 'I':
            self.data[self.cia+1] = float(value)
        else:
            assert len(value) == 6
            self.data[self.cia+2:self.cia+8] = [float(v) for v in value]

    @property
    def uiso(self):
        '''A synonym for :meth:`ADP` to be used for Isotropic atoms.
        Get or set the associated atom's Uiso value.
        Use as ``x = atom.uiso`` to obtain the value and
        ``atom.uiso = x`` to set the value. A
        single float value is returned or used to set.

        .. seealso::
            :meth:`adp_flag`
            :meth:`ADP`
        '''
        if self.adp_flag != 'I':
            raise G2ScriptException(f"Atom {self.label} is not isotropic")
        return self.ADP
    @uiso.setter
    def uiso(self, value):
        if self.adp_flag != 'I':
            raise G2ScriptException(f"Atom {self.label} is not isotropic")
        self.ADP = value

class G2PwdrData(G2ObjectWrapper):
    """Wraps a Powder Data Histogram.
    The object contains these class variables:

        * G2PwdrData.proj: contains a reference to the :class:`G2Project`
          object that contains this histogram
        * G2PwdrData.name: contains the name of the histogram
        * G2PwdrData.data: contains the histogram's associated data in a dict,
          as documented for the :ref:`Powder Diffraction Tree<Powder_table>`.
          The actual histogram values are contained in the 'data' dict item,
          as documented for Data.

    Scripts should not try to create a :class:`G2PwdrData` object directly as
    :meth:`G2PwdrData.__init__` should be invoked from inside :class:`G2Project`.

    """
    def __init__(self, data, proj, name):
        self.data = data
        self.proj = proj
        self.name = name

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ['Limits', 'Sample Parameters', 'Background',
                      'Instrument Parameters']
        return key in valid_keys

    #@property
    #def name(self):
    #    return self.data['data'][-1]

    @property
    def ranId(self):
        return self.data['data'][0]['ranId']

    @property
    def residuals(self):
        '''Provides a dictionary with with the R-factors for this histogram.
        Includes the weighted and unweighted profile terms (R, Rb, wR, wRb, wRmin)
        as well as the Bragg R-values for each phase (ph:H:Rf and ph:H:Rf^2).
        '''
        data = self.data['data'][0]
        return {key: data[key] for key in data
                if key in ['R', 'Rb', 'wR', 'wRb', 'wRmin']
                   or ':' in key}

    @property
    def InstrumentParameters(self):
        '''Provides a dictionary with with the Instrument Parameters
        for this histogram.
        '''
        return self.data['Instrument Parameters'][0]

    @property
    def SampleParameters(self):
        '''Provides a dictionary with with the Sample Parameters
        for this histogram.
        '''
        return self.data['Sample Parameters']

    @property
    def Background(self):
        '''Provides a list with with the Background parameters
        for this histogram.

        :returns: list containing a list and dict with background values
        '''
        return self.data['Background']

    def add_back_peak(self,pos,int,sig,gam,refflags=[]):
        '''Adds a background peak to the Background parameters

        :param float pos: position of peak, a 2theta or TOF value
        :param float int: integrated intensity of background peak, usually large
        :param float sig: Gaussian width of background peak, usually large
        :param float gam: Lorentzian width of background peak, usually unused (small)
        :param list refflags: a list of 1 to 4 boolean refinement flags for
            pos,int,sig & gam, respectively (use [0,1] to refine int only).
            Defaults to [] which means nothing is refined.
        '''
        if 'peaksList' not in self.Background[1]:
            self.Background[1]['peaksList'] = []
        flags = 4*[False]
        for i,f in enumerate(refflags):
            if i>3: break
            flags[i] = bool(f)
        bpk = []
        for i,j in zip((pos,int,sig,gam),flags):
            bpk += [float(i),j]
        self.Background[1]['peaksList'].append(bpk)
        self.Background[1]['nPeaks'] = len(self.Background[1]['peaksList'])

    def del_back_peak(self,peaknum):
        '''Removes a background peak from the Background parameters

        :param int peaknum: the number of the peak (starting from 0)
        '''
        npks = self.Background[1].get('nPeaks',0)
        if peaknum >= npks:
            raise Exception('peak {} not found in histogram {}'.format(peaknum,self.name))
        del self.Background[1]['peaksList'][peaknum]
        self.Background[1]['nPeaks'] = len(self.Background[1]['peaksList'])

    def ref_back_peak(self,peaknum,refflags=[]):
        '''Sets refinement flag for a background peak

        :param int peaknum: the number of the peak (starting from 0)
        :param list refflags: a list of 1 to 4 boolean refinement flags for
            pos,int,sig & gam, respectively. If a flag is not specified
            it defaults to False (use [0,1] to refine int only).
            Defaults to [] which means nothing is refined.
        '''
        npks = self.Background[1].get('nPeaks',0)
        if peaknum >= npks:
            raise Exception('peak {} not found in histogram {}'.format(peaknum,self.name))
        flags = 4*[False]
        for i,f in enumerate(refflags):
            if i>3: break
            flags[i] = bool(f)
        for i,f in enumerate(flags):
            self.Background[1]['peaksList'][peaknum][2*i+1] = f

    @property
    def id(self):
        self.proj.update_ids()
        return self.data['data'][0]['hId']

    @id.setter
    def id(self, val):
        self.data['data'][0]['hId'] = val

    def Limits(self,typ,value=None):
        '''Used to obtain or set the histogram limits. When a value is
        specified, the appropriate limit is set. Otherwise, the value is
        returned. Note that this provides an alternative to setting
        histogram limits with the :meth:`G2Project:do_refinements` or
        :meth:`G2PwdrData.set_refinements` methods.

        :param str typ: a string which must be either 'lower'
          (for 2-theta min or TOF min) or 'upper' (for 2theta max or TOF max).
          Anything else produces an error.
        :param float value: the number to set the limit (in units of degrees
          or TOF (microsec.). If not specified, the command returns the
          selected limit value rather than setting it.
        :returns: The current value of the requested limit
          (when ``value=None``). Units are
          2-theta (degrees) or TOF (microsec).

        Examples::

          h = gpx.histogram(0)
          val = h.Limits('lower')
          h.Limits('upper',75)
        '''

        if value == None:
            if typ == 'lower':
                return self.data['Limits'][1][0]
            elif typ == 'upper':
                return self.data['Limits'][1][1]
            else:
                raise G2ScriptException('Limit errr: must be "lower" or "upper"')
        else:
            if typ == 'lower':
                self.data['Limits'][1][0] = float(value)
            elif typ == 'upper':
                self.data['Limits'][1][1] = float(value)
            else:
                raise G2ScriptException('G2PwdData.Limit error: must be "lower" or "upper"')

    def Excluded(self,value=None):
        '''Used to obtain or set the excluded regions for a histogram.
        When a value is specified, the excluded regions are set.
        Otherwise, the list of excluded region pairs is returned.
        Note that excluded regions may be an empty list or a list
        of regions to be excluded, where each region is provided
        as pair of numbers, where the lower limit comes first.
        Some sample excluded region lists are::

          [[4.5, 5.5], [8.0, 9.0]]

          [[130000.0, 140000.0], [160000.0, 170000.0]]

          []

        The first above describes two excluded regions from 4.5-5.5 and 8-9
        degrees 2-theta. The second is for a TOF pattern and also
        describes two excluded regions, for 130-140 and 160-170 milliseconds.
        The third line would be the case where there are no excluded
        regions.

        :param list value: A list of pairs of excluded region numbers
          (as two-element lists). Some error checking/reformatting is done,
          but users are expected to get this right. Use the GUI to
          create examples or check input. Numbers in the list
          are in units of degrees or TOF (microsec.).

          If a value is not specified, the command returns the list of excluded regions.

        :returns: The list of excluded regions (when ``value=None``). Units are
          2-theta (degrees) or TOF (microsec).

        Example 1::

          h = gpx.histogram(0)  # adds an excluded region (11-13 degrees)
          h.Excluded(h.Excluded() + [[11,13]])

        Example 2::

          h = gpx.histogram(0) # changes the range of the first excluded region
          excl = h.Excluded()
          excl[0] = [120000.0, 160000.0]  # microsec
          h.Excluded(excl)

        Example 3::

          h = gpx.histogram(0)  # deletes all excluded regions
          h.Excluded([])
        '''
        if value is None:
            return self.data['Limits'][2:]

        # set up the excluded regions. Make sure we have a list of list pairs
        s = ''
        listValues = []
        try:
            for i,(l,h) in enumerate(value): # quick bit of validation
                listValues.append([
                    min(float(l),float(h)),
                    max(float(l),float(h)),
                    ])
                if float(l) < self.data['Limits'][0][0]:
                    s += f'Lower limit {l} too low. '
                if float(h) > self.data['Limits'][0][1]:
                    s += f'Upper limit {h} too high. '
        except:
            raise G2ScriptException(f'G2PwdData.Excluded error: incorrectly formatted list or\n\tinvalid value. In: {value}')
        if s:
            raise G2ScriptException(f'G2PwdData.Excluded error(s): {s}')
        self.data['Limits'][2:] = listValues

    def fit_fixed_points(self):
        """Attempts to apply a background fit to the fixed points currently specified."""
        def SetInstParms(Inst):
            dataType = Inst['Type'][0]
            insVary = []
            insNames = []
            insVals = []
            for parm in Inst:
                insNames.append(parm)
                insVals.append(Inst[parm][1])
                if parm in ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha',
                    'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',] and Inst[parm][2]:
                        Inst[parm][2] = False
            instDict = dict(zip(insNames, insVals))
            instDict['X'] = max(instDict['X'], 0.01)
            instDict['Y'] = max(instDict['Y'], 0.01)
            if 'SH/L' in instDict:
                instDict['SH/L'] = max(instDict['SH/L'], 0.002)
            return dataType, instDict, insVary

        import scipy.interpolate as si
        bgrnd = self.data['Background']

        # Need our fixed points in order
        bgrnd[1]['FixedPoints'].sort(key=lambda pair: pair[0])
        X = [x for x, y in bgrnd[1]['FixedPoints']]
        Y = [y for x, y in bgrnd[1]['FixedPoints']]

        limits = self.data['Limits'][1]
        if X[0] > limits[0]:
            X = [limits[0]] + X
            Y = [Y[0]] + Y
        if X[-1] < limits[1]:
            X += [limits[1]]
            Y += [Y[-1]]

        # Some simple lookups
        controls = self.proj['Controls']['data']
        inst, inst2 = self.data['Instrument Parameters']
        pwddata = self.data['data'][1]

        # Construct the data for background fitting
        xBeg = np.searchsorted(pwddata[0], limits[0])
        xFin = np.searchsorted(pwddata[0], limits[1])
        xdata = pwddata[0][xBeg:xFin]
        ydata = si.interp1d(X,Y)(ma.getdata(xdata))

        W = [1]*len(xdata)
        Z = [0]*len(xdata)

        dataType, insDict, insVary = SetInstParms(inst)
        bakType, bakDict, bakVary = G2pwd.SetBackgroundParms(bgrnd)

        # Do the fit
        data = np.array([xdata, ydata, W, Z, Z, Z])
        G2pwd.DoPeakFit('LSQ', [], bgrnd, limits, inst, inst2, data,
                        prevVaryList=bakVary, controls=controls)

        # Post-fit
        parmDict = {}
        bakType, bakDict, bakVary = G2pwd.SetBackgroundParms(bgrnd)
        parmDict.update(bakDict)
        parmDict.update(insDict)
        pwddata[3][xBeg:xFin] *= 0
        pwddata[5][xBeg:xFin] *= 0
        pwddata[4][xBeg:xFin] = G2pwd.getBackground('', parmDict, bakType, dataType, xdata)[0]

        # TODO adjust pwddata? GSASIIpwdGUI.py:1041
        # TODO update background
        self.proj.save()

    def getdata(self,datatype):
        '''Provides access to the histogram data of the selected data type

        :param str datatype: must be one of the following values
          (case is ignored):

           * 'X': the 2theta or TOF values for the pattern
           * 'Q': the 2theta or TOF values for the pattern transformed to Q
           * 'd': the 2theta or TOF values for the pattern transformed to d-space
           * 'Yobs': the observed intensity values
           * 'Yweight': the weights for each data point (1/sigma**2)
           * 'Ycalc': the computed intensity values
           * 'Background': the computed background values
           * 'Residual': the difference between Yobs and Ycalc (obs-calc)

        :returns: an numpy MaskedArray with data values of the requested type

        '''
        enums = ['x', 'yobs', 'yweight', 'ycalc', 'background', 'residual', 'q', 'd']
        if datatype.lower() not in enums:
            raise G2ScriptException("Invalid datatype = "+datatype+" must be one of "+str(enums))
        if datatype.lower() == 'q':
            Inst,Inst2 = self.data['Instrument Parameters']
            return  2 * np.pi / G2lat.Pos2dsp(Inst,self.data['data'][1][0])
        elif datatype.lower() == 'd':
            Inst,Inst2 = self.data['Instrument Parameters']
            return G2lat.Pos2dsp(Inst,self.data['data'][1][0])
        else:
            return self.data['data'][1][enums.index(datatype.lower())]

    def y_calc(self):
        '''Returns the calculated intensity values; better to
        use :meth:`getdata`
        '''
        return self.data['data'][1][3]

    def reflections(self):
        '''Returns a dict with an entry for every phase in the
        current histogram. Within each entry is a dict with keys
        'RefList' (reflection list, see
        :ref:`Powder Reflections <PowderRefl_table>`),
        'Type' (histogram type), 'FF'
        (form factor information), 'Super' (True if this is superspace
        group).
        '''
        return self.data['Reflection Lists']

    def Export(self,fileroot,extension,fmthint=None):
        '''Write the histogram into a file. The path is specified by fileroot and
        extension.

        :param str fileroot: name of the file, optionally with a path (extension is
           ignored)
        :param str extension: includes '.', must match an extension in global
           exportersByExtension['powder'] or a Exception is raised.
        :param str fmthint: If specified, the first exporter where the format
           name (obj.formatName, as shown in Export menu) contains the
           supplied string will be used. If not specified, an error
           will be generated showing the possible choices.
        :returns: name of file that was written
        '''
        LoadG2fil()
        if extension not in exportersByExtension.get('powder',[]):
            print('Defined exporters are:')
            print('  ',list(exportersByExtension.get('powder',[])))
            raise G2ScriptException('No Writer for file type = "'+extension+'"')
        fil = os.path.abspath(os.path.splitext(fileroot)[0]+extension)
        obj = exportersByExtension['powder'][extension]
        if type(obj) is list:
            if fmthint is None:
                print('Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('No format hint for file type = "'+extension+'"')
            for o in obj:
              if fmthint.lower() in o.formatName.lower():
                  obj = o
                  break
            else:
                print('Hint ',fmthint,'not found. Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('Bad format hint for file type = "'+extension+'"')
        self._SetFromArray(obj)
        obj.Writer(self.name,fil)
        return fil

    def _SetFromArray(self,expObj):
        '''Load a histogram into the exporter in preparation for use of
        the .Writer method in the object.

        :param Exporter expObj: Exporter object
        '''
        expObj.Histograms[self.name] =  {}
        expObj.Histograms[self.name]['Data'] = self.data['data'][1]
        for key in 'Instrument Parameters','Sample Parameters','Reflection Lists':
            expObj.Histograms[self.name][key] = self.data[key]

    def plot(self, Yobs=True, Ycalc=True, Background=True, Residual=True):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            G2fil.G2Print('Warning: Unable to import matplotlib, skipping plot')
            return
        data = self.data['data'][1]
        if Yobs:
            plt.plot(data[0], data[1], label='Yobs')
        if Ycalc:
            plt.plot(data[0], data[3], label='Ycalc')
        if Background:
            plt.plot(data[0], data[4], label='Background')
        if Residual:
            plt.plot(data[0], data[5], label="Residual")

    def get_wR(self):
        """returns the overall weighted profile R factor for a histogram

        :returns: a wR value as a percentage or None if not defined
        """
        return self['data'][0].get('wR')

    def _decodeHist(self,hist):
        '''Convert a histogram reference to a histogram name string
        '''
        if isinstance(hist, G2PwdrData):
            return hist.name
        elif isinstance(hist, G2Single):
            return hist.name
        elif hist in [h.name for h in self.proj.histograms()]:
            return hist
        elif type(hist) is int:
            return self.proj.histograms()[hist].name
        else:
            raise G2ScriptException("Invalid histogram reference: "+str(hist))

    def set_background(self, key, value):
        '''Set background parameters (this serves a similar function as in
        :meth:`set_refinements`, but with a simplified interface).

        :param str key: a string that defines the background parameter that will
           be changed. Must appear in the table below.

           =================   ==============   ===========================================
           key name            type of value     meaning of value
           =================   ==============   ===========================================
           fixedHist           int, str,         reference to a histogram in the current
                               None or           project or None to remove the reference.
                               G2PwdrData
           fixedFileMult       float             multiplier applied to intensities in
                                                 the background histogram where a value
                                                 of -1.0 means full subtraction of
                                                 the background histogram.
           =================   ==============   ===========================================

        :param value: a value to set the selected background parameter. The meaning
           and type for this parameter is listed in the table above.

        '''
        bkgPrms, bkgDict = self.data['Background']
        if key == 'fixedHist':
            if value is None:
                bkgDict['background PWDR'][0] = ''
                return
            bkgDict['background PWDR'][0] = self._decodeHist(value)
        elif key == 'fixedFileMult':
            bkgDict['background PWDR'][1] = float(value)
        else:
            raise ValueError("Invalid key in set_background:", key)

    def set_refinements(self, refs):
        """Sets the PWDR histogram refinement parameter 'key' to the specification 'value'.

        :param dict refs: A dictionary of the parameters to be set. See the
                          :ref:`Histogram_parameters_table` table for a description of
                          what these dictionaries should be.

        :returns: None

        """
        do_fit_fixed_points = False
        for key, value in refs.items():
            if key == 'Limits':
                old_limits = self.data['Limits'][1]
                new_limits = value
                if isinstance(new_limits, dict):
                    if 'low' in new_limits:
                        old_limits[0] = new_limits['low']
                    if 'high' in new_limits:
                        old_limits[1] = new_limits['high']
                else:
                    old_limits[0], old_limits[1] = new_limits
            elif key == 'Sample Parameters':
                sample = self.data['Sample Parameters']
                for sparam in value:
                    if sparam not in sample:
                        raise ValueError("Unknown refinement parameter, "
                                         + str(sparam))
                    sample[sparam][1] = True
            elif key == 'Background':
                bkg, peaks = self.data['Background']

                # If True or False, just set the refine parameter
                if value in (True, False):
                    bkg[1] = value
                    return

                if 'type' in value:
                    bkg[0] = value['type']
                if 'refine' in value:
                    bkg[1] = value['refine']
                if 'no. coeffs' in value:
                    cur_coeffs = bkg[2]
                    n_coeffs = value['no. coeffs']
                    if n_coeffs > cur_coeffs:
                        for x in range(n_coeffs - cur_coeffs):
                            bkg.append(0.0)
                    else:
                        for _ in range(cur_coeffs - n_coeffs):
                            bkg.pop()
                    bkg[2] = n_coeffs
                if 'coeffs' in value:
                    bkg[3:] = value['coeffs']
                if 'FixedPoints' in value:
                    peaks['FixedPoints'] = [(float(a), float(b))
                                            for a, b in value['FixedPoints']]
                if value.get('fit fixed points', False):
                    do_fit_fixed_points = True
                if 'peaks' in value:
                    for i,flags in enumerate(value['peaks']):
                        self.ref_back_peak(i,flags)

            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    errmsg = ''
                    try:
                        instrument[iparam][2] = True
                    except IndexError:
                        errmsg = f"Invalid key: {iparam}"
                    if errmsg: raise ValueError(errmsg)
            else:
                raise ValueError("Unknown key:", key)
        # Fit fixed points after the fact - ensure they are after fixed points
        # are added
        if do_fit_fixed_points:
            # Background won't be fit if refinement flag not set
            orig = self.data['Background'][0][1]
            self.data['Background'][0][1] = True
            self.fit_fixed_points()
            # Restore the previous value
            self.data['Background'][0][1] = orig

    def clear_refinements(self, refs):
        """Clears the PWDR refinement parameter 'key' and its associated value.

        :param dict refs: A dictionary of parameters to clear.
          See the :ref:`Histogram_parameters_table` table for what can be specified.
        """
        for key, value in refs.items():
            if key == 'Limits':
                old_limits, cur_limits = self.data['Limits']
                cur_limits[0], cur_limits[1] = old_limits
            elif key == 'Sample Parameters':
                sample = self.data['Sample Parameters']
                for sparam in value:
                    sample[sparam][1] = False
            elif key == 'Background':
                bkg, peaks = self.data['Background']

                # If True or False, just set the refine parameter
                if value in (True, False):
                    bkg[1] = False
                    return

                bkg[1] = False
                if 'FixedPoints' in value:
                    if 'FixedPoints' in peaks:
                        del peaks['FixedPoints']
                if 'peaks' in value:
                    for i in range(len(self.Background[1]['peaksList'])):
                        self.ref_back_peak(i,[])
            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    instrument[iparam][2] = False
            else:
                raise ValueError("Unknown key:", key)

    def add_peak(self,area,dspace=None,Q=None,ttheta=None):
        '''Adds a single peak to the peak list
        :param float area: peak area
        :param float dspace: peak position as d-space (A)
        :param float Q: peak position as Q (A-1)
        :param float ttheta: peak position as 2Theta (deg)

        Note: only one of the parameters: dspace, Q or ttheta may be specified.
        See :ref:`PeakRefine` for an example.
        '''
        if (not dspace) + (not Q) + (not ttheta) != 2:
            G2fil.G2Print('add_peak error: too many or no peak position(s) specified')
            return
        pos = ttheta
        Parms,Parms2 = self.data['Instrument Parameters']
        if Q:
            pos = G2lat.Dsp2pos(Parms,2.*np.pi/Q)
        elif dspace:
            pos = G2lat.Dsp2pos(Parms,dspace)
        peaks = self.data['Peak List']
        peaks['sigDict'] = {}        #no longer valid
        peaks['peaks'].append(G2mth.setPeakparms(Parms,Parms2,pos,area))

    def set_peakFlags(self,peaklist=None,area=None,pos=None,sig=None,gam=None,
                          alp=None,bet=None):
        '''Set refinement flags for peaks

        :param list peaklist: a list of peaks to change flags. If None (default), changes
          are made to all peaks.
        :param bool area: Sets or clears the refinement flag for the peak area value.
          If None (the default), no change is made.
        :param bool pos: Sets or clears the refinement flag for the peak position value.
          If None (the default), no change is made.
        :param bool sig: Sets or clears the refinement flag for the peak sigma (Gaussian width) value.
          If None (the default), no change is made.
        :param bool gam: Sets or clears the refinement flag for the peak gamma (Lorentzian width) value.
          If None (the default), no change is made.
        :param bool alp: Sets or clears the refinement flag for the peak alpha (TOF width) value.
          If None (the default), no change is made.
        :param bool bet: Sets or clears the refinement flag for the peak beta (TOF width) value.
          If None (the default), no change is made.

        Note that when peaks are first created the area flag is on and the other flags are
        initially off.

        Example::

           set_peakFlags(sig=False,gam=True)

        causes the sig refinement flag to be cleared and the gam flag to be set, in both cases for
        all peaks. The position and area flags are not changed from their previous values.
        '''
        peaks = self.data['Peak List']
        if peaklist is None:
            peaklist = range(len(peaks['peaks']))
        for i in peaklist:
            if 'C' in self.data['Instrument Parameters'][0]['Type'][0]:
                for var,j in [(area,3),(pos,1),(sig,5),(gam,7)]:
                    if var is not None:
                        peaks['peaks'][i][j] = var
            else:
                for var,j in [(area,3),(pos,1),(alp,5),(bet,7),(sig,9),(gam,11)]:
                    if var is not None:
                        peaks['peaks'][i][j] = var

    def refine_peaks(self,mode='useIP'):
        '''Causes a refinement of peak position, background and instrument parameters

        :param str mode: this determines how peak widths are determined. If
          the value is 'useIP' (the default) then the width parameter values (sigma, gamma,
          alpha,...) are computed from the histogram's instrument parameters. If the
          value is 'hold', then peak width parameters are not overridden. In
          this case, it is not possible to refine the instrument parameters
          associated with the peak widths and an attempt to do so will result in
          an error.

        :returns: a list of dicts with refinement results. Element 0 has uncertainties
          on refined values (also placed in self.data['Peak List']['sigDict'])
          element 1 has the peak fit result, element 2 has the peak fit uncertainties
          and element 3 has r-factors from the fit.
          (These are generated in :meth:`GSASIIpwd.DoPeakFit`).
        '''
        if mode == 'useIP':
            G2pwd.setPeakInstPrmMode(True)
        elif mode == 'hold':
            G2pwd.setPeakInstPrmMode(False)
        else:
            raise G2ScriptException('Error: invalid mode value in refine_peaks (allowed: "useIP" or "hold")')
        controls = self.proj.data.get('Controls',{})
        controls = controls.get('data',
                            {'deriv type':'analytic','min dM/M':0.001,}     #fill in defaults if needed
                            )
        peaks = self.data['Peak List']
        Parms,Parms2 = self.data['Instrument Parameters']
        background = self.data['Background']
        limits = self.data['Limits'][1]
        bxye = np.zeros(len(self.data['data'][1][1]))
        result = G2pwd.DoPeakFit('LSQ',peaks['peaks'],background,limits,
                                           Parms,Parms2,self.data['data'][1],bxye,[],
                                           oneCycle=False,controls=controls,dlg=None)
        peaks['sigDict'] = result[0]
        return result

    @property
    def Peaks(self):
        '''Provides a dict with the Peak List parameters
        for this histogram.

        :returns: dict with two elements where item
          'peaks' is a list of peaks where each element is
          [pos,pos-ref,area,area-ref,sig,sig-ref,gam,gam-ref],
          where the -ref items are refinement flags and item
          'sigDict' is a dict with possible items 'Back;#',
          'pos#', 'int#', 'sig#', 'gam#'
        '''
        return self.data['Peak List']

    @property
    def PeakList(self):
        '''Provides a list of peaks parameters
        for this histogram.

        :returns: a list of peaks, where each peak is a list containing
          [pos,area,sig,gam]
          (position, peak area, Gaussian width, Lorentzian width)

        '''
        return [i[::2] for i in self.data['Peak List']['peaks']]

    def Export_peaks(self,filename):
        '''Write the peaks file. The path is specified by filename
        extension.

        :param str filename: name of the file, optionally with a path,
            includes an extension
        :returns: name of file that was written
        '''
        import math
        nptand = lambda x: np.tan(x*math.pi/180.)
        fil = os.path.abspath(filename)
        fp = open(filename,'w')
        Inst,Inst2 = self.data['Instrument Parameters']
        Type = Inst['Type'][0]
        if 'T' not in Type:
            wave = G2mth.getWave(Inst)
        else:
            wave = None
        pkdata = self.data['Peak List']
        peaks = pkdata['peaks']
        sigDict = pkdata['sigDict']
        # code taken from GSASIIdataGUI OnExportPeakList
        fp.write("#%s \n" % (self.name+' Peak List'))
        if wave:
            fp.write('#wavelength = %10.6f\n'%(wave))
        if 'T' in Type:
            fp.write('#%9s %10s %10s %12s %10s %10s %10s %10s %10s\n'%('pos','dsp','esd','int','alp','bet','sig','gam','FWHM'))
        else:
            fp.write('#%9s %10s %10s %12s %10s %10s %10s\n'%('pos','dsp','esd','int','sig','gam','FWHM'))
        for ip,peak in enumerate(peaks):
            dsp = G2lat.Pos2dsp(Inst,peak[0])
            if 'T' in Type:  #TOF - more cols
                esds = {'pos':0.,'int':0.,'alp':0.,'bet':0.,'sig':0.,'gam':0.}
                for name in list(esds.keys()):
                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                sig = np.sqrt(peak[8])
                gam = peak[10]
                esddsp = G2lat.Pos2dsp(Inst,esds['pos'])
                FWHM = G2pwd.getgamFW(gam,sig) +(peak[4]+peak[6])*np.log(2.)/(peak[4]*peak[6])     #to get delta-TOF from Gam(peak)
                fp.write("%10.2f %10.5f %10.5f %12.2f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % \
                    (peak[0],dsp,esddsp,peak[2],peak[4],peak[6],peak[8],peak[10],FWHM))
            else:               #CW
                #get esds from sigDict for each peak & put in output - esds for sig & gam from UVWXY?
                esds = {'pos':0.,'int':0.,'sig':0.,'gam':0.}
                for name in list(esds.keys()):
                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                sig = np.sqrt(peak[4]) #var -> sig
                gam = peak[6]
                esddsp = 0.5*esds['pos']*dsp/nptand(peak[0]/2.)
                FWHM = G2pwd.getgamFW(gam,sig)      #to get delta-2-theta in deg. from Gam(peak)
                fp.write("%10.4f %10.5f %10.5f %12.2f %10.5f %10.5f %10.5f \n" % \
                    (peak[0],dsp,esddsp,peak[2],np.sqrt(max(0.0001,peak[4]))/100.,peak[6]/100.,FWHM/100.)) #convert to deg
        fp.close()
        return fil

    def SaveProfile(self,filename):
        '''Writes a GSAS-II (new style) .instprm file
        '''
        data,Parms2 = self.data['Instrument Parameters']
        filename = os.path.splitext(filename)[0]+'.instprm'         # make sure extension is .instprm
        File = open(filename,'w')
        File.write("#GSAS-II instrument parameter file; do not add/delete items!\n")
        for item in data:
            File.write(item+':'+str(data[item][1])+'\n')
        File.close()
        G2fil.G2Print ('Instrument parameters saved to: '+filename)

    def LoadProfile(self,filename,bank=0):
        '''Reads a GSAS-II (new style) .instprm file and overwrites the current
        parameters

        :param str filename: instrument parameter file name, extension ignored if not
          .instprm
        :param int bank: bank number to read, defaults to zero
        '''
        filename = os.path.splitext(filename)[0]+'.instprm'         # make sure extension is .instprm
        File = open(filename,'r')
        S = File.readline()
        newItems = []
        newVals = []
        Found = False
        while S:
            if S[0] == '#':
                if Found:
                    break
                if 'Bank' in S:
                    if bank == int(S.split(':')[0].split()[1]):
                        S = File.readline()
                        continue
                    else:
                        S = File.readline()
                        while S and '#Bank' not in S:
                            S = File.readline()
                        continue
                else:   #a non #Bank file
                    S = File.readline()
                    continue
            Found = True
            [item,val] = S[:-1].split(':')
            newItems.append(item)
            try:
                newVals.append(float(val))
            except ValueError:
                newVals.append(val)
            S = File.readline()
        File.close()
        LoadG2fil()
        self.data['Instrument Parameters'][0] = G2fil.makeInstDict(newItems,newVals,len(newVals)*[False,])

    def EditSimulated(self,Tmin, Tmax, Tstep=None, Npoints=None):
        '''Change the parameters for an existing simulated powder histogram.
        This will reset the previously computed "observed" pattern.

        :param float Tmin: Minimum 2theta or TOF (microsec) for dataset to be simulated
        :param float Tmax: Maximum 2theta or TOF (usec) for dataset to be simulated
        :param float Tstep: Step size in 2theta or TOF (usec) for dataset to be simulated
           Default is to compute this from Npoints.
        :param int Νpoints: the number of data points to be used for computing the
            diffraction pattern. Defaults as None, which sets this to 2500. Do not specify
            both Npoints and Tstep. Due to roundoff the actual nuber of points used may differ
            by +-1 from Npoints. Must be below 25,000.
         '''
        if not self.data['data'][0]['Dummy']:
            raise G2ScriptException("Error: histogram for G2PwdrData.EditSimulated is not simulated")
        if Tmax < Tmin:
            Tmin,Tmax = Tmax,Tmin
        if Tstep is not None and Npoints is not None:
            raise G2ScriptException("Error: Tstep and Npoints both specified")
        elif Tstep is not None:
            Tstep = abs(Tstep)
        elif Npoints is None:
            Npoints = 2500

        if 'T' in self.data['Instrument Parameters'][0]['Type'][0]:
            if Tmax > 200.:
                raise G2ScriptException("Error: Tmax is too large")
            if Npoints:
                N = Npoints
                Tstep = (np.log(Tmax)-np.log(Tmin))/N
            else:
                N = (np.log(Tmax)-np.log(Tmin))/Tstep
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.exp((np.arange(0,N))*Tstep+np.log(Tmin*1000.))
            N = len(x)
            unit = 'millisec'
            limits = [(1000*Tmin, 1000*Tmax), [1000*Tmin, 1000*Tmax]]
        else:
            if Npoints:
                N = Npoints
            else:
                N = int((Tmax-Tmin)/Tstep)+1
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.linspace(Tmin,Tmax,N,True)
            N = len(x)
            unit = 'degrees 2theta'
            limits = [(Tmin, Tmax), [Tmin, Tmax]]
        if N < 3:
            raise G2ScriptException("Error: Range is too small or step is too large, <3 points")
        G2fil.G2Print('Simulating {} points from {} to {} {}'.format(N,Tmin,Tmax,unit))
        self.data['data'][1] = [
            np.array(x), # x-axis values
            np.zeros_like(x), # powder pattern intensities
            np.ones_like(x), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(x), # calc. intensities (zero)
            np.zeros_like(x), # calc. background (zero)
            np.zeros_like(x), # obs-calc profiles
            ]
        self.data['Limits'] = limits

    def getHistEntryList(self, keyname=''):
        """Returns a dict with histogram setting values.

        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys
          is 'SGdata', etc.
        :returns: a set of histogram dict keys.

        See :meth:`G2Phase.getHAPentryList` for a related example.

        .. seealso::
            :meth:`getHistEntryValue`
            :meth:`setHistEntryValue`

        """
        return [i for i in dictDive(self.data,keyname) if i[0] != ['Histograms']]

    def getHistEntryValue(self, keylist):
        """Returns the histogram control value associated with a list of keys.
        Where the value returned is a list, it may be used as the target of
        an assignment (as in
        ``getHistEntryValue(...)[...] = val``)
        to set a value inside a list.

        :param list keylist: a list of dict keys, typically as returned by
          :meth:`getHistEntryList`.

        :returns: a histogram setting; may be a int, float, bool, list,...

        See :meth:`G2Phase.getHAPentryValue` for a related example.

        """
        d = self.data
        for key in keylist:
            d = d[key]
        return d

    def setHistEntryValue(self, keylist, newvalue):
        """Sets a histogram control value associated with a list of keys.

        See :meth:`G2Phase.setHAPentryValue` for a related example.

       :param list keylist: a list of dict keys, typically as returned by
          :meth:`getHistEntryList`.

        :param newvalue: a new value for the hist setting. The type must be
          the same as the initial value, but if the value is a container
          (list, tuple, np.array,...) the elements inside are not checked.

        """
        oldvalue = self.getHistEntryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("getHistEntryValue error: types do not agree for keys {}".format(keylist))
        d = self.data
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue

    def calc_autobkg(self,opt=0,logLam=None):
        """Sets fixed background points using the pybaselines Whittaker
        algorithm.

       :param int opt: 0 for 'arpls' or 1 for 'iarpls'. Default is 0.

       :param float logLam: log_10 of the Lambda value used in the
         pybaselines.whittaker.arpls/.iarpls computation. If None (default)
         is provided, a guess is taken for an appropriate value based
         on the number of points.

       :returns: the array of computed background points
       """
        bkgDict = self.data['Background'][1]
        xydata = self.data['data'][1]
        npts = len(xydata[1])
        bkgDict['autoPrms'] = bkgDict.get('autoPrms',{})
        try:
            opt = int(opt)
        except:
            opt = 0
        bkgDict['autoPrms']['opt'] = opt
        try:
            logLam = float(logLam)
        except:
            logLam =  min(10,float(int(10*np.log10(npts)**1.5)-9.5)/10.)
            print('Using default value of',logLam,'for pybaselines.whittaker.[i]arpls background computation')
        bkgDict['autoPrms']['logLam'] = logLam
        bkgDict['autoPrms']['Mode'] = None
        bkgdata = G2pwd.autoBkgCalc(bkgDict,xydata[1])
        bkgDict['FixedPoints'] = [i for i in zip(
                xydata[0].data[::npts//100],
                bkgdata.data[::npts//100])]
        bkgDict['autoPrms']['Mode'] = 'fixed'
        return bkgdata

class G2Phase(G2ObjectWrapper):
    """A wrapper object around a given phase.
    The object contains these class variables:

        * G2Phase.proj: contains a reference to the :class:`G2Project`
          object that contains this phase
        * G2Phase.name: contains the name of the phase
        * G2Phase.data: contains the phases's associated data in a dict,
          as documented for the :ref:`Phase Tree items<Phase_table>`.

    Scripts should not try to create a :class:`G2Phase` object directly as
    :meth:`G2Phase.__init__` should be invoked from inside :class:`G2Project`.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    def __init__(self, data, name, proj):
        self.data = data
        self.name = name
        self.proj = proj

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ["Cell", "Atoms", "LeBail"]
        return key in valid_keys

    @staticmethod
    def is_valid_HAP_refinement_key(key):
        valid_keys = ["Babinet", "Extinction", "HStrain", "Mustrain",
                      "Pref.Ori.", "Show", "Size", "Use", "Scale", "PhaseFraction"]
        return key in valid_keys

    def atom(self, atomlabel):
        """Returns the atom specified by atomlabel, or None if it does not
        exist.

        :param str atomlabel: The name of the atom (e.g. "O2")
        :returns: A :class:`G2AtomRecord` object
            representing the atom.
        """
        # Consult GSASIIobj.py for the meaning of this
        cx, ct, cs, cia = self.data['General']['AtomPtrs']
        ptrs = [cx, ct, cs, cia]
        atoms = self.data['Atoms']
        for atom in atoms:
            if atom[ct-1] == atomlabel:
                return G2AtomRecord(atom, ptrs, self.proj)

    def atoms(self):
        """Returns a list of atoms present in the current phase.

        :returns: A list of :class:`G2AtomRecord` objects.

        .. seealso::
            :meth:`~G2Phase.atom`
            :class:`G2AtomRecord`
        """
        ptrs = self.data['General']['AtomPtrs']
        return [G2AtomRecord(atom, ptrs, self.proj) for atom in self.data['Atoms']]

    def histograms(self):
        '''Returns a list of histogram names associated with the current phase ordered
        as they appear in the tree (see :meth:`G2Project.histograms`).
        '''
        return [i.name for i in self.proj.histograms() if i.name in self.data.get('Histograms', {})]

    @property
    def composition(self):
        '''Provides a dict where keys are atom types and values are the number of
        atoms of that type in cell (such as {'H': 2.0, 'O': 1.0})
        '''
        out = {}
        for a in self.atoms():
            typ = a.element
            if typ in out:
                out[typ] += a.mult*a.occupancy
            else:
                out[typ] = a.mult*a.occupancy
        return out

    def mu(self,wave):
        '''Provides mu values for a phase at the supplied wavelength in A.
        Uses GSASIImath.XScattDen which seems to be off by an order of
        magnitude, which has been corrected here.
        '''
        vol = self.data['General']['Cell'][7]
        out = {}
        for typ in self.data['General']['NoAtoms']:
            if typ in out:
                out[typ]['Num'] += self.data['General']['NoAtoms'][typ]
            else:
                out[typ] = {}
                out[typ]['Num'] = self.data['General']['NoAtoms'][typ]
                out[typ]['Z'] = 0 # wrong but not needed
        return 10*G2mth.XScattDen(out,vol,wave)[1]

    @property
    def density(self):
        '''Provides a scalar with the density of the phase. In case of a
        powder this assumes a 100% packing fraction.
        '''
        density,mattCoeff = G2mth.getDensity(self.data['General'])
        return density

    @property
    def ranId(self):
        return self.data['ranId']

    @property
    def id(self):
        return self.data['pId']

    @id.setter
    def id(self, val):
        self.data['pId'] = val

    def add_atom(self,x,y,z,element,lbl,occ=1.,uiso=0.01):
        '''Adds an atom to the current phase

        :param float x: atom fractional x coordinate
        :param float y: atom fractional y coordinate
        :param float z: atom fractional z coordinate
        :param str element: an element symbol (capitalization is ignored). Optionally add
          a valence (as in Ba+2)
        :param str lbl: A label for this atom
        :param float occ: A fractional occupancy for this atom (defaults to 1).
        :param float uiso: A Uiso value for this atom (defaults to 0.01).

        :returns: the :class:`~GSASIIscriptable.G2AtomRecord` atom object for the new atom
        '''
        x = float(x)
        y = float(y)
        z = float(z)
        occ = float(occ)
        uiso = float(uiso)

        generalData = self.data['General']
        atomData = self.data['Atoms']
        SGData = generalData['SGData']

        atId = ran.randint(0,sys.maxsize)
        Sytsym,Mult = G2spc.SytSym([x,y,z],SGData)[:2]

        if generalData['Type'] == 'macromolecular':
            atomData.append([0,lbl,'',lbl,element,'',x,y,z,occ,Sytsym,Mult,'I',uiso,0,0,0,0,0,0,atId])
        elif generalData['Type'] in ['nuclear','faulted',]:
            if generalData['Modulated']:
                atomData.append([lbl,element,'',x,y,z,occ,Sytsym,Mult,'I',uiso,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.append([lbl,element,'',x,y,z,occ,Sytsym,Mult,'I',uiso,0,0,0,0,0,0,atId])
        elif generalData['Type'] == 'magnetic':
            if generalData['Modulated']:
                atomData.append([lbl,element,'',x,y,z,occ,0.,0.,0.,Sytsym,Mult,'I',uiso,0,0,0,0,0,0,atId,[],[],
                    {'SS1':{'waveType':'Fourier','Sfrac':[],'Spos':[],'Sadp':[],'Smag':[]}}])
            else:
                atomData.append([lbl,element,'',x,y,z,occ,0.,0.,0.,Sytsym,Mult,'I',uiso,0,0,0,0,0,0,atId])

        SetupGeneral(self.data, None)
        #self.proj.index_ids()  # needed?
        #self.proj.update_ids()  # needed?
        return self.atom(lbl)

    def get_cell(self):
        """Returns a dictionary of the cell parameters, with keys:
            'length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma', 'volume'

        :returns: a dict

        .. seealso::
           :meth:`~G2Phase.get_cell_and_esd`

        """
        cell = self.data['General']['Cell']
        return {'length_a': cell[1], 'length_b': cell[2], 'length_c': cell[3],
                'angle_alpha': cell[4], 'angle_beta': cell[5], 'angle_gamma': cell[6],
                'volume': cell[7]}

    def get_cell_and_esd(self):
        """
        Returns a pair of dictionaries, the first representing the unit cell, the second
        representing the estimated standard deviations of the unit cell.

        :returns: a tuple of two dictionaries

        .. seealso::
           :meth:`~G2Phase.get_cell`

        """
        # translated from GSASIIfiles.ExportBaseclass.GetCell
        try:
            pfx = str(self.id) + '::'
            sgdata = self['General']['SGData']
            covDict = self.proj['Covariance']['data']

            parmDict = dict(zip(covDict.get('varyList',[]),
                                covDict.get('variables',[])))
            sigDict = dict(zip(covDict.get('varyList',[]),
                               covDict.get('sig',[])))

            if covDict.get('covMatrix') is not None:
                sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],covDict['varyList']))

            A, sigA = G2stIO.cellFill(pfx, sgdata, parmDict, sigDict)
            cellSig = G2lat.getCellEsd(pfx, sgdata, A, self.proj['Covariance']['data'])
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            cellDict, cellSigDict = {}, {}
            for i, key in enumerate(['length_a', 'length_b', 'length_c',
                                     'angle_alpha', 'angle_beta', 'angle_gamma',
                                     'volume']):
                cellDict[key] = cellList[i]
                cellSigDict[key] = cellSig[i]
            return cellDict, cellSigDict
        except KeyError:
            cell = self.get_cell()
            return cell, {key: 0.0 for key in cell}

    def export_CIF(self, outputname, quickmode=True):
        """Write this phase to a .cif file named outputname

        :param str outputname: The name of the .cif file to write to
        :param bool quickmode: Currently ignored. Carryover from exports.G2export_CIF"""
        # This code is all taken from exports/G2export_CIF.py
        # Functions copied have the same names

        from GSASII.exports import G2export_CIF as cif # delay as only needed here
#        CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        CIFname = os.path.splitext(self.proj.filename)[0]
        CIFname = os.path.split(CIFname)[1]
        CIFname = ''.join([c if ord(c) < 128 else ''
                           for c in CIFname.replace(' ', '_')])
        # try:
        #     author = self.proj['Controls']['data'].get('Author','').strip()
        # except KeyError:
        #     pass
        # oneblock = True

        covDict = self.proj['Covariance']['data']
        parmDict = dict(zip(covDict.get('varyList',[]),
                            covDict.get('variables',[])))
        sigDict = dict(zip(covDict.get('varyList',[]),
                           covDict.get('sig',[])))

        if covDict.get('covMatrix') is not None:
            sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],covDict['varyList'],noSym=True))

        with open(outputname, 'w') as fp:
            fp.write(' \n' + 70*'#' + '\n')
            cif.WriteCIFitem(fp, 'data_' + CIFname)
            # from exports.G2export_CIF.WritePhaseInfo
            cif.WriteCIFitem(fp, '\n# phase info for '+str(self.name) + ' follows')
            cif.WriteCIFitem(fp, '_pd_phase_name', self.name)
            # TODO get esds
            cellDict = self.get_cell()
#            defsigL = 3*[-0.00001] + 3*[-0.001] + [-0.01] # significance to use when no sigma
#            names = ['length_a','length_b','length_c',
#                     'angle_alpha','angle_beta ','angle_gamma',
#                     'volume']
            for key, val in cellDict.items():
                cif.WriteCIFitem(fp, '_cell_' + key, G2mth.ValEsd(val))

            cif.WriteCIFitem(fp, '_symmetry_cell_setting',
                         self.data['General']['SGData']['SGSys'])

            spacegroup = self.data['General']['SGData']['SpGrp'].strip()
            # regularize capitalization and remove trailing H/R
            spacegroup = spacegroup[0].upper() + spacegroup[1:].lower().rstrip('rh ')
            cif.WriteCIFitem(fp, '_symmetry_space_group_name_H-M', spacegroup)

            # generate symmetry operations including centering and center of symmetry
            SymOpList, offsetList, symOpList, G2oprList, G2opcodes = G2spc.AllOps(
                self.data['General']['SGData'])
            cif.WriteCIFitem(fp, 'loop_\n    _space_group_symop_id\n    _space_group_symop_operation_xyz')
            for i, op in enumerate(SymOpList,start=1):
                cif.WriteCIFitem(fp, '   {:3d}  {:}'.format(i,op.lower()))

            # TODO skipped histograms, exports/G2export_CIF.py:880

            # report atom params
            if self.data['General']['Type'] in ['nuclear','macromolecular']:        #this needs macromolecular variant, etc!
                cif.WriteAtomsNuclear(fp, self.data, self.name, parmDict, sigDict, [])
                # self._WriteAtomsNuclear(fp, parmDict, sigDict)
            else:
                raise G2ScriptException("no export for "+str(self.data['General']['Type'])+" coordinates implemented")
            # report cell contents
            cif.WriteComposition(fp, self.data, self.name, parmDict)
            if not quickmode and self.data['General']['Type'] == 'nuclear':      # report distances and angles
                # WriteDistances(fp,self.name,SymOpList,offsetList,symOpList,G2oprList)
                raise NotImplementedError("only quickmode currently supported")
            if 'Map' in self.data['General'] and 'minmax' in self.data['General']['Map']:
                cif.WriteCIFitem(fp,'\n# Difference density results')
                MinMax = self.data['General']['Map']['minmax']
                cif.WriteCIFitem(fp,'_refine_diff_density_max',G2mth.ValEsd(MinMax[0],-0.009))
                cif.WriteCIFitem(fp,'_refine_diff_density_min',G2mth.ValEsd(MinMax[1],-0.009))


    def set_refinements(self, refs):
        """Sets the phase refinement parameter 'key' to the specification 'value'

        :param dict refs: A dictionary of the parameters to be set. See the
                          :ref:`Phase_parameters_table` table for a description of
                          this dictionary.

        :returns: None"""
        for key, value in refs.items():
            if key == "Cell":
                self.data['General']['Cell'][0] = value

            elif key == "Atoms":
                for atomlabel, atomrefinement in value.items():
                    if atomlabel == 'all':
                        for atom in self.atoms():
                            atom.refinement_flags = atomrefinement
                    else:
                        atom = self.atom(atomlabel)
                        if atom is None:
                            raise ValueError("No such atom: " + atomlabel)
                        atom.refinement_flags = atomrefinement

            elif key == "LeBail":
                hists = self.data['Histograms']
                for hname, hoptions in hists.items():
                    hoptions['LeBail'] = bool(value)
            else:
                raise ValueError("Unknown key:", key)

    def clear_refinements(self, refs):
        """Clears a given set of parameters.

        :param dict refs: The parameters to clear.
          See the :ref:`Phase_parameters_table` table for what can be specified.
        """
        for key, value in refs.items():
            if key == "Cell":
                self.data['General']['Cell'][0] = False
            elif key == "Atoms":
                cx, ct, cs, cia = self.data['General']['AtomPtrs']

                for atomlabel in value:
                    atom = self.atom(atomlabel)
                    # Set refinement to none
                    atom.refinement_flags = ' '
            elif key == "LeBail":
                hists = self.data['Histograms']
                for hname, hoptions in hists.items():
                    hoptions['LeBail'] = False
            else:
                raise ValueError("Unknown key:", key)

    def set_HAP_refinements(self, refs, histograms='all'):
        """Sets the given HAP refinement parameters between the current phase and
        the specified histograms.

        :param dict refs: A dictionary of the parameters to be set. See
                          the :ref:`HAP_parameters_table` table for a description of this
                          dictionary.
        :param histograms: Either 'all' (default) or a list of the histograms by index, name
            or object. The index number is relative to all histograms in the tree, not to
            those in the phase.
            Histograms not associated with the current phase will be ignored.
            whose HAP parameters will be set with this phase. Histogram and phase
            must already be associated.
        :returns: None
        """
        if not self.data.get('Histograms',[]):
            G2fil.G2Print("Error likely: Phase {} has no linked histograms".format(self.name))
            return

        if histograms == 'all':
            histograms = self.data['Histograms'].keys()
        else:
            histograms = [self._decodeHist(h) for h in histograms
                          if self._decodeHist(h) in self.data['Histograms']]
        if not histograms:
            G2fil.G2Print("Warning: Skipping HAP set for phase {}, no selected histograms".format(self.name))
            return
        # remove non-PWDR (HKLF) histograms
        histograms = [i for i in histograms if self.proj.histType(i) == 'PWDR']
        if not histograms: return
        for key, val in refs.items():
            if key == 'Babinet':
                try:
                    sets = list(val)
                except ValueError:
                    sets = ['BabA', 'BabU']
                for param in sets:
                    if param not in ['BabA', 'BabU']:
                        raise ValueError("Not sure what to do with" + param)
                    for h in histograms:
                        self.data['Histograms'][h]['Babinet'][param][1] = True
            elif key == 'Extinction':
                for h in histograms:
                    self.data['Histograms'][h]['Extinction'][1] = bool(val)
            elif key == 'HStrain':
                if isinstance(val,list) or isinstance(val,tuple):
                    for h in histograms:
                        if len(self.data['Histograms'][h]['HStrain'][1]) != len(val):
                            raise Exception('Need {} HStrain terms for phase {} hist {}'
                                .format(len(self.data['Histograms'][h]['HStrain'][1]),self.name,h))
                        for i,v in enumerate(val):
                            self.data['Histograms'][h]['HStrain'][1][i] = bool(v)
                else:
                    for h in histograms:
                        self.data['Histograms'][h]['HStrain'][1] = [bool(val) for p in self.data['Histograms'][h]['HStrain'][1]]
            elif key == 'Mustrain':
                for h in histograms:
                    mustrain = self.data['Histograms'][h]['Mustrain']
                    newType = mustrain[0]
                    direction = None
                    if isinstance(val, (str,bytes)):
                        if val in ['isotropic', 'uniaxial', 'generalized']:
                            newType = val
                        else:
                            raise ValueError("Not a Mustrain type: " + val)
                    elif isinstance(val, dict):
                        newType = val.get('type', newType)
                        direction = val.get('direction', None)

                    if newType:
                        mustrain[0] = newType
                        if newType == 'isotropic':
                            mustrain[2][0] = True == val.get('refine',False)
                            mustrain[5] = [False for p in mustrain[4]]
                        elif newType == 'uniaxial':
                            if 'refine' in val:
                                mustrain[2][0] = False
                                types = val['refine']
                                if isinstance(types, (str,bytes)):
                                    types = [types]
                                elif isinstance(types, bool):
                                    mustrain[2][1] = types
                                    mustrain[2][2] = types
                                    types = []
                                else:
                                    raise ValueError("Not sure what to do with: "
                                                     + str(types))
                            else:
                                types = []

                            for unitype in types:
                                if unitype == 'equatorial':
                                    mustrain[2][0] = True
                                elif unitype == 'axial':
                                    mustrain[2][1] = True
                                else:
                                    msg = 'Invalid uniaxial mustrain type'
                                    raise ValueError(msg + ': ' + unitype)
                        else:  # newtype == 'generalized'
                            mustrain[2] = [False for p in mustrain[1]]
                            if 'refine' in val:
                                mustrain[5] = [True == val['refine']]*len(mustrain[5])

                    if direction:
                        if len(direction) != 3:
                            raise ValueError("Expected hkl, found", direction)
                        direction = [int(n) for n in direction]
                        mustrain[3] = direction
            elif key == 'Size':
                newSize = None
                if 'value' in val:
                    newSize = float(val['value'])
                for h in histograms:
                    size = self.data['Histograms'][h]['Size']
                    newType = size[0]
                    direction = None
                    if isinstance(val, (str,bytes)):
                        if val in ['isotropic', 'uniaxial', 'ellipsoidal']:
                            newType = val
                        else:
                            raise ValueError("Not a valid Size type: " + val)
                    elif isinstance(val, dict):
                        newType = val.get('type', size[0])
                        direction = val.get('direction', None)

                    if newType:
                        size[0] = newType
                        refine = bool(val.get('refine'))
                        if newType == 'isotropic' and refine is not None:
                            size[2][0] = bool(refine)
                            if newSize: size[1][0] = newSize
                        elif newType == 'uniaxial' and refine is not None:
                            size[2][1] = bool(refine)
                            size[2][2] = bool(refine)
                            if newSize: size[1][1] = size[1][2] =newSize
                        elif newType == 'ellipsoidal' and refine is not None:
                            size[5] = [bool(refine) for p in size[5]]
                            if newSize: size[4] = [newSize for p in size[4]]

                    if direction:
                        if len(direction) != 3:
                            raise ValueError("Expected hkl, found", direction)
                        direction = [int(n) for n in direction]
                        size[3] = direction
            elif key == 'Pref.Ori.':
                for h in histograms:
                    self.data['Histograms'][h]['Pref.Ori.'][2] = bool(val)
            elif key == 'Show':
                for h in histograms:
                    self.data['Histograms'][h]['Show'] = bool(val)
            elif key == 'Use':
                for h in histograms:
                    self.data['Histograms'][h]['Use'] = bool(val)
            elif key == 'Scale' or key == 'PhaseFraction':
                for h in histograms:
                    self.data['Histograms'][h]['Scale'][1] = bool(val)
            else:
                G2fil.G2Print('Warning: Unknown HAP key: '+key)

    def clear_HAP_refinements(self, refs, histograms='all'):
        """Clears the given HAP refinement parameters between this phase and
        the given histograms.

        :param dict refs: A dictionary of the parameters to be cleared.
            See the the :ref:`HAP_parameters_table` table for what can be specified.
        :param histograms: Either 'all' (default) or a list of the histograms by index, name
            or object.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
            Histograms not associated with the current phase will be ignored.
            whose HAP parameters will be set with this phase. Histogram and phase
            must already be associated
        :returns: None
        """
        if histograms == 'all':
            histograms = self.data['Histograms'].keys()
        else:
            histograms = [self._decodeHist(h) for h in histograms
                          if self._decodeHist(h) in self.data['Histograms']]

        for key, val in refs.items():
            for h in histograms:
                if key == 'Babinet':
                    try:
                        sets = list(val)
                    except ValueError:
                        sets = ['BabA', 'BabU']
                    for param in sets:
                        if param not in ['BabA', 'BabU']:
                            raise ValueError("Not sure what to do with" + param)
                        for h in histograms:
                            self.data['Histograms'][h]['Babinet'][param][1] = False
                elif key == 'Extinction':
                    for h in histograms:
                        self.data['Histograms'][h]['Extinction'][1] = False
                elif key == 'HStrain':
                    for h in histograms:
                        self.data['Histograms'][h]['HStrain'][1] = [False for p in self.data['Histograms'][h]['HStrain'][1]]
                elif key == 'Mustrain':
                    for h in histograms:
                        mustrain = self.data['Histograms'][h]['Mustrain']
                        mustrain[2] = [False for p in mustrain[2]]
                        mustrain[5] = [False for p in mustrain[4]]
                elif key == 'Pref.Ori.':
                    for h in histograms:
                        self.data['Histograms'][h]['Pref.Ori.'][2] = False
                elif key == 'Show':
                    for h in histograms:
                        self.data['Histograms'][h]['Show'] = False
                elif key == 'Size':
                    for h in histograms:
                        size = self.data['Histograms'][h]['Size']
                        size[2] = [False for p in size[2]]
                        size[5] = [False for p in size[5]]
                elif key == 'Use':
                    for h in histograms:
                        self.data['Histograms'][h]['Use'] = False
                elif key == 'Scale' or key == 'PhaseFraction':
                    for h in histograms:
                        self.data['Histograms'][h]['Scale'][1] = False
                else:
                    G2fil.G2Print('Warning: Unknown HAP key: '+key)

    def _decodeHist(self,hist):
        '''Convert a histogram reference to a histogram name string
        '''
        if isinstance(hist, G2PwdrData):
            return hist.name
        elif isinstance(hist, G2Single):
            return hist.name
        elif hist in self.data['Histograms']:
            return hist
        elif type(hist) is int:
            return self.proj.histograms()[hist].name
        else:
            raise G2ScriptException("Invalid histogram reference: "+str(hist))

    def getHAPvalues(self, histname):
        """Returns a dict with HAP values for the selected histogram

        :param histogram: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
        :returns: HAP value dict
        """
        return self.data['Histograms'][self._decodeHist(histname)]

    def copyHAPvalues(self, sourcehist, targethistlist='all', skip=[], use=None):
        """Copies HAP parameters for one histogram to a list of other histograms.
        Use skip or use to select specific entries to be copied or not used.

        :param sourcehist: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram to copy
            parameters from.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`),
            a histogram name or the index number of the histogram.
            If the string 'all' (default), then all histograms in the phase
            are used.
        :param list skip: items in the HAP dict that should not be
            copied. The default is an empty list, which causes all items
            to be copied. To see a list of items in the dict, use
            :meth:`getHAPvalues`.
            Don't use with :attr:`use`.
        :param list use: specifies the items in the HAP dict should be
            copied. The default is None, which causes all items
            to be copied.
            Don't use with :attr:`skip`.

        examples::

            ph0.copyHAPvalues(0,[1,2,3])
            ph0.copyHAPvalues(0,use=['HStrain','Size'])

        The first example copies all HAP parameters from the first histogram to
        the second, third and fourth histograms (as listed in the project tree).
        The second example copies only the 'HStrain' (Dij parameters and
        refinement flags) and the 'Size' (crystallite size settings, parameters
        and refinement flags) from the first histogram to all histograms.
        """
        sourcehist = self._decodeHist(sourcehist)
        if targethistlist == 'all':
            targethistlist = self.histograms()

        copydict = copy.deepcopy(self.data['Histograms'][sourcehist])
        for item in skip:
            if item == 'PhaseFraction': item = 'Scale'
            if item in list(copydict.keys()):
                del copydict[item]
            else:
                G2fil.G2Print('items in HAP dict are: {}'.format(
                    list(self.data['Histograms'][sourcehist])))
                raise Exception('HAP skip list entry {} invalid'.format(item))
        if use:
            for item in list(copydict.keys()):
                if item == 'PhaseFraction': item = 'Scale'
                if item not in use:
                    del copydict[item]
        txt = ', '.join(copydict.keys())
        G2fil.G2Print(f'copyHAPvalues for phase {self.name}')
        G2fil.G2Print(f'Copying item(s): {txt}\n from histogram: {sourcehist}')
        txt = '\n\t'.join([self._decodeHist(h) for h in targethistlist])
        G2fil.G2Print(f' to histogram(s):\n\t{txt}')
        for h in targethistlist:
            h = self._decodeHist(h)
            if h not in self.data['Histograms']:
                G2fil.G2Print('Unexpected Warning: histogram {} not in phase {}'.format(h,self.name))
                continue
            self.data['Histograms'][h].update(copy.deepcopy(copydict))

    def setSampleProfile(self, histname, parmType, mode, val1, val2=None, axis=None, LGmix=None):
        """Sets sample broadening parameters for a histogram associated with the
        current phase. This currently supports isotropic and uniaxial broadening
        modes only.

        :param histogram: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
        :param str parmType: should be 'size' or 'microstrain' (can be abbreviated to 's' or 'm')
        :param str mode: should be 'isotropic' or 'uniaxial' (can be abbreviated to 'i' or 'u')
        :param float val1: value for isotropic size (in :math:`\\mu m`) or
           microstrain (unitless, :math:`\\Delta Q/Q \\times 10^6`) or the equatorial value in the uniaxial case
        :param float val2: value for axial size (in :math:`\\mu m`) or
           axial microstrain (unitless, :math:`\\Delta Q/Q \\times 10^6`)
           in uniaxial case; not used for isotropic
        :param list axis: tuple or list with three values indicating the preferred direction
          for uniaxial broadening; not used for isotropic
        :param float LGmix: value for broadening type (1=Lorentzian, 0=Gaussian or a value
          between 0 and 1. Default value (None) is ignored.

        Examples::

            phase0.setSampleProfile(0,'size','iso',1.2)
            phase0.setSampleProfile(0,'micro','isotropic',1234)
            phase0.setSampleProfile(0,'m','u',1234,4567,[1,1,1],.5)
            phase0.setSampleProfile(0,'s','uni',1.2,2.3,[0,0,1])
        """
        if parmType.lower().startswith('s'):
            key = 'Size'
        elif parmType.lower().startswith('m'):
            key = 'Mustrain'
        else:
            G2fil.G2Print('setSampleProfile Error: value for parmType of {} is not size or microstrain'.
                              format(parmType))
            raise Exception('Invalid parameter in setSampleProfile')
        if mode.lower().startswith('i'):
            iso = True
        elif mode.lower().startswith('u'):
            iso = False
            if val2 is None:
                G2fil.G2Print('setSampleProfile Error: value for val2 is required with mode of uniaxial')
                raise Exception('Invalid val2 parameter in setSampleProfile')
            if axis is None:
                G2fil.G2Print('setSampleProfile Error: value for axis is required with mode of uniaxial')
                raise Exception('Invalid axis parameter in setSampleProfile')
        else:
            G2fil.G2Print('setSampleProfile Error: value for mode of {} is not isotropic or uniaxial'.
                              format(mode))
            raise Exception('Invalid parameter in setSampleProfile')

        d = self.data['Histograms'][self._decodeHist(histname)][key]
        if iso:
            d[0] = 'isotropic'
            d[1][0] = float(val1)
            if LGmix is not None: d[1][2] = float(LGmix)
        else:
            d[3] = [int(axis[0]),int(axis[1]),int(axis[2])]
            d[0] = 'uniaxial'
            d[1][0] = float(val1)
            d[1][1] = float(val2)
            if LGmix is not None: d[1][2] = float(LGmix)

    def HAPvalue(self, param=None, newValue=None, targethistlist='all'):
        """Retrieves or sets individual HAP parameters for one histogram or
        multiple histograms.

        :param str param: is a parameter name, which can be 'Scale' or
          'PhaseFraction' (either can be used for phase
          fraction), 'Use', 'Extinction', 'LeBail', 'PO' 
          (for Preferred Orientation).
          If not specified or invalid
          an exception is generated showing the list of valid parameters.
          At present, only these HAP parameters cannot be accessed with 
          this function: 'Size', 'Mustrain', 'HStrain', 'Babinet'. 
          This might be addressed in the future. 
          Some of these values can be set via
          :meth:`G2Phase.set_HAP_refinements`.
        :param newValue: the value to use when setting the HAP parameter for the
          appropriate histogram(s). Will be converted to the proper type or
          an exception will be generated if not possible. If not specified,
          and only one histogram is selected, the value is retrieved and
          returned.
          When param='PO' then this value is interpreted as the following:

            if the value is 0 or an even integer, then the preferred 
            orientation model is set to "Spherical harmonics". If the value is 
            1, then "March-Dollase" is used. Any other value generates an error.

        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`),
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
            If the string 'all' (default), then all histograms in the phase
            are used.

            targethistlist must correspond to a single histogram if a value
            is to be returned (i.e. when argument newValue is not specified).

        :returns: the value of the parameter, when argument newValue is not specified.

        .. seealso::
            :meth:`~G2Phase.set_HAP_refinements`

        Example::

            val = ph0.HAPvalue('Scale')
            val = ph0.HAPvalue('PhaseFraction',targethistlist=[0])
            ph0.HAPvalue('Scale',2.5)

        The first command returns the phase fraction if only one histogram
        is associated with the current phase, or raises an exception.
        The second command returns the phase fraction from the first histogram
        associated with the current phase. The third command sets the phase
        fraction for all histograms associated with the current phase.

        """
        doSet = not newValue is None

        if targethistlist == 'all':
            targethistlist = self.histograms()
        boolParam = ('Use','LeBail')
        refFloatParam = ('Scale','Extinction')
        useBool = False
        useFloat = False
        useInt = False
        if param == 'PhaseFraction': param = 'Scale'
        if param in boolParam:
            useBool = True
        elif param in refFloatParam:
            useFloat = True
        elif param in ['PO']:
            useInt = True
        else:
            s = ''
            for i in boolParam+refFloatParam+['PhaseFraction','PO']:
                if s != '': s += ', '
                s += f'"{i}"'
            raise G2ScriptException(f'Invalid parameter. Valid choices are: {s}')
        if not doSet and len(targethistlist) > 1:
            raise G2ScriptException(f'Unable to report value from {len(targethistlist)} histograms')
        for h in targethistlist:
            h = self._decodeHist(h)
            if h not in self.data['Histograms']:
                G2fil.G2Print('Warning: histogram {} not in phase {}'.format(h,self.name))
                continue
            if not doSet and useBool:
                return self.data['Histograms'][h][param]
            elif not doSet and useFloat:
                return self.data['Histograms'][h][param][0]
            elif useBool:
                self.data['Histograms'][h][param] = bool(newValue)
            elif useFloat:
                self.data['Histograms'][h][param][0] = float(newValue)
            elif useInt:
                if newValue is None:
                    return self.data['Histograms'][h]['Pref.Ori.']
                try:
                    intval = int(newValue)
                except:
                    intval = None
                if intval == 1:
                    self.data['Histograms'][h]['Pref.Ori.'][0] = 'MD'
                elif intval is not None and 2*int(intval//2) == intval:
                    SGData = self.data['General']['SGData']
                    cofNames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',intval,False)     #cylindrical & no M

                    self.data['Histograms'][h]['Pref.Ori.'][0] = 'SH'
                    self.data['Histograms'][h]['Pref.Ori.'][4] = intval
                    olddict = self.data['Histograms'][h]['Pref.Ori.'][5]
                    newdict = dict(zip(cofNames,len(cofNames)*[0.]))
                    # retain any old values in existing dict
                    newdict.update({i:olddict[i] for i in olddict if i in newdict})
                    self.data['Histograms'][h]['Pref.Ori.'][5] = newdict
                else:
                    raise G2ScriptException(f'Preferred orientation value of {newValue} is invalid')
            else:
                print('unexpected action')

    def setHAPvalues(self, HAPdict, targethistlist='all', skip=[], use=None):
        """Copies HAP parameters for one histogram to a list of other histograms.
        Use skip or use to select specific entries to be copied or not used.
        Note that ``HStrain`` and sometimes ``Mustrain`` values can be specific to
        a Laue class and should be copied with care between phases of different
        symmetry. A "sanity check" on the number of Dij terms is made if ``HStrain``
        values are copied.

        :param dict HAPdict: is a dict returned by :meth:`getHAPvalues` containing
            HAP parameters.
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`),
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to
            those in the phase.
            If the string 'all' (default), then all histograms in the phase
            are used.
        :param list skip: items in the HAP dict that should not be
            copied. The default is an empty list, which causes all items
            to be copied. To see a list of items in the dict, use
            :meth:`getHAPvalues`. Don't use with :attr:`use`.
        :param list use: specifies the items in the HAP dict should be
            copied. The default is None, which causes all items
            to be copied. Don't use with :attr:`skip`.

        Example::

            HAPdict = ph0.getHAPvalues(0)
            ph1.setHAPvalues(HAPdict,use=['HStrain','Size'])

        This copies the Dij (hydrostatic strain) HAP parameters and the
        crystallite size broadening terms from the first histogram in
        phase ``ph0`` to all histograms in phase ``ph1``.
        """
        if targethistlist == 'all':
            targethistlist = self.histograms()
        copydict = copy.deepcopy(HAPdict)
        for item in skip:
            if item == 'PhaseFraction': item = 'Scale'
            if item in list(copydict.keys()):
                del copydict[item]
            # else:
            #     G2fil.G2Print('items in HAP dict are: {}'.format(
            #         list(self.data['Histograms'][sourcehist])))
            #     raise Exception('HAP skip list entry {} invalid'.format(item))
        if use:
            for item in list(copydict.keys()):
                if item == 'PhaseFraction': item = 'Scale'
                if item not in use:
                    del copydict[item]

        first = True
        for h in targethistlist:
            h = self._decodeHist(h)
            if h not in self.data['Histograms']:
                G2fil.G2Print('Warning: histogram {} not in phase {}'.format(h,self.name))
                continue
            if first:
                first = False
                if 'HStrain' in self.data['Histograms'][h] and 'HStrain' in copydict:
                    if len(copydict['HStrain'][0]) != len(self.data['Histograms'][h]['HStrain'][0]):
                        G2fil.G2Print('Error: HStrain has differing numbers of terms. Input: {}, phase {}: {}'.
                                  format(len(copydict['HStrain'][0]),
                                        self.name,len(self.data['Histograms'][h]['HStrain'][0])))
                        raise Exception('HStrain has differing numbers of terms.')
            self.data['Histograms'][h].update(copy.deepcopy(copydict))
        G2fil.G2Print('Copied item(s) {} from dict'.format(list(copydict.keys())))
        G2fil.G2Print(' to histogram(s) {}'.format([self._decodeHist(h) for h in targethistlist]))

    def getPhaseEntryList(self, keyname=''):
        """Returns a dict with control values.

        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys
          is 'SGdata', etc.
        :returns: a set of phase dict keys. Note that HAP items, while
          technically part of the phase entries, are not included.

        See :meth:`getHAPentryList` for a related example.

        .. seealso::
            :meth:`getPhaseEntryValue`
            :meth:`setPhaseEntryValue`

        """
        return [i for i in dictDive(self.data,keyname) if i[0] != ['Histograms']]

    def getPhaseEntryValue(self, keylist):
        """Returns the value associated with a list of keys.
        Where the value returned is a list, it may be used as the target of
        an assignment (as in
        ``getPhaseEntryValue(...)[...] = val``)
        to set a value inside a list.

        :param list keylist: a list of dict keys, typically as returned by
          :meth:`getPhaseEntryList`.

        :returns: a phase setting; may be a int, float, bool, list,...

        See :meth:`getHAPentryValue` for a related example.

        """
        d = self.data
        for key in keylist:
            d = d[key]
        return d

    def setPhaseEntryValue(self, keylist, newvalue):
        """Sets a phase control value associated with a list of keys.

        :param list keylist: a list of dict keys, typically as returned by
          :meth:`getPhaseEntryList`.

        :param newvalue: a new value for the phase setting. The type must be
          the same as the initial value, but if the value is a container
          (list, tuple, np.array,...) the elements inside are not checked.

        See :meth:`setHAPentryValue` for a related example.

        """
        oldvalue = self.getPhaseEntryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("getPhaseEntryValue error: types do not agree for keys {}".format(keylist))
        d = self.data
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue

    def getHAPentryList(self, histname=None, keyname=''):
        """Returns a dict with HAP values. Optionally a histogram
        may be selected.

        :param histname: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to
            those in the phase. If no histogram is specified, all histograms
            are selected.
        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys
          is 'SGdata', etc.
        :returns: a set of HAP dict keys.

        Example:

        >>> p.getHAPentryList(0,'Scale')
        [(['PWDR test Bank 1', 'Scale'], list, [1.0, False])]

        .. seealso::
            :meth:`getHAPentryValue`
            :meth:`setHAPentryValue`

        """
        if histname:
            h = [self._decodeHist(histname)]
        else:
            h = []
        return dictDive(self.data['Histograms'],keyname,h)

    def getHAPentryValue(self, keylist):
        """Returns the HAP value associated with a list of keys. Where the
        value returned is a list, it may be used as the target of
        an assignment (as in
        ``getHAPentryValue(...)[...] = val``)
        to set a value inside a list.

        :param list keylist: a list of dict keys, typically as returned by
          :meth:`getHAPentryList`. Note the first entry is a histogram name.
          Example: ``['PWDR hist1.fxye Bank 1', 'Scale']``

        :returns: HAP value

        Example:

        >>> sclEnt = p.getHAPentryList(0,'Scale')[0]
        >>> sclEnt
        [(['PWDR test Bank 1', 'Scale'], list, [1.0, False])]
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, False]
        >>> p.getHAPentryValue(sclEnt[0])[1] = True
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, True]

        """
        d = self.data['Histograms']
        for key in keylist:
            d = d[key]
        return d

    def setHAPentryValue(self, keylist, newvalue):
        """Sets an HAP value associated with a list of keys.

        :param list keylist: a list of dict keys, typically as returned by
          :meth:`getHAPentryList`. Note the first entry is a histogram name.
          Example: ``['PWDR hist1.fxye Bank 1', 'Scale']``

        :param newvalue: a new value for the HAP setting. The type must be
          the same as the initial value, but if the value is a container
          (list, tuple, np.array,...) the elements inside are not checked.

        Example:

        >>> sclEnt = p.getHAPentryList(0,'Scale')[0]
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, False]
        >>> p.setHAPentryValue(sclEnt[0], (1, True))
        GSASIIscriptable.G2ScriptException: setHAPentryValue error: types do not agree for keys ['PWDR test.fxye Bank 1', 'Scale']
        >>> p.setHAPentryValue(sclEnt[0], [1, True])
        >>> p.getHAPentryValue(sclEnt[0])
        [1, True]

        """
        oldvalue = self.getHAPentryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("setHAPentryValue error: types do not agree for keys {}".format(keylist))
        d = self.data['Histograms']
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue

    def _getBondRest(self,nam):
        if 'Restraints' not in self.proj.data:
            raise G2ScriptException(f"{nam} error: Restraints entry not in data tree")
        errmsg = ''
        try:
            return self.proj.data['Restraints']['data'][self.name]['Bond']
        except:
            raise G2ScriptException(f"{nam} error: Bonds for phase not in Restraints")

    def setDistRestraintWeight(self, factor=1):
        '''Sets the weight for the bond distance restraint(s) to factor

        :param float factor: the weighting factor for this phase's restraints. Defaults
          to 1 but this value is typically much larger (10**2 to 10**4)

        .. seealso::
           :meth:`G2Phase.addDistRestraint`
        '''
        bondRestData = self._getBondRest('setDistRestraintWeight')
        bondRestData['wtFactor'] = float(factor)

    def clearDistRestraint(self):
        '''Deletes any previously defined bond distance restraint(s) for the selected phase

        .. seealso::
           :meth:`G2Phase.addDistRestraint`

        '''
        bondRestData = self._getBondRest('clearDistRestraint')
        bondRestData['Bonds'] = []

    def addDistRestraint(self, origin, target, bond, factor=1.1, ESD=0.01):
        '''Adds bond distance restraint(s) for the selected phase

        This works by search for interatomic distances between atoms in the
        origin list and the target list (the two lists may be the same but most
        frequently will not) with a length between bond/factor and bond*factor.
        If a distance is found in that range, it is added to the restraints
        if it was not already found.

        :param list origin: a list of atoms, each atom may be an atom
           object, an index or an atom label
        :param list target: a list of atoms, each atom may be an atom
           object, an index or an atom label
        :param float bond: the target bond length in A for the located atom
        :param float factor: a tolerance factor used when searching for
           bonds (defaults to 1.1)
        :param float ESD: the uncertainty for the bond (defaults to 0.01)

        :returns: returns the number of new restraints that are found

        As an example::

            gpx = G2sc.G2Project('restr.gpx')
            ph = gpx.phases()[0]
            ph.clearDistRestraint()
            origin = [a for a in ph.atoms() if a.element == 'Si']
            target = [i for i,a in enumerate(ph.atoms()) if a.element == 'O']
            c = ph.addDistRestraint(origin, target, 1.64)
            print(c,'new restraints found')
            ph.setDistRestraintWeight(1000)
            gpx.save('restr-mod.gpx')

        This example locates the first phase in a project file, clears any previous
        restraints. Then it places restraints on bonds between Si and O atoms at
        1.64 A. Each restraint is weighted 1000 times in comparison to
        (obs-calc)/sigma for a data point. To show how atom selection can
        work, the origin atoms are identified here
        by atom object while the target atoms are identified by atom index.
        The methods are interchangeable. If atom labels are unique, then::

           origin = [a.label for a in ph.atoms() if a.element == 'Si']

        would also work identically.

        '''
        bondRestData = self._getBondRest('addDistRestraint')
        originList = []
        for a in origin:
            if type(a) is G2AtomRecord:
                ao = a
            elif type(a) is int:
                ao = self.atoms()[a]
            elif type(a) is str:
                ao = self.atom(a)
                if ao is None:
                    raise G2ScriptException(
                        f'addDistRestraint error: Origin atom matching label "{a}" not found')
            else:
                raise G2ScriptException(
                    f'addDistRestraint error: Origin atom input {a} has unknown type ({type(a)})')
            originList.append([ao.ranId, ao.element, list(ao.coordinates)])
        targetList = []
        for a in target:
            if type(a) is G2AtomRecord:
                ao = a
            elif type(a) is int:
                ao = self.atoms()[a]
            elif type(a) is str:
                ao = self.atom(a)
                if ao is None:
                    raise G2ScriptException(
                        f'addDistRestraint error: Target atom matching label "{a}" not found')
            else:
                raise G2ScriptException(
                    f'addDistRestraint error: Target atom input {a} has unknown type ({type(a)})')
            targetList.append([ao.ranId, ao.element, list(ao.coordinates)])

        GType = self.data['General']['Type']
        SGData = self.data['General']['SGData']
        Amat,Bmat = G2lat.cell2AB(self.data['General']['Cell'][1:7])
        bondlst = G2mth.searchBondRestr(originList,targetList,bond,factor,GType,SGData,Amat,ESD)
        count = 0
        for newBond in bondlst:
            if newBond not in bondRestData['Bonds']:
                count += 1
                bondRestData['Bonds'].append(newBond)
        return count

class G2SeqRefRes(G2ObjectWrapper):
    '''Wrapper for a Sequential Refinement Results tree entry, containing the
    results for a refinement

    Scripts should not try to create a :class:`G2SeqRefRes` object directly as
    this object will be created when a .gpx project file is read.

    As an example::

        import os
        PathWrap = lambda fil: os.path.join('/Users/toby/Scratch/SeqTut2019Mar',fil)
        gpx = G2sc.G2Project(PathWrap('scr4.gpx'))
        seq = gpx.seqref()
        lbl = ('a','b','c','alpha','beta','gamma','Volume')
        for j,h in enumerate(seq.histograms()):
            cell,cellU,uniq = seq.get_cell_and_esd(1,h)
            print(h)
            print([cell[i] for i in list(uniq)+[6]])
            print([cellU[i] for i in list(uniq)+[6]])
            print('')
        print('printed',[lbl[i] for i in list(uniq)+[6]])

    .. seealso::
        :meth:`G2Project.seqref`
    '''
    def __init__(self, data, proj):
        self.data = data
        self.proj = proj
        self.newCellDict = {} # dict with recp. cell tensor & Dij name
        # newAtomDict = {} # dict with atom positions; relative & absolute
        for name in self.data['histNames']:
            self.newCellDict.update(self.data[name].get('newCellDict',{}))
            #newAtomDict.update(self.data[name].get('newAtomDict',{}) # dict with atom positions; relative & absolute
        #ESDlookup = {self.newCellDict[item][0]:item for item in self.newCellDict}
        #Dlookup = {item:self.newCellDict[item][0] for item in self.newCellDict}
        # Possible error: the next might need to be data[histNames[0]]['varyList']
        #atomLookup = {newAtomDict[item][0]:item for item in newAtomDict if item in self.data['varyList']}
        #Dlookup.update({atomLookup[parm]:parm for parm in atomLookup}
        #ESDlookup.update({parm:atomLookup[parm] for parm in atomLookup})

#    @property
    def histograms(self):
        '''returns a list of histograms in the squential fit
        '''
        return self.data['histNames']

    def get_cell_and_esd(self,phase,hist):
        '''Returns a vector of cell lengths and esd values

        :param phase: A phase, which may be specified as a phase object
          (see :class:`G2Phase`), the phase name (str) or the index number (int)
          of the phase in the project, numbered starting from 0.
        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered as in in the project tree
          starting from 0.
        :returns: cell,cellESD,uniqCellIndx where cell (list)
          with the unit cell parameters (a,b,c,alpha,beta,gamma,Volume);
          cellESD are the standard uncertainties on the 7 unit cell
          parameters; and uniqCellIndx is a tuple with indicies for the
          unique (non-symmetry determined) unit parameters (e.g.
          [0,2] for a,c in a tetragonal cell)
        '''
        def striphist(var,insChar=''):
            'strip a histogram number from a var name'
            sv = var.split(':')
            if len(sv) <= 1: return var
            if sv[1]:
                sv[1] = insChar
            return ':'.join(sv)

        uniqCellLookup = [
        [['m3','m3m'],(0,)],
        [['3R','3mR'],(0,3)],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],(0,2)],
        [['mmm'],(0,1,2)],
        [['2/m'+'a'],(0,1,2,3)],
        [['2/m'+'b'],(0,1,2,4)],
        [['2/m'+'c'],(0,1,2,5)],
        [['-1'],(0,1,2,3,4,5)],
        ]

        seqData,histData = self.RefData(hist)
        hId = histData['data'][0]['hId']
        phasedict = self.proj.phase(phase).data
        pId = phasedict['pId']
        pfx = str(pId)+'::' # prefix for A values from phase
        phfx = '%d:%d:'%(pId,hId) # Dij prefix
        # get unit cell & symmetry for phase
        RecpCellTerms = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        zeroDict = {pfx+'A'+str(i):0.0 for i in range(6)}
        SGdata = phasedict['General']['SGData']
        # determine the cell items not defined by symmetry
        laue = SGdata['SGLaue'][:]
        if laue == '2/m':
            laue += SGdata['SGUniq']
        for symlist,celllist in uniqCellLookup:
            if laue in symlist:
                uniqCellIndx = celllist
                break
        else: # should not happen
            uniqCellIndx = list(range(6))
        initialCell = {}
        for i in uniqCellIndx:
            initialCell[str(pId)+'::A'+str(i)] =  RecpCellTerms[i]

        esdLookUp = {}
        dLookup = {}
        # note that varyList keys are p:h:Dij while newCellDict keys are p::Dij
        for nKey in seqData['newCellDict']:
            p = nKey.find('::')+1
            vKey = nKey[:p] + str(hId) + nKey[p:]
            if vKey in seqData['varyList']:
                esdLookUp[self.newCellDict[nKey][0]] = nKey
                dLookup[nKey] = self.newCellDict[nKey][0]
        covData = {'varyList': [dLookup.get(striphist(v),v) for v in seqData['varyList']],
                'covMatrix': seqData['covMatrix']}
        A = RecpCellTerms[:] # make copy of starting A values

        for i,j in enumerate(('D11','D22','D33','D12','D13','D23')):
            var = pfx+'A'+str(i)
            Dvar = phfx+j
            # apply Dij value if non-zero
            if Dvar in seqData['parmDict']:
                A[i] += seqData['parmDict'][Dvar]
            # override with fit result if is Dij varied
            try:
                A[i] = seqData['newCellDict'][esdLookUp[var]][1] # get refined value
            except KeyError:
                pass
        Albls = [pfx+'A'+str(i) for i in range(6)]
        cellDict = dict(zip(Albls,A))


        A,zeros = G2stIO.cellFill(pfx,SGdata,cellDict,zeroDict)
        # convert to direct cell
        c = G2lat.A2cell(A)
        vol = G2lat.calc_V(A)
        cE = G2lat.getCellEsd(pfx,SGdata,A,covData)
        return list(c)+[vol],cE,uniqCellIndx

    def get_VaryList(self,hist):
        '''Returns a list of the refined variables in the
        last refinement cycle for the selected histogram

        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered starting from 0.
        :returns: a list of variables or None if no refinement has been
          performed.
        '''
        try:
            seqData,histData = self.RefData(hist)
            return seqData['varyList']
        except:
            return

    def get_ParmList(self,hist):
        '''Returns a list of all the parameters defined in the
        last refinement cycle for the selected histogram

        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered as in the project tree
          starting from 0.
        :returns: a list of parameters or None if no refinement has been
          performed.
        '''
        try:
            seqData,histData = self.RefData(hist)
            return list(seqData['parmDict'].keys())
        except:
            return

    def get_Variable(self,hist,var):
        '''Returns the value and standard uncertainty (esd) for a variable
        parameters, as defined for the selected histogram
        in the last sequential refinement cycle

        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered as in the project tree
          starting from 0.
        :param str var: a variable name of form '<p>:<h>:<name>', such as
          ':0:Scale'
        :returns: (value,esd) if the parameter is refined or
          (value, None) if the variable is in a constraint or is not
          refined or None if the parameter is not found.
        '''
        try:
            seqData,histData = self.RefData(hist)
            val = seqData['parmDict'][var]
        except:
            return
        try:
            pos = seqData['varyList'].index(var)
            esd = np.sqrt(seqData['covMatrix'][pos,pos])
            return (val,esd)
        except ValueError:
            return (val,None)

    def get_Covariance(self,hist,varList):
        '''Returns the values and covariance matrix for a series of variable
        parameters, as defined for the selected histogram
        in the last sequential refinement cycle

        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered as in the project tree
          starting from 0.
        :param tuple varList: a list of variable names of form '<p>:<h>:<name>'
        :returns: (valueList,CovMatrix) where valueList contains the (n) values
          in the same order as varList (also length n) and CovMatrix is a
          (n x n) matrix. If any variable name is not found in the varyList
          then None is returned.

        Use this code, where sig provides standard uncertainties for
        parameters and where covArray provides the correlation between
        off-diagonal terms::

            sig = np.sqrt(np.diag(covMatrix))
            xvar = np.outer(sig,np.ones_like(sig))
            covArray = np.divide(np.divide(covMatrix,xvar),xvar.T)

        '''
        try:
            seqData,histData = self.RefData(hist)
        except:
            G2fil.G2Print('Warning: Histogram {} not found in the sequential fit'.format(hist))
            return
        missing = [i for i in varList if i not in seqData['varyList']]
        if missing:
            G2fil.G2Print('Warning: Variable(s) {} were not found in the varyList'.format(missing))
            return None
        vals = [seqData['parmDict'][i] for i in varList]
        cov = G2mth.getVCov(varList,seqData['varyList'],seqData['covMatrix'])
        return (vals,cov)

    def RefData(self,hist):
        '''Provides access to the output from a particular histogram

        :param hist: Specify a histogram or using the histogram name (str)
          or the index number (int) of the histogram in the sequential
          refinement (not the project), numbered as in the project tree
          starting from 0.
        :returns: a list of dicts where the first element has sequential
          refinement results and the second element has the contents of
          the histogram tree items.
        '''
        errmsg = ''
        try:
            hist = self.data['histNames'][hist]
        except IndexError:
            errmsg = 'Histogram #{hist} is out of range from the Sequential Refinement'
        except TypeError:
            pass
        if errmsg: raise Exception(errmsg)
        if hist not in self.data['histNames']:
            raise Exception('Histogram {} is not included in the Sequential Refinement'
                                .format(hist))
        return self.data[hist],self.proj.histogram(hist).data

class G2PDF(G2ObjectWrapper):
    """Wrapper for a PDF tree entry, containing the information needed to
    compute a PDF and the S(Q), G(r) etc. after the computation is done.
    Note that in a GSASIIscriptable script, instances of G2PDF will be created by
    calls to :meth:`G2Project.add_PDF` or :meth:`G2Project.pdf`.
    Scripts should not try to create a :class:`G2PDF` object directly.

    Example use of :class:`G2PDF`::

       gpx.add_PDF('250umSiO2.pdfprm',0)
       pdf.set_formula(['Si',1],['O',2])
       pdf.set_background('Container',1,-0.21)
       for i in range(5):
           if pdf.optimize(): break
       pdf.calculate()
       pdf.export(gpx.filename,'S(Q), pdfGUI')
       gpx.save('pdfcalc.gpx')

    .. seealso::
        :meth:`G2Project.pdf`
        :meth:`G2Project.pdfs`
    """
    def __init__(self, data, name, proj):
        self.data = data
        self.name = name
        self.proj = proj
    def set_background(self,btype,histogram,mult=-1.,refine=False):
        '''Sets a histogram to be used as the 'Sample Background',
        the 'Container' or the 'Container Background.'

        :param str btype: Type of background to set, must contain
          the string 'samp' for Sample Background', 'cont' and 'back'
          for the 'Container Background' or only 'cont' for the
          'Container'. Note that capitalization and extra characters
          are ignored, so the full strings (such as 'Sample
          Background' & 'Container Background') can be used.

        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :param float mult: a multiplier for the histogram; defaults
          to -1.0
        :param bool refine: a flag to enable refinement (only
          implemented for 'Sample Background'); defaults to False
        '''
        if 'samp' in btype.lower():
            key = 'Sample Bkg.'
        elif 'cont' in btype.lower() and 'back' in btype.lower():
            key = 'Container Bkg.'
        elif 'cont' in btype.lower():
            key = 'Container'
        else:
            raise Exception('btype = {} is invalid'.format(btype))
        self.data['PDF Controls'][key]['Name'] = self.proj.histogram(histogram).name
        self.data['PDF Controls'][key]['Mult'] = mult
        self.data['PDF Controls'][key]['Refine'] = refine

    def set_formula(self,*args):
        '''Set the chemical formula for the PDF computation.
        Use pdf.set_formula(['Si',1],['O',2]) for SiO2.

        :param list item1: The element symbol and number of atoms in formula for first element
        :param list item2: The element symbol and number of atoms in formula for second element,...

        repeat parameters as needed for all elements in the formula.
        '''
        powderHist = self.proj.histogram(self.data['PDF Controls']['Sample']['Name'])
        inst = powderHist.data['Instrument Parameters'][0]
        ElList = self.data['PDF Controls']['ElList']
        ElList.clear()
        sumVol = 0.
        for elem,mult in args:
            ElList[elem] = G2elem.GetElInfo(elem,inst)
            ElList[elem]['FormulaNo'] = mult
            Avol = (4.*np.pi/3.)*ElList[elem]['Drad']**3
            sumVol += Avol*ElList[elem]['FormulaNo']
        self.data['PDF Controls']['Form Vol'] = max(10.0,sumVol)

    def calculate(self,xydata=None,limits=None,inst=None):
        '''Compute the PDF using the current parameters. Results are set
        in the PDF object arrays (self.data['PDF Controls']['G(R)'] etc.).
        Note that if ``xydata``, is specified, the background histograms(s)
        will not be accessed from the project file associated with the current
        PDF entry. If ``limits`` and ``inst`` are both specified, no histograms
        need be in the current project. However, the self.data['PDF Controls']
        sections ('Sample', 'Sample Bkg.','Container Bkg.') must be
        non-blank for the corresponding items to be used from``xydata``.

        :param dict xydata: an array containing the Sample's I vs Q, and
          any or none of the Sample Background, the Container scattering and
          the Container Background. If xydata is None (default), the values are
          taken from histograms, as named in the PDF's self.data['PDF Controls']
          entries with keys 'Sample', 'Sample Bkg.','Container Bkg.' &
          'Container'.
        :param list limits: upper and lower Q values to be used for PDF
          computation. If None (default), the values are
          taken from the Sample histogram's .data['Limits'][1] values.
        :param dict inst: The Sample histogram's instrument parameters
          to be used for PDF computation. If None (default), the values are
          taken from the Sample histogram's .data['Instrument Parameters'][0]
          values.
        '''
        data = self.data['PDF Controls']
        if xydata is None:
            xydata = {}
            for key in 'Sample Bkg.','Container Bkg.','Container','Sample':
                name = data[key]['Name'].strip()
                if name:
                    xydata[key] = self.proj.histogram(name).data['data']
        if limits is None:
            name = data['Sample']['Name'].strip()
            limits = self.proj.histogram(name).data['Limits'][1]
        if inst is None:
            name = data['Sample']['Name'].strip()
            inst = self.proj.histogram(name).data['Instrument Parameters'][0]
        G2pwd.CalcPDF(data,inst,limits,xydata)
        data['I(Q)'] = xydata['IofQ']
        data['S(Q)'] = xydata['SofQ']
        data['F(Q)'] = xydata['FofQ']
        data['G(R)'] = xydata['GofR']
        data['g(r)'] = xydata['gofr']

    def optimize(self,showFit=True,maxCycles=5,
                     xydata=None,limits=None,inst=None):
        '''Optimize the low R portion of G(R) to minimize selected
        parameters. Note that this updates the parameters in the settings
        (self.data['PDF Controls']) but does not update the PDF object
        arrays (self.data['PDF Controls']['G(R)'] etc.) with the computed
        values, use :meth:`calculate` after a fit to do that.

        :param bool showFit: if True (default) the optimized parameters
          are shown before and after the fit, as well as the RMS value
          in the minimized region.
        :param int maxCycles: the maximum number of least-squares cycles;
          defaults to 5.
        :returns: the result from the optimizer as True or False, depending
          on if the refinement converged.
        :param dict xydata: an array containing the Sample's I vs Q, and
          any or none of the Sample Background, the Container scattering and
          the Container Background. If xydata is None (default), the values are
          taken from histograms, as named in the PDF's self.data['PDF Controls']
          entries with keys 'Sample', 'Sample Bkg.','Container Bkg.' &
          'Container'.
        :param list limits: upper and lower Q values to be used for PDF
          computation. If None (default), the values are
          taken from the Sample histogram's .data['Limits'][1] values.
        :param dict inst: The Sample histogram's instrument parameters
          to be used for PDF computation. If None (default), the values are
          taken from the Sample histogram's .data['Instrument Parameters'][0]
          values.
        '''
        data = self.data['PDF Controls']
        if xydata is None:
            xydata = {}
            for key in 'Sample Bkg.','Container Bkg.','Container','Sample':
                name = data[key]['Name'].strip()
                if name:
                    xydata[key] = self.proj.histogram(name).data['data']
        if limits is None:
            name = data['Sample']['Name'].strip()
            limits = self.proj.histogram(name).data['Limits'][1]
        if inst is None:
            name = data['Sample']['Name'].strip()
            inst = self.proj.histogram(name).data['Instrument Parameters'][0]

        res = G2pwd.OptimizePDF(data,xydata,limits,inst,showFit,maxCycles)
        return res['success']

    def export(self,fileroot,formats):
        '''Write out the PDF-related data (G(r), S(Q),...) into files

        :param str fileroot: name of file(s) to be written. The extension
          will be ignored and set to .iq, .sq, .fq or .gr depending
          on the formats selected.
        :param str formats: string specifying the file format(s) to be written,
          should contain at least one of the following keywords:
          I(Q), S(Q), F(Q), G(r) and/or PDFgui (capitalization and
          punctuation is ignored). Note that G(r) and PDFgui should not
          be specifed together.
        '''
        PDFsaves = 5*[False]
        PDFentry = self.name
        name = self.data['PDF Controls']['Sample']['Name'].strip()
        limits = self.proj.histogram(name).data['Limits']
        inst = self.proj.histogram(name).data['Instrument Parameters'][0]
        for i,lbl in enumerate(['I(Q)', 'S(Q)', 'F(Q)', 'G(r)', 'PDFgui']):
            PDFsaves[i] = lbl.lower() in formats.lower()
        G2fil.PDFWrite(PDFentry,fileroot,PDFsaves,self.data['PDF Controls'],inst,limits)

class G2Single(G2ObjectWrapper):
    """Wrapper for a HKLF tree entry, containing a single crystal histogram
    Note that in a GSASIIscriptable script, instances of G2Single will be
    created by calls to :meth:`G2Project.histogram`,
    :meth:`G2Project.histograms`, or :meth:`G2Project.add_single_histogram`.
    Scripts should not try to create a :class:`G2Single` object directly.

    This object contains these class variables:
        * G2Single.proj: contains a reference to the :class:`G2Project`
          object that contains this histogram
        * G2Single.name: contains the name of the histogram
        * G2Single.data: contains the histogram's associated data in a dict,
          as documented for the :ref:`Single Crystal Tree Item<Xtal_table>`.
          This contains the actual histogram values, as documented for Data.

    Example use of :class:`G2Single`::

        gpx0 = G2sc.G2Project('HTO_base.gpx')
        gpx0.add_single_histogram('HTO_xray/xtal1/xs2555a.hkl',0,fmthint='Shelx HKLF 4')
        gpx0.save('HTO_scripted.gpx')

    This opens an existing GSAS-II project file and adds a single
    crystal dataset that is linked to the first phase and saves it
    under a new name.

    .. seealso::
        :meth:`~G2Project.add_single_histogram`
        :meth:`~G2Project.histogram`
        :meth:`~G2Project.histograms`
        :meth:`~G2Project.link_histogram_phase`
    """
    def __init__(self, data, proj, name):
        self.data = data
        self.name = name
        self.proj = proj

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ["Single xtal"]
        return key in valid_keys

    def set_refinements(self, refs):
        """Sets the HKLF histogram refinement parameter 'key' to the
        specification 'value'.

        :param dict refs: A dictionary of the parameters to be set. See the
                          :ref:`Histogram_parameters_table` table for a description of
                          what these dictionaries should be.

        Example::

            hist.set_refinements({'Scale':True,'Es':False,'Flack':True})
        """
        # find the phase associated with this histogram (there is at most one)
        for p in self.proj.phases():
            if self.name in p.data['Histograms']:
                phase = p
                break
        else:  # no phase matches
            raise ValueError(f'set_refinements error: no phase found for {self.name!r}')
        d = phase.data['Histograms'][self.name]
        for key, value in refs.items():
            if key in d:  # Scale & Flack
                d[key][1] = value
            elif key in d['Babinet']:  # BabA & BabU
                d['Babinet'][key][1] = value
            elif key in d['Extinction'][2]:  # Eg, Es, Ep
                d['Extinction'][2][key][1] = value
            else:
                raise ValueError(f'set_refinements: {key} is unknown')

    def clear_refinements(self, refs):
        """Clears the HKLF refinement parameter 'key' and its associated value.

        :param dict refs: A dictionary of parameters to clear.
          See the :ref:`Histogram_parameters_table` table for what can be specified.

        Example::

            hist.clear_refinements(['Scale','Es','Flack'])
            hist.clear_refinements({'Scale':True,'Es':False,'Flack':True})

        Note that the two above commands are equivalent: the values specified
        in the dict in the second command are ignored.
        """
        # find the phase associated with this histogram (there is at most one)
        for p in self.proj.phases():
            if self.name in p.data['Histograms']:
                phase = p
                break
        else:  # no phase matches
            raise ValueError(f'clear_refinements error: no phase found for {self.name!r}')
        d = phase.data['Histograms'][self.name]
        for key in refs:
            if key in d:  # Scale & Flack
                d[key][1] = False
            elif key in d['Babinet']:  # BabA & BabU
                d['Babinet'][key][1] = False
            elif key in d['Extinction'][2]:  # Eg, Es, Ep
                d['Extinction'][2][key][1] = False
            else:
                raise ValueError(f'clear_refinements: {key} is unknown')

    def Export(self,fileroot,extension,fmthint=None):
        '''Write the HKLF histogram into a file. The path is specified by
        fileroot and extension.

        :param str fileroot: name of the file, optionally with a path (extension is
           ignored)
        :param str extension: includes '.', must match an extension in global
           exportersByExtension['single'] or a Exception is raised.
        :param str fmthint: If specified, the first exporter where the format
           name (obj.formatName, as shown in Export menu) contains the
           supplied string will be used. If not specified, an error
           will be generated showing the possible choices.
        :returns: name of file that was written
        '''
        LoadG2fil()
        if extension not in exportersByExtension.get('single',[]):
            print('Defined exporters are:')
            print('  ',list(exportersByExtension.get('single',[])))
            raise G2ScriptException('No Writer for file type = "'+extension+'"')
        fil = os.path.abspath(os.path.splitext(fileroot)[0]+extension)
        obj = exportersByExtension['single'][extension]
        if type(obj) is list:
            if fmthint is None:
                print('Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('No format hint for file type = "'+extension+'"')
            for o in obj:
              if fmthint.lower() in o.formatName.lower():
                  obj = o
                  break
            else:
                print('Hint ',fmthint,'not found. Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('Bad format hint for file type = "'+extension+'"')
        #self._SetFromArray(obj)
        obj.Writer(self,fil)
        return fil

blkSize = 128
'''Integration block size; 128 or 256 seems to be optimal for CPU use, but 128 uses
less memory, must be <=1024 (for polymask/histogram3d)
'''

def calcMaskMap(imgprms,mskprms):
    '''Computes a set of blocked mask arrays for a set of image controls and mask parameters.
    This capability is also provided with :meth:`G2Image.IntMaskMap`.
    '''
    return G2img.MakeUseMask(imgprms,mskprms,blkSize)

def calcThetaAzimMap(imgprms):
    '''Computes the set of blocked arrays for theta-azimuth mapping from
    a set of image controls, which can be cached and reused for
    integration of multiple images with the same calibration parameters.
    This capability is also provided with :meth:`G2Image.IntThetaAzMap`.
    '''
    return G2img.MakeUseTA(imgprms,blkSize)

class G2Image(G2ObjectWrapper):
    '''Wrapper for an IMG tree entry, containing an image and associated metadata.

    Note that in a GSASIIscriptable script, instances of G2Image will be created by
    calls to :meth:`G2Project.add_image` or :meth:`G2Project.images`.
    Scripts should not try to create a :class:`G2Image` object directly as
    :meth:`G2Image.__init__` should be invoked from inside :class:`G2Project`.

    The object contains these class variables:

        * G2Image.proj: contains a reference to the :class:`G2Project`
          object that contains this image
        * G2Image.name: contains the name of the image
        * G2Image.data: contains the image's associated data in a dict,
          as documented for the :ref:`Image Data Structure<Image_table>`.
        * G2Image.image: optionally contains a cached the image to
          save time in reloading. This is saved only when cacheImage=True
          is specified when :meth:`G2Project.add_image` is called.

    Example use of G2Image:

    >>> gpx = G2sc.G2Project(newgpx='itest.gpx')
    >>> imlst = gpx.add_image(idata,fmthint="TIF")
    >>> imlst[0].loadControls('stdSettings.imctrl')
    >>> imlst[0].setCalibrant('Si    SRM640c')
    >>> imlst[0].loadMasks('stdMasks.immask')
    >>> imlst[0].Recalibrate()
    >>> imlst[0].setControl('outAzimuths',3)
    >>> pwdrList = imlst[0].Integrate()

    More detailed image processing examples are shown in the
    :ref:`ImageProc` section of this chapter.

    '''
    # parameters in that can be accessed via setControl. This may need future attention
    ControlList = {
        'int': ['calibskip', 'pixLimit', 'edgemin', 'outChannels',
                    'outAzimuths'],
        'float': ['cutoff', 'setdist', 'wavelength', 'Flat Bkg',
                      'azmthOff', 'tilt', 'calibdmin', 'rotation',
                      'distance', 'DetDepth'],
        'bool': ['setRings', 'setDefault', 'centerAzm', 'fullIntegrate',
                     'DetDepthRef', 'showLines'],
        'str': ['SampleShape', 'binType', 'formatName', 'color',
                    'type', ],
        'list': ['GonioAngles', 'IOtth', 'LRazimuth', 'Oblique', 'PolaVal',
                   'SampleAbs', 'center', 'ellipses', 'linescan',
                    'pixelSize', 'range', 'ring', 'rings', 'size', ],
        'dict': ['varyList'],
        }
    '''Defines the items known to exist in the Image Controls tree section
    and the item's data types. A few are not included here
    ('background image', 'dark image', 'Gain map', and 'calibrant') because
    these items have special set routines,
    where references to entries are checked to make sure their values are
    correct.
    '''

    def __init__(self, data, name, proj, image=None):
        self.data = data
        self.name = name
        self.proj = proj
        self.image = image

    def clearImageCache(self):
        '''Clears a cached image, if one is present
        '''
        self.image = None

    def setControl(self,arg,value):
        '''Set an Image Controls parameter in the current image.
        If the parameter is not found an exception is raised.

        :param str arg: the name of a parameter (dict entry) in the
          image. The parameter must be found in :data:`ControlList`
          or an exception is raised.
        :param value: the value to set the parameter. The value is
          cast as the appropriate type from :data:`ControlList`.
        '''
        for typ in self.ControlList:
            if arg in self.ControlList[typ]: break
        else:
            G2fil.G2Print('Allowed args:\n',[nam for nam,typ in self.findControl('')])
            raise Exception('arg {} not defined in G2Image.setControl'
                                .format(arg))
        errmsg = ''
        try:
            if typ == 'int':
                self.data['Image Controls'][arg] = int(value)
            elif typ == 'float':
                self.data['Image Controls'][arg] = float(value)
            elif typ == 'bool':
                self.data['Image Controls'][arg] = bool(value)
            elif typ == 'str':
                self.data['Image Controls'][arg] = str(value)
            elif typ == 'list':
                self.data['Image Controls'][arg] = list(value)
            elif typ == 'dict':
                self.data['Image Controls'][arg] = dict(value)
            else:
                raise Exception('Unknown type {} for arg {} in  G2Image.setControl'
                                    .format(typ,arg))
        except:
            errmsg = 'Error formatting value {value} as type {typ} for arg {arg} in  G2Image.setControl'

        if errmsg: raise Exception(errmsg)

    def getControl(self,arg):
        '''Return an Image Controls parameter in the current image.
        If the parameter is not found an exception is raised.

        :param str arg: the name of a parameter (dict entry) in the
          image.
        :returns: the value as a int, float, list,...
        '''
        if arg in self.data['Image Controls']:
            return self.data['Image Controls'][arg]
        G2fil.G2Print(self.findControl(''))
        raise Exception('arg {} not defined in G2Image.getControl'.format(arg))

    def findControl(self,arg=''):
        '''Finds the Image Controls parameter(s) in the current image
        that match the string in arg. Default is '' which returns all
        parameters.

            Example:

            >>> findControl('calib')
            [['calibskip', 'int'], ['calibdmin', 'float'], ['calibrant', 'str']]

        :param str arg: a string containing part of the name of a
          parameter (dict entry) in the image's Image Controls.
        :returns: a list of matching entries in form
          [['item','type'], ['item','type'],...] where each 'item' string
          contains the sting in arg.
        '''
        matchList = []
        for typ in self.ControlList:
            for item in self.ControlList[typ]:
                if arg in item:
                    matchList.append([item,typ])
        return matchList

    def setCalibrant(self,calib):
        '''Set a calibrant for the current image

        :param str calib: specifies a calibrant name which must be one of
          the entries in file ImageCalibrants.py. This is validated and
          an error provides a list of valid choices.
        '''
        from . import ImageCalibrants as calFile
        if calib in calFile.Calibrants.keys():
            self.data['Image Controls']['calibrant'] = calib
            return
        G2fil.G2Print('Calibrant {} is not valid. Valid calibrants'.format(calib))
        for i in calFile.Calibrants.keys():
            if i: G2fil.G2Print('\t"{}"'.format(i))

    def setControlFile(self,typ,imageRef,mult=None):
        '''Set a image to be used as a background/dark/gain map image

        :param str typ: specifies image type, which must be one of:
           'background image', 'dark image', 'gain map'; N.B. only the first
           four characters must be specified and case is ignored.
        :param imageRef: A reference to the desired image. Either the Image
          tree name (str), the image's index (int) or
          a image object (:class:`G2Image`)
        :param float mult: a multiplier to be applied to the image (not used
          for 'Gain map'; required for 'background image', 'dark image'
        '''
        if 'back' in typ.lower():
            key = 'background image'
            mult = float(mult)
        elif 'dark' in typ.lower():
            key = 'dark image'
            mult = float(mult)
        elif 'gain' in typ.lower():
            #key = 'Gain map'
            if mult is not None:
                G2fil.G2Print('Warning: Ignoring multiplier for Gain map')
            mult = None
        else:
            raise Exception("Invalid typ {} for setControlFile".format(typ))
        imgNam = self.proj.image(imageRef).name
        if mult is None:
            self.data['Image Controls']['Gain map'] = imgNam
        else:
            self.data['Image Controls'][key] = [imgNam,mult]

    def loadControls(self,filename=None,imgDict=None):
        '''load controls from a .imctrl file

        :param str filename: specifies a file to be read, which should end
          with .imctrl (defaults to None, meaning parameters are input
          with imgDict.)
        :param dict imgDict: contains a set of image parameters (defaults to
          None, meaning parameters are input with filename.)
        '''
        if filename:
            File = open(filename,'r')
            Slines = File.readlines()
            File.close()
            G2fil.LoadControls(Slines,self.data['Image Controls'])
            G2fil.G2Print('file {} read into {}'.format(filename,self.name))
        elif imgDict:
            self.data['Image Controls'].update(imgDict)
            G2fil.G2Print('Image controls set')
        else:
            raise Exception("loadControls called without imgDict or filename specified")

    def saveControls(self,filename):
        '''write current controls values to a .imctrl file

        :param str filename: specifies a file to write, which should end
          with .imctrl
        '''
        G2fil.WriteControls(filename,self.data['Image Controls'])
        G2fil.G2Print('file {} written from {}'.format(filename,self.name))

    def getControls(self,clean=False):
        '''returns current Image Controls as a dict

        :param bool clean: causes the calbration information to be deleted
        '''
        ImageControls = copy.deepcopy(self.data['Image Controls'])
        if clean:
            ImageControls['showLines'] = True
            ImageControls['ring'] = []
            ImageControls['rings'] = []
            ImageControls['ellipses'] = []
            ImageControls['setDefault'] = False
            for i in 'range','size','GonioAngles':
                if i in ImageControls: del ImageControls[i]
        return ImageControls

    def setControls(self,controlsDict):
        '''uses dict from :meth:`getControls` to set Image Controls for current image
        '''
        self.data['Image Controls'].update(copy.deepcopy(controlsDict))


    def loadMasks(self,filename,ignoreThreshold=False):
        '''load masks from a .immask file

        :param str filename: specifies a file to be read, which should end
          with .immask
        :param bool ignoreThreshold: If True, masks are loaded with
          threshold masks. Default is False which means any Thresholds
          in the file are ignored.
        '''
        G2fil.readMasks(filename,self.data['Masks'],ignoreThreshold)
        G2fil.G2Print('file {} read into {}'.format(filename,self.name))

    def initMasks(self):
        '''Initialize Masks, including resetting the Thresholds values
        '''
        self.data['Masks'] = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[]}
        if self.image is not None:
            ImageZ = self.image
        else:
            ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        Imin = max(0.,np.min(ImageZ))
        Imax = np.max(ImageZ)
        self.data['Masks']['Thresholds'] = [(0,Imax),[Imin,Imax]]

    def getMasks(self):
        '''load masks from an IMG tree entry
        '''
        return self.data['Masks']

    def setMasks(self,maskDict,resetThresholds=False):
        '''load masks dict (from :meth:`getMasks`) into current IMG record

        :param dict maskDict: specifies a dict with image parameters,
          from :meth:`getMasks`
        :param bool resetThresholds: If True, Threshold Masks in the
          dict are ignored. The default is False which means Threshold
          Masks are retained.
        '''
        self.data['Masks'] = copy.deepcopy(maskDict)
        if resetThresholds:
            if self.image is not None:
                ImageZ = self.image
            else:
                ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
            Imin = max(0.,np.min(ImageZ))
            Imax = np.max(ImageZ)
            self.data['Masks']['Thresholds'] = [(0,Imax),[Imin,Imax]]

    def IntThetaAzMap(self):
        '''Computes the set of blocked arrays for 2theta-azimuth mapping from
        the controls settings of the current image for image integration.
        The output from this is optionally supplied as input to
        :meth:`~G2Image.Integrate`. Note that if not supplied, image
        integration will compute this information as it is needed, but this
        is a relatively slow computation so time can be saved by caching and
        reusing this computation for other images that have the
        same calibration parameters as the current image.
        '''
        return G2img.MakeUseTA(self.getControls(),blkSize)

    def IntMaskMap(self):
        '''Computes a series of masking arrays for the current image (based on
        mask input, but not calibration parameters or the image intensities).
        See :meth:`GSASIIimage.MakeMaskMap` for more details. The output from
        this is optionally supplied as input to :meth:`~G2Image.Integrate`).

        Note this is not the same as pixel mask
        searching (:meth:`~G2Image.GeneratePixelMask`).
        '''
        return G2img.MakeUseMask(self.getControls(),self.getMasks(),blkSize)

    def MaskThetaMap(self):
        '''Computes the theta mapping matrix from the controls settings
        of the current image to be used for pixel mask computation
        in :meth:`~G2Image.GeneratePixelMask`.
        This is optional, as if not supplied, mask computation will compute
        this, but this is a relatively slow computation and the
        results computed here can be reused for other images that have the
        same calibration parameters.
        '''
        ImShape = self.getControls()['size']
        return G2img.Make2ThetaAzimuthMap(self.getControls(), (0, ImShape[0]), (0, ImShape[1]))[0]

    def MaskFrameMask(self):
        '''Computes a Frame mask from map input for the current image to be
        used for a pixel mask computation in
        :meth:`~G2Image.GeneratePixelMask`.
        This is optional, as if not supplied, mask computation will compute
        this, but this is a relatively slow computation and the
        results computed here can be reused for other images that have the
        same calibration parameters.
        '''
        Controls = self.getControls()
        Masks = self.getMasks()
        frame = Masks['Frames']
        if self.image is not None:
            ImageZ = self.image
        else:
            ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        tam = ma.make_mask_none(ImageZ.shape)
        if frame:
            tam = ma.mask_or(tam,ma.make_mask(np.abs(G2img.polymask(Controls,frame)-255)))
        return tam

    def getVary(self,*args):
        '''Return the refinement flag(s) for calibration of
        Image Controls parameter(s) in the current image.
        If the parameter is not found, an exception is raised.

        :param str arg: the name of a refinement parameter in the
          varyList for the image. The name should be one of
          'dep', 'det-X', 'det-Y', 'dist', 'phi', 'tilt', or 'wave'
        :param str arg1: the name of a parameter (dict entry) as before,
          optional

        :returns: a list of bool value(s)
        '''
        res = []
        for arg in args:
            if arg in self.data['Image Controls']['varyList']:
                res.append(self.data['Image Controls']['varyList'][arg])
            else:
                raise Exception('arg {} not defined in G2Image.getVary'.format(arg))
        return res

    def setVary(self,arg,value):
        '''Set a refinement flag for Image Controls parameter in the
        current image that is used for fitting calibration parameters.
        If the parameter is not '*' or found, an exception is raised.

        :param str arg: the name of a refinement parameter in the
          varyList for the image. The name should be one of
          'dep', 'det-X', 'det-Y', 'dist', 'phi', 'tilt', or 'wave',
          or it may be a list or tuple of names,
          or it may be '*' in which all parameters are set accordingly.
        :param value: the value to set the parameter. The value is
          cast as bool.
        '''
        if arg == '*':
            for a in self.data['Image Controls']['varyList']:
                self.data['Image Controls']['varyList'][a] = bool(value)
            return
        if not isinstance(arg,tuple) and not isinstance(arg,list):
            arg = [arg]
        for a in arg:
            if a in self.data['Image Controls']['varyList']:
                self.data['Image Controls']['varyList'][a] = bool(value)
            else:
                raise Exception('arg {} not defined in G2Image.setVary'.format(a))

    def Recalibrate(self):
        '''Invokes a recalibration fit (same as Image Controls/Calibration/Recalibrate
        menu command). Note that for this to work properly, the calibration
        coefficients (center, wavelength, distance & tilts) must be fairly close.
        This may produce a better result if run more than once.
        '''
        LoadG2fil()
        if self.image is not None:
            ImageZ = self.image
        else:
            ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        G2img.ImageRecalibrate(None,ImageZ,self.data['Image Controls'],self.data['Masks'])

    def Integrate(self,name=None,MaskMap=None,ThetaAzimMap=None):
        '''Invokes an image integration (same as Image Controls/Integration/Integrate
        menu command). All parameters will have previously been set with Image Controls
        so no input is needed here. However, the optional parameters MaskMap
        and ThetaAzimMap may be supplied to save computing these items more than
        once, speeding integration of multiple images with the same
        image/mask parameters.

        Note that if integration is performed on an
        image more than once, histogram entries may be overwritten. Use the name
        parameter to prevent this if desired.

        :param str name: base name for created histogram(s). If None (default),
          the histogram name is taken from the image name.
        :param list MaskMap: from :func:`IntMaskMap`
        :param list ThetaAzimMap: from :meth:`G2Image.IntThetaAzMap`
        :returns: a list of created histogram (:class:`G2PwdrData` or :class:`G2SmallAngle`) objects.
        '''
        if self.image is not None:
            ImageZ = self.image
        else:
            ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        # Apply a Pixel Mask to image, if present
        Masks = self.getMasks()
        if Masks.get('SpotMask',{'spotMask':None})['spotMask'] is not None:
            ImageZ = ma.array(ImageZ,mask=Masks['SpotMask']['spotMask'])
        # do integration
        ints,azms,Xvals,cancel = G2img.ImageIntegrate(ImageZ,
                self.data['Image Controls'],self.data['Masks'],blkSize=blkSize,
                useMask=MaskMap,useTA=ThetaAzimMap)
        # code from here on based on G2IO.SaveIntegration, but places results in the current
        # project rather than tree
        X = Xvals[:-1]
        N = len(X)

        data = self.data['Image Controls']
        Comments = self.data['Comments']
        # make name from image, unless overridden
        if name:
            if not name.startswith(data['type']+' '):
                name = data['type']+' '+name
        else:
            name = self.name.replace('IMG ',data['type']+' ')
        if 'PWDR' in name:
            if 'target' in data:
                names = ['Type','Lam1','Lam2','I(L2)/I(L1)','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
                codes = [0 for i in range(14)]
            else:
                names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
                codes = [0 for i in range(12)]
        elif 'SASD' in name:
            names = ['Type','Lam','Zero','Azimuth']
            codes = [0 for i in range(4)]
            X = 4.*np.pi*npsind(X/2.)/data['wavelength']    #convert to q
        Xminmax = [X[0],X[-1]]
        Azms = []
        dazm = 0.
        if data['fullIntegrate'] and data['outAzimuths'] == 1:
            Azms = [45.0,]                              #a poor man's average?
        else:
            for i,azm in enumerate(azms[:-1]):
                if azm > 360. and azms[i+1] > 360.:
                    Azms.append(G2img.meanAzm(azm%360.,azms[i+1]%360.))
                else:
                    Azms.append(G2img.meanAzm(azm,azms[i+1]))
            dazm = np.min(np.abs(np.diff(azms)))/2.
        # pull out integration results and make histograms for each
        IntgOutList = []
        for i,azm in enumerate(azms[:-1]):
            Aname = name+" Azm= %.2f"%((azm+dazm)%360.)
            # MT dict to contain histogram
            HistDict = {}
            histItems = [Aname]
            Sample = G2obj.SetDefaultSample()       #set as Debye-Scherrer
            Sample['Gonio. radius'] = data['distance']
            Sample['Omega'] = data['GonioAngles'][0]
            Sample['Chi'] = data['GonioAngles'][1]
            Sample['Phi'] = data['GonioAngles'][2]
            Sample['Azimuth'] = (azm+dazm)%360.    #put here as bin center
            polariz = 0.99    #set default polarization for synchrotron radiation!
            for item in Comments:
                if 'polariz' in item:
                    try:
                        polariz = float(item.split('=')[1])
                    except:
                        polariz = 0.99
                for key in ('Temperature','Pressure','Time','FreePrm1','FreePrm2','FreePrm3','Omega',
                    'Chi','Phi'):
                    if key.lower() in item.lower():
                        try:
                            Sample[key] = float(item.split('=')[1])
                        except:
                            pass
                # if 'label_prm' in item.lower():
                #     for num in ('1','2','3'):
                #         if 'label_prm'+num in item.lower():
                #             Controls['FreePrm'+num] = item.split('=')[1].strip()
            if 'PWDR' in Aname:
                if 'target' in data:    #from lab x-ray 2D imaging data
                    waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
                                 'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'MoKa':[0.70926,0.713543],
                                 'AgKa':[0.559363,0.563775]}
                    wave1,wave2 = waves[data['target']]
                    parms = ['PXC',wave1,wave2,0.5,0.0,polariz,290.,-40.,30.,6.,-14.,0.0,0.0001,Azms[i]]
                else:
                    parms = ['PXC',data['wavelength'],0.0,polariz,1.0,-0.10,0.4,0.30,1.0,0.0,0.0001,Azms[i]]
            elif 'SASD' in Aname:
                Sample['Trans'] = data['SampleAbs'][0]
                parms = ['LXC',data['wavelength'],0.0,Azms[i]]
            Y = ints[i]
            Ymin = np.min(Y)
            Ymax = np.max(Y)
            W = np.where(Y>0.,1./Y,1.e-6)                    #probably not true
            section = 'Comments'
            histItems += [section]
            HistDict[section] = Comments
            section = 'Limits'
            histItems += [section]
            HistDict[section] = copy.deepcopy([tuple(Xminmax),Xminmax])
            if 'PWDR' in Aname:
                section = 'Background'
                histItems += [section]
                HistDict[section] = [['chebyschev-1',1,3,1.0,0.0,0.0],
                    {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}]
            inst = [dict(zip(names,zip(parms,parms,codes))),{}]
            for item in inst[0]:
                inst[0][item] = list(inst[0][item])
            section = 'Instrument Parameters'
            histItems += [section]
            HistDict[section] = inst
            if 'PWDR' in Aname:
                section = 'Sample Parameters'
                histItems += [section]
                HistDict[section] = Sample
                section = 'Peak List'
                histItems += [section]
                HistDict[section] = {'sigDict':{},'peaks':[]}
                section = 'Index Peak List'
                histItems += [section]
                HistDict[section] = [[],[]]
                section = 'Unit Cells List'
                histItems += [section]
                HistDict[section] = []
                section = 'Reflection Lists'
                histItems += [section]
                HistDict[section] = {}
            elif 'SASD' in Aname:
                section = 'Substances'
                histItems += [section]
                HistDict[section] = G2pwd.SetDefaultSubstances()
                section = 'Sample Parameters'
                histItems += [section]
                HistDict[section] = Sample
                section = 'Models'
                histItems += [section]
                HistDict[section] = G2pwd.SetDefaultSASDModel()
            valuesdict = { # same for SASD and PWDR
                'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),'Offset':[0.0,0.0],'delOffset':0.02*Ymax,
                'refOffset':-0.1*Ymax,'refDelt':0.1*Ymax,'Yminmax':[Ymin,Ymax]}
            # if Aname is already in the project replace it
            for j in self.proj.names:
                if j[0] == Aname:
                    G2fil.G2Print('Replacing "{}" in project'.format(Aname))
                    break
            else:
                G2fil.G2Print('Adding "{}" to project'.format(Aname))
                self.proj.names.append(histItems)
            HistDict['data'] = [valuesdict,
                [np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]]
            self.proj.data[Aname] = HistDict
            if 'PWDR' in Aname:
                IntgOutList.append(self.proj.histogram(Aname))
            elif 'SASD' in Aname:
                IntgOutList.append(self.proj.SAS(Aname))
            else:
                print(f'Created histogram of unknown type {Aname}')
        return IntgOutList

    def TestFastPixelMask(self):
        '''Tests to see if the fast (C) code for pixel masking is installed.

        :returns: A value of True is returned if fast pixel masking is
          available. Otherwise False is returned.
        '''
        return G2img.TestFastPixelMask()

    def GeneratePixelMask(self,esdMul=3.0,ttmin=0.,ttmax=180.,
                              FrameMask=None,ThetaMap=None,
                              fastmode=True,combineMasks=False):
        '''Generate a Pixel mask with True at the location of pixels that are
        statistical outliers (in comparison with others with the same 2theta
        value.) The process for this is that a median is computed for pixels
        within a small 2theta window and then the median difference is computed
        from magnitude of the difference for those pixels from that median. The
        medians are used for this rather than a standard deviation as the
        computation used here is less sensitive to outliers.
        (See :func:`GSASIIimage.AutoPixelMask` and
        :func:`scipy.stats.median_abs_deviation` for more details.)

        Mask is placed into the G2image object where it will be
        accessed during integration. Note that this increases the .gpx file
        size significantly; use :meth:`~G2Image.clearPixelMask` to delete
        this, if it need not be saved.

        This code is based on :func:`GSASIIimage.FastAutoPixelMask`
        but has been modified to recycle expensive computations
        where possible.

        :param float esdMul: Significance threshold applied to remove
          outliers. Default is 3. The larger this number, the fewer
          "glitches" that will be removed.
        :param float ttmin: A lower 2theta limit to be used for pixel
          searching. Pixels outside this region may be considered for
          establishing the medians, but only pixels with 2theta >= :attr:`ttmin`
          are masked. Default is 0.
        :param float ttmax: An upper 2theta limit to be used for pixel
          searching. Pixels outside this region may be considered for
          establishing the medians, but only pixels with 2theta < :attr:`ttmax`
          are masked. Default is 180.
        :param np.array FrameMask: An optional precomputed Frame mask
          (from :func:`~G2Image.MaskFrameMask`). Compute this once for
          a series of similar images to reduce computational time.
        :param np.array ThetaMap: An optional precomputed array that
          defines 2theta for each pixel, computed in
          :func:`~G2Image.MaskThetaMap`. Compute this once for
          a series of similar images to reduce computational time.
        :param bool fastmode: If True (default) fast Pixel map
          searching is done if the C module is available. If the
          module is not available or this is False, the pure Python
          implementatruion is used. It is not clear why False is
          ever needed.
        :param bool combineMasks: When True, the current Pixel mask
          will be combined with any previous Pixel map. If False (the
          default), the Pixel map from the current search will
          replace any previous ones. The reason for use of this as
          True would be where different :attr:`esdMul` values are
          used for different regions of the image (by setting
          :attr:`ttmin` & :attr:`ttmax`) so that the outlier level
          can be tuned by combining different searches.
        '''
        import math
        sind = lambda x: math.sin(x*math.pi/180.)
        if self.image is not None:
            Image = self.image
        else:
            Image = _getCorrImage(Readers['Image'],self.proj,self)
        Controls = self.getControls()
        Masks = self.getMasks()
        if FrameMask is None:
            frame = Masks['Frames']
            tam = ma.make_mask_none(Image.shape)
            if frame:
                tam = ma.mask_or(tam,ma.make_mask(np.abs(G2img.polymask(Controls,frame)-255)))
        else:
            tam = FrameMask
        if ThetaMap is None:
            TA = G2img.Make2ThetaAzimuthMap(Controls, (0, Image.shape[0]), (0, Image.shape[1]))[0]
        else:
            TA = ThetaMap
        LUtth = np.array(Controls['IOtth'])
        wave = Controls['wavelength']
        dsp0 = wave/(2.0*sind(LUtth[0]/2.0))
        dsp1 = wave/(2.0*sind(LUtth[1]/2.0))
        x0 = G2img.GetDetectorXY2(dsp0,0.0,Controls)[0]
        x1 = G2img.GetDetectorXY2(dsp1,0.0,Controls)[0]
        if not np.any(x0) or not np.any(x1):
            raise Exception
        numChans = int(1000*(x1-x0)/Controls['pixelSize'][0])//2
        if G2img.TestFastPixelMask() and fastmode:
            if GSASIIpath.binaryPath:
                import fmask
            else:
                from . import fmask
            G2fil.G2Print(f'Fast mask: Spots greater or less than {esdMul:.1f} of median abs deviation are masked')
            outMask = np.zeros_like(tam,dtype=bool).ravel()
            TThs = np.linspace(LUtth[0], LUtth[1], numChans, False)
            errmsg = ''
            try:
                fmask.mask(esdMul, tam.ravel(), TA.ravel(),
                                        Image.ravel(), TThs, outMask, ttmin, ttmax)
            except Exception as msg:
                errmsg = msg
            if errmsg:
                print('Exception in fmask.mask\n\t',errmsg)
                raise Exception(errmsg)
            outMask = outMask.reshape(Image.shape)
        else: # slow search, no sense using cache to save time
            Masks['SpotMask']['SearchMin'] = ttmin
            Masks['SpotMask']['SearchMax'] = ttmax
            outMask = G2img.AutoPixelMask(Image, Masks, Controls, numChans)
        if Masks['SpotMask'].get('spotMask') is not None and combineMasks:
            Masks['SpotMask']['spotMask'] |= outMask
        else:
            Masks['SpotMask']['spotMask'] = outMask

    def clearPixelMask(self):
        '''Removes a pixel map from an image, to reduce the .gpx file
        size & memory use
        '''
        self.getMasks()['SpotMask']['spotMask'] = None

class G2SmallAngle(G2ObjectWrapper):
    """Wrapper for SASD histograms (and hopefully, in the future, other
    small angle histogram types).

    Note that in a GSASIIscriptable script, instances of G2SmallAngle will be
    created by calls to
    :meth:`~G2Project.SAS`, :meth:`~G2Project.SASs`,
    or by :meth:`G2Project.Integrate`.
    Also, someday :meth:`G2Project.add_SAS`.
    Scripts should not try to create a :class:`G2SmallAngle` object directly.

    This object contains these class variables:
        * G2SmallAngle.proj: contains a reference to the :class:`G2Project`
          object that contains this histogram
        * G2SmallAngle.name: contains the name of the histogram
        * G2SmallAngle.data: contains the histogram's associated data
          in a dict with keys 'Comments', 'Limits', 'Instrument Parameters',
          'Substances', 'Sample Parameters' and 'Models'.
          Further documentation on SASD entries needs to be written.

    .. seealso::
        :meth:`~G2Project.add_SAS`
        :meth:`~G2Project.SAS`
        :meth:`~G2Project.SASs`
        :meth:`~G2Project.Integrate`
    """
    def __init__(self, data, proj, name):
        self.data = data
        self.name = name
        self.proj = proj

##########################
# Command Line Interface #
##########################
# Each of these takes an argparse.Namespace object as their argument,
# representing the parsed command-line arguments for the relevant subcommand.
# The argument specification for each is in the subcommands dictionary (see
# below)

commandhelp={}
commandhelp["create"] = "creates a GSAS-II project, optionally adding histograms and/or phases"
def create(args):
    """Implements the create command-line subcommand. This creates a GSAS-II project, optionally adding histograms and/or phases::

  usage: GSASIIscriptable.py create [-h] [-d HISTOGRAMS [HISTOGRAMS ...]]
                                  [-i IPARAMS [IPARAMS ...]]
                                  [-p PHASES [PHASES ...]]
                                  filename

positional arguments::

  filename              the project file to create. should end in .gpx

optional arguments::

  -h, --help            show this help message and exit
  -d HISTOGRAMS [HISTOGRAMS ...], --histograms HISTOGRAMS [HISTOGRAMS ...]
                        list of datafiles to add as histograms
  -i IPARAMS [IPARAMS ...], --iparams IPARAMS [IPARAMS ...]
                        instrument parameter file, must be one for every
                        histogram
  -p PHASES [PHASES ...], --phases PHASES [PHASES ...]
                        list of phases to add. phases are automatically
                        associated with all histograms given.

    """
    proj = G2Project(gpxfile=args.filename)

    hist_objs = []
    if args.histograms:
        for h,i in zip(args.histograms,args.iparams):
            G2fil.G2Print("Adding histogram from",h,"with instparm ",i)
            hist_objs.append(proj.add_powder_histogram(h, i))

    if args.phases:
        for p in args.phases:
            G2fil.G2Print("Adding phase from",p)
            proj.add_phase(p, histograms=hist_objs)
        G2fil.G2Print('Linking phase(s) to histogram(s):')
        for h in hist_objs:
            G2fil.G2Print ('   '+h.name)

    proj.save()

commandhelp["add"] = "adds histograms and/or phases to GSAS-II project"
def add(args):
    """Implements the add command-line subcommand. This adds histograms and/or phases to GSAS-II project::

  usage: GSASIIscriptable.py add [-h] [-d HISTOGRAMS [HISTOGRAMS ...]]
                               [-i IPARAMS [IPARAMS ...]]
                               [-hf HISTOGRAMFORMAT] [-p PHASES [PHASES ...]]
                               [-pf PHASEFORMAT] [-l HISTLIST [HISTLIST ...]]
                               filename


positional arguments::

  filename              the project file to open. Should end in .gpx

optional arguments::

  -h, --help            show this help message and exit
  -d HISTOGRAMS [HISTOGRAMS ...], --histograms HISTOGRAMS [HISTOGRAMS ...]
                        list of datafiles to add as histograms
  -i IPARAMS [IPARAMS ...], --iparams IPARAMS [IPARAMS ...]
                        instrument parameter file, must be one for every
                        histogram
  -hf HISTOGRAMFORMAT, --histogramformat HISTOGRAMFORMAT
                        format hint for histogram import. Applies to all
                        histograms
  -p PHASES [PHASES ...], --phases PHASES [PHASES ...]
                        list of phases to add. phases are automatically
                        associated with all histograms given.
  -pf PHASEFORMAT, --phaseformat PHASEFORMAT
                        format hint for phase import. Applies to all phases.
                        Example: -pf CIF
  -l HISTLIST [HISTLIST ...], --histlist HISTLIST [HISTLIST ...]
                        list of histgram indices to associate with added
                        phases. If not specified, phases are associated with
                        all previously loaded histograms. Example: -l 2 3 4

    """
    proj = G2Project(args.filename)

    if args.histograms:
        for h,i in zip(args.histograms,args.iparams):
            G2fil.G2Print("Adding histogram from",h,"with instparm ",i)
            proj.add_powder_histogram(h, i, fmthint=args.histogramformat)

    if args.phases:
        if not args.histlist:
            histlist = proj.histograms()
        else:
            histlist = [proj.histogram(i) for i in args.histlist]

        for p in args.phases:
            G2fil.G2Print("Adding phase from",p)
            proj.add_phase(p, histograms=histlist, fmthint=args.phaseformat)

        if not args.histlist:
            G2fil.G2Print('Linking phase(s) to all histogram(s)')
        else:
            G2fil.G2Print('Linking phase(s) to histogram(s):')
            for h in histlist:
                G2fil.G2Print('   '+h.name)
    proj.save()


commandhelp["dump"] = "Shows the contents of a GSAS-II project"
def dump(args):
    """Implements the dump command-line subcommand, which shows the contents of a GSAS-II project::

       usage: GSASIIscriptable.py dump [-h] [-d] [-p] [-r] files [files ...]

positional arguments::

  files

optional arguments::

  -h, --help        show this help message and exit
  -d, --histograms  list histograms in files, overrides --raw
  -p, --phases      list phases in files, overrides --raw
  -r, --raw         dump raw file contents, default

    """
    if not args.histograms and not args.phases:
        args.raw = True
    if args.raw:
        import IPython.lib.pretty as pretty

    for fname in args.files:
        if args.raw:
            proj, nameList = LoadDictFromProjFile(fname)
            print("file:", fname)
            print("names:", nameList)
            for key, val in proj.items():
                print(key, ":")
                pretty.pprint(val)
        else:
            proj = G2Project(fname)
            if args.histograms:
                hists = proj.histograms()
                for h in hists:
                    print(fname, "hist", h.id, h.name)
            if args.phases:
                phase_list = proj.phases()
                for p in phase_list:
                    print(fname, "phase", p.id, p.name)


commandhelp["browse"] = "Load a GSAS-II project and then open a IPython shell to browse it"
def IPyBrowse(args):
    """Load a .gpx file and then open a IPython shell to browse it::

  usage: GSASIIscriptable.py browse [-h] files [files ...]

positional arguments::

  files       list of files to browse

optional arguments::

  -h, --help  show this help message and exit

    """
    for fname in args.files:
        proj, nameList = LoadDictFromProjFile(fname)
        msg = "\nfname {} loaded into proj (dict) with names in nameList".format(fname)
        GSASIIpath.IPyBreak_base(msg)
        break


commandhelp["refine"] = '''
Conducts refinements on GSAS-II projects according to a list of refinement
steps in a JSON dict
'''
def refine(args):
    """Implements the refine command-line subcommand:
    conducts refinements on GSAS-II projects according to a JSON refinement dict::

        usage: GSASIIscriptable.py refine [-h] gpxfile [refinements]

positional arguments::

  gpxfile      the project file to refine
  refinements  json file of refinements to apply. if not present refines file
               as-is

optional arguments::

  -h, --help   show this help message and exit

    """
    proj = G2Project(args.gpxfile)
    if args.refinements is None:
        proj.refine()
    else:
        import json
        with open(args.refinements) as refs:
            refs = json.load(refs)
        if type(refs) is not dict:
            raise G2ScriptException("Error: JSON object must be a dict.")
        if "code" in refs:
            print("executing code:\n|  ",'\n|  '.join(refs['code']))
            exec('\n'.join(refs['code']))
        proj.do_refinements(refs['refinements'])


commandhelp["export"] = "Export phase as CIF"
def export(args):
    """Implements the export command-line subcommand: Exports phase as CIF::

      usage: GSASIIscriptable.py export [-h] gpxfile phase exportfile

positional arguments::

  gpxfile     the project file from which to export
  phase       identifier of phase to export
  exportfile  the .cif file to export to

optional arguments::

  -h, --help  show this help message and exit

    """
    proj = G2Project(args.gpxfile)
    phase = proj.phase(args.phase)
    phase.export_CIF(args.exportfile)


def _args_kwargs(*args, **kwargs):
    return args, kwargs

# A dictionary of the name of each subcommand, and a tuple
# of its associated function and a list of its arguments
# The arguments are passed directly to the add_argument() method
# of an argparse.ArgumentParser

subcommands = {"create":
               (create, [_args_kwargs('filename',
                                      help='the project file to create. should end in .gpx'),

                         _args_kwargs('-d', '--histograms',
                                      nargs='+',
                                      help='list of datafiles to add as histograms'),

                         _args_kwargs('-i', '--iparams',
                                      nargs='+',
                                      help='instrument parameter file, must be one'
                                           ' for every histogram'
                                      ),

                         _args_kwargs('-p', '--phases',
                                      nargs='+',
                                      help='list of phases to add. phases are '
                                           'automatically associated with all '
                                           'histograms given.')]),
               "add": (add, [_args_kwargs('filename',
                                      help='the project file to open. Should end in .gpx'),

                         _args_kwargs('-d', '--histograms',
                                      nargs='+',
                                      help='list of datafiles to add as histograms'),

                         _args_kwargs('-i', '--iparams',
                                      nargs='+',
                                      help='instrument parameter file, must be one'
                                           ' for every histogram'
                                      ),

                         _args_kwargs('-hf', '--histogramformat',
                                      help='format hint for histogram import. Applies to all'
                                           ' histograms'
                                      ),

                         _args_kwargs('-p', '--phases',
                                      nargs='+',
                                      help='list of phases to add. phases are '
                                           'automatically associated with all '
                                           'histograms given.'),

                         _args_kwargs('-pf', '--phaseformat',
                                      help='format hint for phase import. Applies to all'
                                           ' phases. Example: -pf CIF'
                                      ),

                         _args_kwargs('-l', '--histlist',
                                      nargs='+',
                                      help='list of histgram indices to associate with added'
                                           ' phases. If not specified, phases are'
                                           ' associated with all previously loaded'
                                           ' histograms. Example: -l 2 3 4')]),

               "dump": (dump, [_args_kwargs('-d', '--histograms',
                                     action='store_true',
                                     help='list histograms in files, overrides --raw'),

                               _args_kwargs('-p', '--phases',
                                            action='store_true',
                                            help='list phases in files, overrides --raw'),

                               _args_kwargs('-r', '--raw',
                                      action='store_true', help='dump raw file contents, default'),

                               _args_kwargs('files', nargs='+')]),

               "refine":
               (refine, [_args_kwargs('gpxfile', help='the project file to refine'),
                         _args_kwargs('refinements',
                                      help='JSON file of refinements to apply. if not present'
                                           ' refines file as-is',
                                      default=None,
                                      nargs='?')]),

               "export": (export, [_args_kwargs('gpxfile',
                                                help='the project file from which to export'),
                                   _args_kwargs('phase', help='identifier of phase to export'),
                                   _args_kwargs('exportfile', help='the .cif file to export to')]),
               "browse": (IPyBrowse, [_args_kwargs('files', nargs='+',
                                                   help='list of files to browse')])}


def main():
    '''The command-line interface for calling GSASIIscriptable as a shell command,
    where it is expected to be called as::

       python GSASIIscriptable.py <subcommand> <file.gpx> <options>

    The following subcommands are defined:

        * create, see :func:`create`
        * add, see :func:`add`
        * dump, see :func:`dump`
        * refine, see :func:`refine`
        * export, :func:`export`
        * browse, see :func:`IPyBrowse`

    .. seealso::
        :func:`create`
        :func:`add`
        :func:`dump`
        :func:`refine`
        :func:`export`
        :func:`IPyBrowse`
    '''
    parser = argparse.ArgumentParser(description=
        "Use of "+os.path.split(__file__)[1]+" Allows GSAS-II actions from command line."
        )
    subs = parser.add_subparsers()

    # Create all of the specified subparsers
    for name, (func, args) in subcommands.items():
        new_parser = subs.add_parser(name,help=commandhelp.get(name),
                                     description='Command "'+name+'" '+commandhelp.get(name))
        for listargs, kwds in args:
            new_parser.add_argument(*listargs, **kwds)
        new_parser.set_defaults(func=func)

    # Parse and trigger subcommand
    result = parser.parse_args()
    try:
        result.func(result)
    except:
        print('Use {} -h or --help'.format(sys.argv[0]))

def dictDive(d,search="",keylist=[],firstcall=True,l=None):
    '''Recursive routine to scan a nested dict. Reports a list of keys
    and the associated type and value for that key.

    :param dict d: a dict that will be scanned
    :param str search: an optional search string. If non-blank,
       only entries where one of the keys constains search (case ignored)
    :param list keylist: a list of keys to apply to the dict.
    :param bool firstcall: do not specify
    :param list l: do not specify

    :returns: a list of keys located by this routine
       in form [([keylist], type, value),...] where if keylist is ['a','b','c']
       then d[['a']['b']['c'] will have the value.

    This routine can be called in a number of ways, as are shown in a few
    examples:

    >>> for i in G2sc.dictDive(p.data['General'],'paw'): print(i)
    ...
    (['Pawley dmin'], <class 'float'>, 1.0)
    (['doPawley'], <class 'bool'>, False)
    (['Pawley dmax'], <class 'float'>, 100.0)
    (['Pawley neg wt'], <class 'float'>, 0.0)
    >>>
    >>> for i in G2sc.dictDive(p.data,'paw',['General']): print(i)
    ...
    (['General', 'Pawley dmin'], <class 'float'>, 1.0)
    (['General', 'doPawley'], <class 'bool'>, False)
    (['General', 'Pawley dmax'], <class 'float'>, 100.0)
    (['General', 'Pawley neg wt'], <class 'float'>, 0.0)
    >>>
    >>> for i in G2sc.dictDive(p.data,'',['General','doPawley']): print(i)
    ...
    (['General', 'doPawley'], <class 'bool'>, False)

    '''
    if firstcall:
        for k in keylist:
            d = d[k]
        l = []
#        first = True
    if type(d) is dict:
        [dictDive(d[key],search,keylist+[key],False,l) for key in d]
    elif search:
        for key in keylist:
            if search.lower() in key.lower():
                l.append((keylist,type(d),d))
                break
        return
    else:
        l.append((keylist,type(d),d))
    if firstcall:
        return l

if __name__ == '__main__':
    #fname='/tmp/corundum-template.gpx'
    #prj = G2Project(fname)
    main()
