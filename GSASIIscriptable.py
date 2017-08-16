#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2017-04-12 15:12:45 -0500 (Wed, 12 Apr 2017) $
# $Author: vondreele $
# $Revision: 2777 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIscriptable.py $
# $Id: GSASIIIO.py 2777 2017-04-12 20:12:45Z vondreele $
########### SVN repository information ###################
"""
*GSASIIscriptable: Scripting Tools*
-----------------------------------

Routines for reading, writing, modifying and creating GSAS-II project (.gpx) files.

Supports a command line interface as well.

Look at :class:`G2Project` to start.

=====================
Refinement parameters
=====================

There are three classes of refinement parameters:

    * Histogram. Turned on and off through :func:`~G2PwdrData.set_refinements`
      and :func:`~G2PwdrData.clear_refinements`
    * Phase. Turned on and off through :func:`~G2Phase.set_refinements`
      and :func:`~G2Phase.clear_refinements`
    * Histogram-and-phase (HAP). Turned on and off through
      :func:`~G2Phase.set_HAP_refinements` and :func:`~G2Phase.clear_HAP_refinements`
"""
from __future__ import division, print_function # needed?
import os.path as ospath
import datetime as dt
import sys
import cPickle
import imp
import copy
import os
import random as ran

import numpy.ma as ma
import scipy.interpolate as si
import numpy as np
import scipy as sp

import GSASIIpath
GSASIIpath.SetBinaryPath(False) # would rather have this in __name__ == '__main__' stanza
import GSASIIobj as G2obj
import GSASIIpwd as G2pwd
import GSASIIstrMain as G2strMain
import GSASIIspc as G2spc
import GSASIIElem as G2elem

# Delay imports to not slow down small scripts
G2fil = None
PwdrDataReaders = []
PhaseReaders = []


def LoadG2fil():
    """Delay importing this module, it is slow"""
    global G2fil
    global PwdrDataReaders
    global PhaseReaders
    if G2fil is None:
        import GSASIIfiles
        G2fil = GSASIIfiles
        PwdrDataReaders = G2fil.LoadImportRoutines("pwd", "Powder_Data")
        PhaseReaders = G2fil.LoadImportRoutines("phase", "Phase")


def LoadDictFromProjFile(ProjFile):
    '''Read a GSAS-II project file and load items to dictionary
    :param str ProjFile: GSAS-II project (name.gpx) full file name
    :returns: Project,nameList, where

      * Project (dict) is a representation of gpx file following the GSAS-II tree struture
        for each item: key = tree name (e.g. 'Controls','Restraints',etc.), data is dict
        data dict = {'data':item data whch may be list, dict or None,'subitems':subdata (if any)}
      * nameList (list) has names of main tree entries & subentries used to reconstruct project file

    Example for fap.gpx::

      Project = {                 #NB:dict order is not tree order
        u'Phases':{'data':None,'fap':{phase dict}},
        u'PWDR FAP.XRA Bank 1':{'data':[histogram data list],'Comments':comments,'Limits':limits, etc},
        u'Rigid bodies':{'data': {rigid body dict}},
        u'Covariance':{'data':{covariance data dict}},
        u'Controls':{'data':{controls data dict}},
        u'Notebook':{'data':[notebook list]},
        u'Restraints':{'data':{restraint data dict}},
        u'Constraints':{'data':{constraint data dict}}]
        }
      nameList = [                #NB: reproduces tree order
        [u'Notebook',],
        [u'Controls',],
        [u'Covariance',],
        [u'Constraints',],
        [u'Restraints',],
        [u'Rigid bodies',],
        [u'PWDR FAP.XRA Bank 1',
             u'Comments',
             u'Limits',
             u'Background',
             u'Instrument Parameters',
             u'Sample Parameters',
             u'Peak List',
             u'Index Peak List',
             u'Unit Cells List',
             u'Reflection Lists'],
        [u'Phases', u'fap']
        ]
    '''
    # Let IOError be thrown if file does not exist
    # if not ospath.exists(ProjFile):
    #     print ('\n*** Error attempt to open project file that does not exist:\n   '+
    #         str(ProjFile))
    #     return
    file = open(ProjFile,'rb')
    # print('loading from file: {}'.format(ProjFile))
    Project = {}
    nameList = []
    try:
        while True:
            try:
                data = cPickle.load(file)
            except EOFError:
                break
            datum = data[0]
            Project[datum[0]] = {'data':datum[1]}
            nameList.append([datum[0],])
            for datus in data[1:]:
                Project[datum[0]][datus[0]] = datus[1]
                nameList[-1].append(datus[0])
        file.close()
        # print('project load successful')
    except:
        raise IOError("Error reading file "+str(ProjFile)+". This is not a GSAS-II .gpx file")
    finally:
        file.close()
    return Project,nameList

def SaveDictToProjFile(Project,nameList,ProjFile):
    '''Save a GSAS-II project file from dictionary/nameList created by LoadDictFromProjFile

    :param dict Project: representation of gpx file following the GSAS-II
        tree structure as described for LoadDictFromProjFile
    :param list nameList: names of main tree entries & subentries used to reconstruct project file
    :param str ProjFile: full file name for output project.gpx file (including extension)
    '''
    file = open(ProjFile,'wb')
    # print('save to file: {}'.format(ProjFile))
    try:
        for name in nameList:
            data = []
            item = Project[name[0]]
            data.append([name[0],item['data']])
            for item2 in name[1:]:
                data.append([item2,item[item2]])
            cPickle.dump(data,file,1)
    finally:
        file.close()
    # print('project save successful')

def ImportPowder(reader,filename):
    '''Use a reader to import a powder diffraction data file

    :param str reader: a scriptable reader
    :param str filename: full name of powder data file; can be "multi-Bank" data

    :returns list rdlist: list of reader objects containing powder data, one for each
        "Bank" of data encountered in file. Items in reader object of interest are:

          * rd.comments: list of str: comments found on powder file
          * rd.dnames: list of str: data nammes suitable for use in GSASII data tree NB: duplicated in all rd entries in rdlist
          * rd.powderdata: list of numpy arrays: pos,int,wt,zeros,zeros,zeros as needed for a PWDR entry in  GSASII data tree.
    '''
    rdfile,rdpath,descr = imp.find_module(reader)
    rdclass = imp.load_module(reader,rdfile,rdpath,descr)
    rd = rdclass.GSAS_ReaderClass()
    if not rd.scriptable:
        print(u'**** ERROR: '+reader+u' is not a scriptable reader')
        return None
    fl = open(filename,'rb')
    rdlist = []
    if rd.ContentsValidator(fl):
        fl.seek(0)
        repeat = True
        rdbuffer = {} # create temporary storage for file reader
        block = 0
        while repeat: # loop if the reader asks for another pass on the file
            block += 1
            repeat = False
            rd.objname = ospath.basename(filename)
            flag = rd.Reader(filename,fl,None,buffer=rdbuffer,blocknum=block,)
            if flag:
                rdlist.append(copy.deepcopy(rd)) # save the result before it is written over
                if rd.repeat:
                    repeat = True
        return rdlist
    print(rd.errors)
    return None

def SetDefaultDData(dType,histoName,NShkl=0,NDij=0):
    '''Create an initial Histogram dictionary

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    '''
    if dType in ['SXC','SNC']:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {'Tbar':0.1,'Cos2TM':0.955,
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}],
            'Flack':[0.0,False]}
    elif dType == 'SNT':
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}]}
    elif 'P' in dType:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,False],
            'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},[],0.1],
            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                [1.,1.,1.,0.,0.,0.],6*[False,]],
            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                NShkl*[0.01,],NShkl*[False,]],
            'HStrain':[NDij*[0.0,],NDij*[False,]],
            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}


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
    if 'modulated' in data['General']['Type']:
        data['General']['Modulated'] = True
        data['General']['Type'] = 'nuclear'


def SetupGeneral(data, dirname):
    """Helps initialize phase data.

    From GSASIIphsGui.py, function of the same name. Minor changes for imports etc.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    mapDefault = {'MapType':'','RefList':'','Resolution':0.5,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}
    generalData = data['General']
    atomData = data['Atoms']
    generalData['AtomTypes'] = []
    generalData['Isotopes'] = {}

    if 'Isotope' not in generalData:
        generalData['Isotope'] = {}
    if 'Data plot type' not in generalData:
        generalData['Data plot type'] = 'Mustrain'
    if 'POhkl' not in generalData:
        generalData['POhkl'] = [0,0,1]
    if 'Map' not in generalData:
        generalData['Map'] = mapDefault.copy()
    if 'Flip' not in generalData:
        generalData['Flip'] = {'RefList':'','Resolution':0.5,'Norm element':'None',
            'k-factor':0.1,'k-Max':20.,}
    if 'testHKL' not in generalData['Flip']:
        generalData['Flip']['testHKL'] = [[0,0,2],[2,0,0],[1,1,1],[0,2,0],[1,2,3]]
    if 'doPawley' not in generalData:
        generalData['doPawley'] = False     #ToDo: change to ''
    if 'Pawley dmin' not in generalData:
        generalData['Pawley dmin'] = 1.0
    if 'Pawley neg wt' not in generalData:
        generalData['Pawley neg wt'] = 0.0
    if 'Algolrithm' in generalData.get('MCSA controls',{}) or \
        'MCSA controls' not in generalData:
        generalData['MCSA controls'] = {'Data source':'','Annealing':[50.,0.001,50],
        'dmin':2.0,'Algorithm':'log','Jump coeff':[0.95,0.5],'boltzmann':1.0,
        'fast parms':[1.0,1.0,1.0],'log slope':0.9,'Cycles':1,'Results':[],'newDmin':True}
    if 'AtomPtrs' not in generalData:
        generalData['AtomPtrs'] = [3,1,7,9]
        if generalData['Type'] == 'macromolecular':
            generalData['AtomPtrs'] = [6,4,10,12]
        elif generalData['Type'] == 'magnetic':
            generalData['AtomPtrs'] = [3,1,10,12]
    if generalData['Type'] in ['modulated',]:
        generalData['Modulated'] = True
        generalData['Type'] = 'nuclear'
        if 'Super' not in generalData:
            generalData['Super'] = 1
            generalData['SuperVec'] = [[0,0,.1],False,4]
            generalData['SSGData'] = {}
        if '4DmapData' not in generalData:
            generalData['4DmapData'] = mapDefault.copy()
            generalData['4DmapData'].update({'MapType':'Fobs'})
    if 'Modulated' not in generalData:
        generalData['Modulated'] = False
    if 'HydIds' not in generalData:
        generalData['HydIds'] = {}
    cx,ct,cs,cia = generalData['AtomPtrs']
    generalData['NoAtoms'] = {}
    generalData['BondRadii'] = []
    generalData['AngleRadii'] = []
    generalData['vdWRadii'] = []
    generalData['AtomMass'] = []
    generalData['Color'] = []
    if generalData['Type'] == 'magnetic':
        generalData['MagDmin'] = generalData.get('MagDmin',1.0)
        landeg = generalData.get('Lande g',[])
    generalData['Mydir'] = dirname
    badList = {}
    for iat,atom in enumerate(atomData):
        atom[ct] = atom[ct].lower().capitalize()              #force to standard form
        if generalData['AtomTypes'].count(atom[ct]):
            generalData['NoAtoms'][atom[ct]] += atom[cx+3]*float(atom[cs+1])
        elif atom[ct] != 'UNK':
            Info = G2elem.GetAtomInfo(atom[ct])
            if not Info:
                if atom[ct] not in badList:
                    badList[atom[ct]] = 0
                badList[atom[ct]] += 1
                atom[ct] = 'UNK'
                continue
            atom[ct] = Info['Symbol'] # N.B. symbol might be changed by GetAtomInfo
            generalData['AtomTypes'].append(atom[ct])
            generalData['Z'] = Info['Z']
            generalData['Isotopes'][atom[ct]] = Info['Isotopes']
            generalData['BondRadii'].append(Info['Drad'])
            generalData['AngleRadii'].append(Info['Arad'])
            generalData['vdWRadii'].append(Info['Vdrad'])
            if atom[ct] in generalData['Isotope']:
                if generalData['Isotope'][atom[ct]] not in generalData['Isotopes'][atom[ct]]:
                    isotope = generalData['Isotopes'][atom[ct]].keys()[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]]['Mass'])
            else:
                generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                if 'Nat. Abund.' not in generalData['Isotopes'][atom[ct]]:
                    isotope = generalData['Isotopes'][atom[ct]].keys()[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Mass'])
            generalData['NoAtoms'][atom[ct]] = atom[cx+3]*float(atom[cs+1])
            generalData['Color'].append(Info['Color'])
            if generalData['Type'] == 'magnetic':
                if len(landeg) < len(generalData['AtomTypes']):
                    landeg.append(2.0)
    if generalData['Type'] == 'magnetic':
        generalData['Lande g'] = landeg[:len(generalData['AtomTypes'])]

    if badList:
        msg = 'Warning: element symbol(s) not found:'
        for key in badList:
            msg += '\n\t' + key
            if badList[key] > 1:
                msg += ' (' + str(badList[key]) + ' times)'
        raise Exception("Phase error:\n" + msg)
        # wx.MessageBox(msg,caption='Element symbol error')
    F000X = 0.
    F000N = 0.
    for i,elem in enumerate(generalData['AtomTypes']):
        F000X += generalData['NoAtoms'][elem]*generalData['Z']
        isotope = generalData['Isotope'][elem]
        F000N += generalData['NoAtoms'][elem]*generalData['Isotopes'][elem][isotope]['SL'][0]
    generalData['F000X'] = F000X
    generalData['F000N'] = F000N
    import GSASIImath as G2mth
    generalData['Mass'] = G2mth.getMass(generalData)


############################################
############ CIF export helpers ############
############################################
## These functions are translated from
## exports/G2export_CIF.py
## They are used in G2Phase.export_CIF



def make_empty_project(author=None, filename=None):
    """Creates an dictionary in the style of GSASIIscriptable, for an empty
    project.

    If no author name or filename is supplied, 'no name' and
    <current dir>/test_output.gpx are used , respectively.

    Returns: project dictionary, name list

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    if not filename:
        filename = os.path.join(os.getcwd(), 'test_output.gpx')
    else:
        filename = os.path.abspath(filename)
    gsasii_version = str(GSASIIpath.GetVersionNumber())
    LoadG2fil()
    import matplotlib as mpl
    python_library_versions = G2fil.get_python_versions([mpl, np, sp])

    controls_data = dict(G2obj.DefaultControls)
    controls_data['LastSavedAs'] = unicode(filename)
    controls_data['LastSavedUsing'] = gsasii_version
    controls_data['PythonVersions'] = python_library_versions
    if author:
        controls_data['Author'] = author

    output = {'Constraints': {'data': {'HAP': [], 'Hist': [], 'Phase': [],
                                       'Global': []}},
              'Controls': {'data': controls_data},
              u'Covariance': {'data': {}},
              u'Notebook': {'data': ['']},
              u'Restraints': {'data': {}},
              u'Rigid bodies': {'data': {'RBIds': {'Residue': [], 'Vector': []},
                                'Residue': {'AtInfo': {}},
                                'Vector':  {'AtInfo': {}}}}}

    names = [[u'Notebook'], [u'Controls'], [u'Covariance'],
             [u'Constraints'], [u'Restraints'], [u'Rigid bodies']]

    return output, names


class G2ImportException(Exception):
    pass


def import_generic(filename, readerlist):
    """Attempt to import a filename, using a list of reader objects.

    Returns the first reader object which worked."""
    # Translated from OnImportGeneric method in GSASII.py
    primaryReaders, secondaryReaders = [], []
    for reader in readerlist:
        flag = reader.ExtensionValidator(filename)
        if flag is None:
            secondaryReaders.append(reader)
        elif flag:
            primaryReaders.append(reader)
    if not secondaryReaders and not primaryReaders:
        raise G2ImportException("Could not read file: ", filename)

    with open(filename, 'Ur') as fp:
        rd_list = []

        for rd in primaryReaders + secondaryReaders:
            # Initialize reader
            rd.selections = []
            rd.dnames = []
            rd.ReInitialize()
            # Rewind file
            fp.seek(0)
            rd.errors = ""
            if not rd.ContentsValidator(fp):
                # Report error
                pass
            if len(rd.selections) > 1:
                # Select data?
                # GSASII.py:543
                raise G2ImportException("Not sure what data to select")

            block = 0
            rdbuffer = {}
            fp.seek(0)
            repeat = True
            while repeat:
                repeat = False
                block += 1
                rd.objname = os.path.basename(filename)
                flag = rd.Reader(filename, fp, buffer=rdbuffer, blocknum=block)
                if flag:
                    # Omitting image loading special cases
                    rd.readfilename = filename
                    rd_list.append(copy.deepcopy(rd))
                    repeat = rd.repeat
                else:
                    raise G2ImportException("{}.Reader() returned:".format(rd),
                                            flag)

            if rd_list:
                if rd.warnings:
                    print("Read warning by", rd.formatName, "reader:",
                          rd.warnings, file=sys.stderr)
                return rd_list
    raise G2ImportException("No reader could read file: " + filename)


def load_iprms(instfile, reader):
    """Loads instrument parameters from a file, and edits the
    given reader.

    Returns a 2-tuple of (Iparm1, Iparm2) parameters
    """
    LoadG2fil()
    ext = os.path.splitext(instfile)[1]

    if ext.lower() == '.instprm':
        # New GSAS File, load appropriately
        with open(instfile) as f:
            lines = f.readlines()
        bank = reader.powderentry[2]
        numbanks = reader.numbanks
        iparms = G2fil.ReadPowderInstprm(lines, bank, numbanks, reader)

        reader.instfile = instfile
        reader.instmsg = 'GSAS-II file' + instfile
        return iparms
    elif ext.lower() not in ('.prm', '.inst', '.ins'):
        raise ValueError('Expected .prm file, found: ', instfile)

    # It's an old GSAS file, load appropriately
    Iparm = {}
    with open(instfile, 'Ur') as fp:
        for line in fp:
            if '#' in line:
                continue
            Iparm[line[:12]] = line[12:-1]
    ibanks = int(Iparm.get('INS   BANK  ', '1').strip())
    if ibanks == 1:
        reader.instbank = 1
        reader.powderentry[2] = 1
        reader.instfile = instfile
        reader.instmsg = instfile + ' bank ' + str(reader.instbank)
        return G2fil.SetPowderInstParms(Iparm, reader)
    # TODO handle >1 banks
    raise NotImplementedError("Check GSASIIfiles.py: ReadPowderInstprm")

def load_pwd_from_reader(reader, instprm, existingnames=[]):
    """Loads powder data from a reader object, and assembles it into a GSASII data tree.

    :returns: (name, tree) - 2-tuple of the histogram name (str), and data

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    HistName = 'PWDR ' + G2obj.StripUnicode(reader.idstring, '_')
    HistName = unicode(G2obj.MakeUniqueLabel(HistName, existingnames))

    try:
        Iparm1, Iparm2 = instprm
    except ValueError:
        Iparm1, Iparm2 = load_iprms(instprm, reader)

    Ymin = np.min(reader.powderdata[1])
    Ymax = np.max(reader.powderdata[1])
    valuesdict = {'wtFactor': 1.0,
                  'Dummy': False,
                  'ranId': ran.randint(0, sys.maxint),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'qPlot': False, 'dPlot': False, 'sqrtPlot': False,
                  'Yminmax': [Ymin, Ymax]}
    reader.Sample['ranId'] = valuesdict['ranId']

    # Ending keys:
    # [u'Reflection Lists',
    #  u'Limits',
    #  'data',
    #  u'Index Peak List',
    #  u'Comments',
    #  u'Unit Cells List',
    #  u'Sample Parameters',
    #  u'Peak List',
    #  u'Background',
    #  u'Instrument Parameters']
    Tmin = np.min(reader.powderdata[0])
    Tmax = np.max(reader.powderdata[0])

    default_background = [['chebyschev', False, 3, 1.0, 0.0, 0.0],
                          {'nDebye': 0, 'debyeTerms': [], 'nPeaks': 0, 'peaksList': []}]

    output_dict = {u'Reflection Lists': {},
                   u'Limits': reader.pwdparms.get('Limits', [(Tmin, Tmax), [Tmin, Tmax]]),
                   u'data': [valuesdict, reader.powderdata, HistName],
                   u'Index Peak List': [[], []],
                   u'Comments': reader.comments,
                   u'Unit Cells List': [],
                   u'Sample Parameters': reader.Sample,
                   u'Peak List': {'peaks': [], 'sigDict': {}},
                   u'Background': reader.pwdparms.get('Background', default_background),
                   u'Instrument Parameters': [Iparm1, Iparm2],
                   }

    names = [u'Comments',
             u'Limits',
             u'Background',
             u'Instrument Parameters',
             u'Sample Parameters',
             u'Peak List',
             u'Index Peak List',
             u'Unit Cells List',
             u'Reflection Lists']

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
            for i in xrange(len(from_)):
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
    """Represents an entire GSAS-II project."""
    def __init__(self, gpxfile=None, author=None, filename=None):
        """Loads a GSAS-II project from a specified filename.

        :param str gpxfile: Existing .gpx file to be loaded. If nonexistent,
            creates an empty project.
        :param str author: Author's name. Optional.
        :param str filename: The filename the project should be saved to in
            the future. If both filename and gpxfile are present, the project is
            loaded from the gpxfile, then set to  save to filename in the future"""
        if gpxfile is None:
            filename = os.path.abspath(os.path.expanduser(filename))
            self.filename = filename
            self.data, self.names = make_empty_project(author=author, filename=filename)
        elif isinstance(gpxfile, str):
            # TODO set author, filename
            self.filename = os.path.abspath(os.path.expanduser(gpxfile))
            self.data, self.names = LoadDictFromProjFile(gpxfile)
            self.index_ids()
        else:
            raise ValueError("Not sure what to do with gpxfile")

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
        # TODO update LastSavedUsing ?
        if filename:
            filename = os.path.abspath(os.path.expanduser(filename))
            self.data['Controls']['data']['LastSavedAs'] = filename
            self.filename = filename
        elif not self.filename:
            raise AttributeError("No file name to save to")
        SaveDictToProjFile(self.data, self.names, self.filename)

    def add_powder_histogram(self, datafile, iparams):
        """Loads a powder data histogram into the project.

        Automatically checks for an instrument parameter file, or one can be
        provided.

        :param str datafile: The powder data file to read, a filename.
        :param str iparams: The instrument parameters file, a filename.

        :returns: A :class:`G2PwdrData` object representing
            the histogram
        """
        LoadG2fil()
        datafile = os.path.abspath(os.path.expanduser(datafile))
        iparams = os.path.abspath(os.path.expanduser(iparams))
        pwdrreaders = import_generic(datafile, PwdrDataReaders)
        histname, new_names, pwdrdata = load_pwd_from_reader(
                                          pwdrreaders[0], iparams,
                                          [h.name for h in self.histograms()])
        if histname in self.data:
            print("Warning - redefining histogram", histname)
        else:
            if self.names[-1][0] == 'Phases':
                self.names.insert(-1, new_names)
            else:
                self.names.append(new_names)
        self.data[histname] = pwdrdata
        return self.histogram(histname)

    def add_phase(self, phasefile, phasename=None, histograms=[]):
        """Loads a phase into the project from a .cif file

        :param str phasefile: The CIF file from which to import the phase.
        :param str phasename: The name of the new phase, or None for the default
        :param list histograms: The names of the histograms to associate with
            this phase

        :returns: A :class:`G2Phase` object representing the
            new phase.
        """
        LoadG2fil()
        phasefile = os.path.abspath(os.path.expanduser(phasefile))

        # TODO handle multiple phases in a file
        phasereaders = import_generic(phasefile, PhaseReaders)
        phasereader = phasereaders[0]

        phasename = phasename or phasereader.Phase['General']['Name']
        phaseNameList = [p.name for p in self.phases()]
        phasename = G2obj.MakeUniqueLabel(phasename, phaseNameList)
        phasereader.Phase['General']['Name'] = phasename

        if 'Phases' not in self.data:
            self.data[u'Phases'] = { 'data': None }
        assert phasename not in self.data['Phases'], "phase names should be unique"
        self.data['Phases'][phasename] = phasereader.Phase

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
        generalData = data['General']
        SGData = generalData['SGData']
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        Super = generalData.get('Super', 0)
        if Super:
            SuperVec = np.array(generalData['SuperVec'][0])
        else:
            SuperVec = []
        UseList = data['Histograms']

        for histname in histograms:
            if histname.startswith('HKLF '):
                raise NotImplementedError("Does not support HKLF yet")
            elif histname.startswith('PWDR '):
                hist = self.histogram(histname)
                hist['Reflection Lists'][generalData['Name']] = {}
                UseList[histname] = SetDefaultDData('PWDR', histname, NShkl=NShkl, NDij=NDij)
                for key, val in [('Use', True), ('LeBail', False),
                                 ('newLeBail', True),
                                 ('Babinet', {'BabA': [0.0, False],
                                              'BabU': [0.0, False]})]:
                    if key not in UseList[histname]:
                        UseList[histname][key] = val
            else:
                raise NotImplementedError("Unexpected histogram" + histname)

        for obj in self.names:
            if obj[0] == 'Phases':
                phasenames = obj
                break
        else:
            phasenames = [u'Phases']
            self.names.append(phasenames)
        phasenames.append(unicode(phasename))

        # TODO should it be self.filename, not phasefile?
        SetupGeneral(data, os.path.dirname(phasefile))
        self.index_ids()

        return self.phase(phasename)

    def reload(self):
        """Reload self from self.filename"""
        data, names = LoadDictFromProjFile(self.filename)
        self.names = names
        # Need to deep copy the new data file data into the current tree,
        # so that any existing G2Phase, or G2PwdrData objects will still be
        # valid
        _deep_copy_into(from_=data, into=self.data)

    def refine(self, newfile=None, printFile=None, makeBack=False):
        # index_ids will automatically save the project
        self.index_ids()
        # TODO G2strMain does not properly use printFile
        G2strMain.Refine(self.filename, makeBack=makeBack)
        # Reload yourself
        self.reload()

    def histogram(self, histname):
        """Returns the histogram named histname, or None if it does not exist.

        :param histname: The name of the histogram (str), or ranId or index.
        :returns: A :class:`G2PwdrData` object, or None if
            the histogram does not exist

        .. seealso::
            :func:`~GSASIIscriptable.G2Project.histograms`
            :func:`~GSASIIscriptable.G2Project.phase`
            :func:`~GSASIIscriptable.G2Project.phases`
            """
        if histname in self.data:
            return G2PwdrData(self.data[histname], self)
        for key, val in G2obj.HistIdLookup.items():
            name, ranId = val
            # histname can be either ranId (key) or index (val)
            if ranId == histname or key == str(histname):
                return self.histogram(name)

    def histograms(self):
        """Return a list of all histograms, as
        :class:`G2PwdrData` objects

        .. seealso::
            :func:`~GSASIIscriptable.G2Project.histograms`
            :func:`~GSASIIscriptable.G2Project.phase`
            :func:`~GSASIIscriptable.G2Project.phases`
            """
        output = []
        for obj in self.names:
            if len(obj) > 1 and obj[0] != u'Phases':
                output.append(self.histogram(obj[0]))
        return output

    def phase(self, phasename):
        """
        Gives an object representing the specified phase in this project.

        :param str phasename: The name of the desired phase. Either the name
            (str), the phase's ranId, or the phase's index
        :returns: A :class:`G2Phase` object
        :raises: KeyError

        .. seealso::
            :func:`~GSASIIscriptable.G2Project.histograms`
            :func:`~GSASIIscriptable.G2Project.phase`
            :func:`~GSASIIscriptable.G2Project.phases`
            """
        phases = self.data['Phases']
        if phasename in phases:
            return G2Phase(phases[phasename], phasename, self)
        for key, val in G2obj.PhaseIdLookup.items():
            name, ranId = val
            # phasename can be either ranId (key) or index (val)
            if ranId == phasename or key == str(phasename):
                return self.phase(name)

    def phases(self):
        """
        Returns a list of all the phases in the project.

        :returns: A list of :class:`G2Phase` objects

        .. seealso::
            :func:`~GSASIIscriptable.G2Project.histogram`
            :func:`~GSASIIscriptable.G2Project.histograms`
            :func:`~GSASIIscriptable.G2Project.phase`
            """
        for obj in self.names:
            if obj[0] == 'Phases':
                return [self.phase(p) for p in obj[1:]]
        return []

    def do_refinements(self, refinements, histogram='all', phase='all',
                       outputnames=None):
        """Conducts a series of refinements.

        :param list refinements: A list of dictionaries defining refinements
        :param str histogram: Name of histogram for refinements to be applied
            to, or 'all'
        :param str phase: Name of phase for refinements to be applied to, or
            'all'
        """
        if outputnames:
            if len(refinements) != len(outputnames):
                raise ValueError("Should have same number of outuputs to"
                                 "refinements")
        else:
            outputnames = [None for r in refinements]

        for output, refinement in zip(outputnames, refinements):
            self.set_refinement(refinement, histogram)
            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinement:
                temp = {'set': refinement['once']}
                self.set_refinement(temp, histogram, phase)

            if output:
                self.save(output)

            self.refine()  # newFile=output)

            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinement:
                temp = {'clear': refinement['once']}
                self.set_refinement(temp, histogram, phase)

    def set_refinement(self, refinement, histogram='all', phase='all'):
        """Apply specified refinements to a given histogram(s) or phase(s).

        Refinement parameters are categorize in three groups:

        1. Histogram parameters
        2. Phase parameters
        3. Histogram-and-Phase (HAP) parameters

        :param dict refinement: The refinements to be conducted
        :param histogram: Either a name of a histogram (str), a list of
            histogram names, or 'all' (default)
        :param phase: Either a name of a phase (str), a list of phase names, or
            'all' (default)

        .. seealso::
            :func:`~G2PwdrData.set_refinements`
            :func:`~G2PwdrData.clear_refinements`
            :func:`~G2Phase.set_refinements`
            :func:`~G2Phase.clear_refinements`
            :func:`~G2Phase.set_HAP_refinements`
            :func:`~G2Phase.clear_HAP_refinements`"""

        if histogram == 'all':
            hists = self.histograms()
        elif isinstance(histogram, str) or isinstance(histogram, int):
            hists = [self.histogram(histogram)]
        else:
            hists = [self.histogram(name) for name in histogram]

        if phase == 'all':
            phases = self.phases()
        elif isinstance(phase, str) or isinstance(phase, int):
            phases = [self.phase(phase)]
        else:
            phases = [self.phase(name) for name in phase]


        # TODO: HAP parameters:
        #   Babinet
        #   Extinction
        #   HStrain
        #   Mustrain
        #   Pref. Ori
        #   Size

        pwdr_set = {}
        phase_set = {}
        hap_set = {}
        for key, val in refinement.get('set', {}).items():
            # Apply refinement options
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_set[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_set[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                hap_set[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            hist.set_refinements(pwdr_set)
        for phase in phases:
            phase.set_refinements(phase_set)
        for phase in phases:
            phase.set_HAP_refinements(hap_set, hists)

        pwdr_clear = {}
        phase_clear = {}
        hap_clear = {}
        for key, val in refinement.get('clear', {}).items():
            # Clear refinement options
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_clear[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_clear[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                hap_set[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            hist.clear_refinements(pwdr_clear)
        for phase in phases:
            phase.clear_refinements(phase_clear)
        for phase in phases:
            phase.clear_HAP_refinements(hap_clear, hists)

    def index_ids(self):
        import GSASIIstrIO as G2strIO
        self.save()
        return G2strIO.GetUsedHistogramsAndPhases(self.filename)

    def add_constraint_raw(self, cons_scope, constr):
        """Adds a constraint of type consType to the project.
        cons_scope should be one of "Hist", "Phase", "HAP", or "Global".

        WARNING it does not check the constraint is well-constructed"""
        constrs = self['Constraints']['data']
        if 'Global' not in constrs:
            constrs['Global'] = []
        constrs[cons_scope].append(constr)

    def hold_many(self, vars, type):
        """Apply holds for all the variables in vars, for constraint of a given type.

        type is passed directly to add_constraint_raw as consType

        :param list vars: A list of variables to hold. Either :class:`GSASIIobj.G2VarObj` objects,
            string variable specifiers, or arguments for :meth:`make_var_obj`
        :param str type: A string constraint type specifier. See
        :class:`~GSASIIscriptable.G2Project.add_constraint_raw`"""
        for var in vars:
            if isinstance(var, str):
                var = self.make_var_obj(var)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var)
            self.add_constraint_raw(type, [[1.0, var], None, None, 'h'])

    def make_var_obj(self, phase=None, hist=None, varname=None, atomId=None,
                     reloadIdx=True):
        """Wrapper to create a G2VarObj. Takes either a string representaiton ("p:h:name:a")
        or individual names of phase, histogram, varname, and atomId.

        Automatically converts string phase, hist, or atom names into the ID required
        by G2VarObj."""

        if reloadIdx:
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
        elif phase in self['Phases']:
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


class G2AtomRecord(G2ObjectWrapper):
    """Wrapper for an atom record. Has convenient accessors via @property.


    Available accessors: label, type, refinement_flags, coordinates,
        occupancy, ranId, id, adp_flag, uiso

    Example:

    >>> atom = some_phase.atom("O3")
    >>> # We can access the underlying data format
    >>> atom.data
    ['O3', 'O-2', '', ... ]
    >>> # We can also use wrapper accessors
    >>> atom.coordinates
    (0.33, 0.15, 0.5)
    >>> atom.refinement_flags
    u'FX'
    >>> atom.ranId
    4615973324315876477
    >>> atom.occupancy
    1.0
    """
    def __init__(self, data, indices, proj):
        self.data = data
        self.cx, self.ct, self.cs, self.cia = indices
        self.proj = proj

    @property
    def label(self):
        return self.data[self.ct-1]

    @property
    def type(self):
        return self.data[self.ct]

    @property
    def refinement_flags(self):
        return self.data[self.ct+1]

    @refinement_flags.setter
    def refinement_flags(self, other):
        # Automatically check it is a valid refinement
        for c in other:
            if c not in ' FXU':
                raise ValueError("Invalid atom refinement: ", other)
        self.data[self.ct+1] = unicode(other)

    @property
    def coordinates(self):
        return tuple(self.data[self.cx:self.cx+3])

    @property
    def occupancy(self):
        return self.data[self.cx+3]

    @occupancy.setter
    def occupancy(self, val):
        self.data[self.cx+3] = float(val)

    @property
    def ranId(self):
        return self.data[self.cia+8]

    @property
    def adp_flag(self):
        # Either 'I' or 'A'
        return self.data[self.cia]

    @property
    def uiso(self):
        if self.adp_flag == 'I':
            return self.data[self.cia+1]
        else:
            return self.data[self.cia+2:self.cia+8]


class G2PwdrData(G2ObjectWrapper):
    """Wraps a Powder Data Histogram."""
    def __init__(self, data, proj):
        self.data = data
        self.proj = proj

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ['Limits', 'Sample Parameters', 'Background',
                      'Instrument Parameters']
        return key in valid_keys

    @property
    def name(self):
        return self['data'][-1]

    @property
    def ranId(self):
        return self['data'][0]['ranId']

    @property
    def residuals(self):
        data = self['data'][0]
        return {key: data[key]
                for key in ['R', 'Rb', 'wR', 'wRb', 'wRmin']}

    @property
    def id(self):
        return G2obj.HistRanIdLookup[self.ranId]

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
                if parm in ['U','V','W','X','Y','SH/L','I(L2)/I(L1)','alpha',
                    'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',] and Inst[parm][2]:
                        Inst[parm][2] = False
            instDict = dict(zip(insNames, insVals))
            instDict['X'] = max(instDict['X'], 0.01)
            instDict['Y'] = max(instDict['Y'], 0.01)
            if 'SH/L' in instDict:
                instDict['SH/L'] = max(instDict['SH/L'], 0.002)
            return dataType, instDict, insVary

        bgrnd = self['Background']

        # Need our fixed points in order
        bgrnd[1]['FixedPoints'].sort(key=lambda pair: pair[0])
        X = [x for x, y in bgrnd[1]['FixedPoints']]
        Y = [y for x, y in bgrnd[1]['FixedPoints']]

        limits = self['Limits'][1]
        if X[0] > limits[0]:
            X = [limits[0]] + X
            Y = [Y[0]] + Y
        if X[-1] < limits[1]:
            X += [limits[1]]
            Y += [Y[-1]]

        # Some simple lookups
        controls = self.proj['Controls']['data']
        inst, inst2 = self['Instrument Parameters']
        pwddata = self['data'][1]

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

    def y_calc(self):
        return self['data'][1][3]

    def plot(self, Yobs=True, Ycalc=True, Background=True, Residual=True):
        try:
            import matplotlib.pyplot as plt
            data = self['data'][1]
            if Yobs:
                plt.plot(data[0], data[1], label='Yobs')
            if Ycalc:
                plt.plot(data[0], data[3], label='Ycalc')
            if Background:
                plt.plot(data[0], data[4], label='Background')
            if Residual:
                plt.plot(data[0], data[5], label="Residual")
        except ImportError:
            pass

    def set_refinements(self, refs):
        """Sets the refinement parameter 'key' to the specification 'value'

        :param dict refs: A dictionary of the parameters to be set.

        :returns: None"""
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

            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    try:
                        instrument[iparam][2] = True
                    except IndexError:
                        raise ValueError("Invalid key:", iparam)
            else:
                raise ValueError("Unknown key:", key)
        # Fit fixed points after the fact - ensure they are after fixed points
        # are added
        if do_fit_fixed_points:
            # Background won't be fit if refinement flag not set
            orig = self['Background'][0][1]
            self['Background'][0][1] = True
            self.fit_fixed_points()
            # Restore the previous value
            self['Background'][0][1] = orig

    def clear_refinements(self, refs):
        """Clears the refinement parameter 'key' and its associated value.

        :param dict refs: A dictionary of parameters to clear."""
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
            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    instrument[iparam][2] = False
            else:
                raise ValueError("Unknown key:", key)


class G2Phase(G2ObjectWrapper):
    """A wrapper object around a given phase.

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
                      "Pref.Ori.", "Show", "Size", "Use", "Scale"]
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
        """Returns a list of atoms present in the phase.

        :returns: A list of :class:`G2AtomRecord` objects.

        .. seealso::
            :func:`~GSASIIscriptable.G2Phase.atom`
            :class:`G2AtomRecord`
        """
        ptrs = self.data['General']['AtomPtrs']
        output = []
        atoms = self.data['Atoms']
        for atom in atoms:
            output.append(G2AtomRecord(atom, ptrs, self.proj))
        return output

    def histograms(self):
        output = []
        for hname in self.data.get('Histograms', {}).keys():
            output.append(self.proj.histogram(hname))
        return output

    @property
    def ranId(self):
        return self.data['ranId']

    @property
    def id(self):
        return G2obj.PhaseRanIdLookup[self.ranId]

    def get_cell(self):
        """Returns a dictionary of the cell parameters, with keys:
            'length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma', 'volume'

        :returns: a dict"""
        cell = self['General']['Cell']
        return {'length_a': cell[1], 'length_b': cell[2], 'length_c': cell[3],
                'angle_alpha': cell[4], 'angle_beta': cell[5], 'angle_gamma': cell[6],
                'volume': cell[7]}

    def export_CIF(self, outputname, quickmode=True):
        """Write this phase to a .cif file named outputname

        :param str outputname: The name of the .cif file to write to
        :param bool quickmode: Currently ignored. Carryover from exports.G2export_CIF"""
        # This code is all taken from exports/G2export_CIF.py
        # Functions copied have the same names
        import GSASIImath as G2mth
        import GSASIImapvars as G2mv
        from exports import G2export_CIF as cif

        CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        CIFname = os.path.splitext(self.proj.filename)[0]
        CIFname = os.path.split(CIFname)[1]
        CIFname = ''.join([c if ord(c) < 128 else ''
                           for c in CIFname.replace(' ', '_')])
        try:
            author = self.proj['Controls'].get('Author','').strip()
        except KeyError:
            pass
        oneblock = True

        covDict = self.proj['Covariance']
        parmDict = dict(zip(covDict.get('varyList',[]),
                            covDict.get('variables',[])))
        sigDict = dict(zip(covDict.get('varyList',[]),
                           covDict.get('sig',[])))

        if covDict.get('covMatrix') is not None:
            sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],
                                              covDict['varyList'],
                                              parmDict))

        with open(outputname, 'w') as fp:
            fp.write(' \n' + 70*'#' + '\n')
            cif.WriteCIFitem(fp, 'data_' + CIFname)
            # from exports.G2export_CIF.WritePhaseInfo
            cif.WriteCIFitem(fp, '\n# phase info for '+str(self.name) + ' follows')
            cif.WriteCIFitem(fp, '_pd_phase_name', self.name)
            # TODO get esds
            cellDict = self.get_cell()
            defsigL = 3*[-0.00001] + 3*[-0.001] + [-0.01] # significance to use when no sigma
            names = ['length_a','length_b','length_c',
                     'angle_alpha','angle_beta ','angle_gamma',
                     'volume']
            for key, val in cellDict.items():
                cif.WriteCIFitem(fp, '_cell_' + key, G2mth.ValEsd(val))

            cif.WriteCIFitem(fp, '_symmetry_cell_setting',
                         self['General']['SGData']['SGSys'])

            spacegroup = self['General']['SGData']['SpGrp'].strip()
            # regularize capitalization and remove trailing H/R
            spacegroup = spacegroup[0].upper() + spacegroup[1:].lower().rstrip('rh ')
            cif.WriteCIFitem(fp, '_symmetry_space_group_name_H-M', spacegroup)

            # generate symmetry operations including centering and center of symmetry
            SymOpList, offsetList, symOpList, G2oprList, G2opcodes = G2spc.AllOps(
                self['General']['SGData'])
            cif.WriteCIFitem(fp, 'loop_\n    _space_group_symop_id\n    _space_group_symop_operation_xyz')
            for i, op in enumerate(SymOpList,start=1):
                cif.WriteCIFitem(fp, '   {:3d}  {:}'.format(i,op.lower()))

            # TODO skipped histograms, exports/G2export_CIF.py:880

            # report atom params
            if self['General']['Type'] in ['nuclear','macromolecular']:        #this needs macromolecular variant, etc!
                cif.WriteAtomsNuclear(fp, self.data, self.name, parmDict, sigDict, [])
                # self._WriteAtomsNuclear(fp, parmDict, sigDict)
            else:
                raise Exception,"no export for "+str(self['General']['Type'])+" coordinates implemented"
            # report cell contents
            cif.WriteComposition(fp, self.data, self.name, parmDict)
            if not quickmode and self['General']['Type'] == 'nuclear':      # report distances and angles
                # WriteDistances(fp,self.name,SymOpList,offsetList,symOpList,G2oprList)
                raise NotImplementedError("only quickmode currently supported")
            if 'Map' in self['General'] and 'minmax' in self['General']['Map']:
                cif.WriteCIFitem(fp,'\n# Difference density results')
                MinMax = self['General']['Map']['minmax']
                cif.WriteCIFitem(fp,'_refine_diff_density_max',G2mth.ValEsd(MinMax[0],-0.009))
                cif.WriteCIFitem(fp,'_refine_diff_density_min',G2mth.ValEsd(MinMax[1],-0.009))


    def set_refinements(self, refs):
        """Sets the refinement parameter 'key' to the specification 'value'

        :param dict refs: A dictionary of the parameters to be set.

        :returns: None"""
        for key, value in refs.items():
            if key == "Cell":
                self.data['General']['Cell'][0] = True
            elif key == "Atoms":
                cx, ct, cs, cia = self.data['General']['AtomPtrs']

                for atomlabel, atomrefinement in value.items():
                    if atomlabel == 'all':
                        for atom in self.atoms():
                            atom.refinement_flags = atomrefinement
                    else:
                        atom = self.atom(atomlabel)
                        atom.refinement_flags = atomrefinement
            elif key == "LeBail":
                hists = self.data['Histograms']
                for hname, hoptions in hists.items():
                    if 'LeBail' not in hoptions:
                        hoptions['newLeBail'] = True
                    hoptions['LeBail'] = bool(value)
            else:
                raise ValueError("Unknown key:", key)

    def clear_refinements(self, refs):
        """Clears a given set of parameters.

        :param dict refs: The parameters to clear"""
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
                    if 'LeBail' not in hoptions:
                        hoptions['newLeBail'] = True
                    hoptions['LeBail'] = False
            else:
                raise ValueError("Unknown key:", key)

    def set_HAP_refinements(self, refs, histograms='all'):
        """Sets the given HAP refinement parameters between this phase and
        the given histograms

        :param dict refs: A dictionary of the parameters to be set.
        :param histograms: Either 'all' (default) or a list of the histograms
            whose HAP parameters will be set with this phase. Histogram and phase
            must already be associated

        :returns: None
        """
        if histograms == 'all':
            histograms = self['Histograms'].values()
        else:
            histograms = [self['Histograms'][h.name] for h in histograms]

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
                        for hist in histograms:
                            hist['Babinet'][param][1] = True
                elif key == 'Extinction':
                    for h in histograms:
                        h['Extinction'][1] = True
                elif key == 'HStrain':
                    for h in histograms:
                        hist['HStrain'][1] = [True for p in hist['Hstrain'][0]]
                elif key == 'Mustrain':
                    for h in histograms:
                        mustrain = h['Mustrain']
                        newType = None
                        if isinstance(val, str):
                            if val in ['isotropic', 'uniaxial', 'generalized']:
                                newType = val
                            else:
                                raise ValueError("Not a Mustrain type: " + val)
                        elif isinstance(val, dict):
                            newType = val.get('type', None)
                            direction = val.get('direction', None)

                        if newType:
                            mustrain[0] = newType
                            if newType == 'isotropic':
                                mustrain[2][0] = True
                                mustrain[5] = [False for p in mustrain[4]]
                            elif newType == 'uniaxial':
                                pass
                            else:  # newtype == 'generalized'
                                mustrain[2] = [False for p in mustrain[1]]
                        if direction:
                            direction = [int(n) for n in direction]
                            if len(direction) != 3:
                                raise ValueError("Expected hkl, found", direction)
                            mustrain[3] = direction
                elif key == 'Pref.Ori.':
                    for h in histograms:
                        h['Pref.Ori.'][2] = True
                elif key == 'Show':
                    for h in histograms:
                        h['Show'] = True
                elif key == 'Size':
                    raise NotImplementedError()
                elif key == 'Use':
                    for h in histograms:
                        h['Use'] = True
                elif key == 'Scale':
                    for h in histograms:
                        h['Scale'][1] = False

    def clear_HAP_refinements(self, refs, histograms='all'):
        """Clears the given HAP refinement parameters between this phase and
        the given histograms

        :param dict refs: A dictionary of the parameters to be cleared.
        :param histograms: Either 'all' (default) or a list of the histograms
            whose HAP parameters will be cleared with this phase. Histogram and
            phase must already be associated

        :returns: None
        """
        if histograms == 'all':
            histograms = self['Histograms'].values()
        else:
            histograms = [self['Histograms'][h.name] for h in histograms]

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
                        for hist in histograms:
                            hist['Babinet'][param][1] = False
                elif key == 'Extinction':
                    for h in histograms:
                        h['Extinction'][1] = False
                elif key == 'HStrain':
                    for h in histograms:
                        hist['HStrain'][1] = [False for p in hist['Hstrain'][0]]
                elif key == 'Mustrain':
                    for h in histograms:
                        mustrain = h['Mustrain']
                        mustrain[2] = [False for p in mustrain[1]]
                        mustrain[5] = [False for p in mustrain[4]]
                elif key == 'Pref.Ori.':
                    for h in histograms:
                        h['Pref.Ori.'][2] = False
                elif key == 'Show':
                    for h in histograms:
                        h['Show'] = False
                elif key == 'Size':
                    raise NotImplementedError()
                elif key == 'Use':
                    for h in histograms:
                        h['Use'] = False
                elif key == 'Scale':
                    for h in histograms:
                        h['Scale'][1] = False


##########################
# Command Line Interface #
##########################


# TODO SUBPARSERS

def create(*args):
    """The create subcommand.

    Should be passed all the command-line arguments after `create`"""
    import argparse
    parser = argparse.ArgumentParser(prog=' '.join([sys.argv[0], sys.argv[1]]))
    parser.add_argument('filename',
                        help='the project file to create. should end in .gpx')
    parser.add_argument('-g', '--histograms', nargs='+',
                        help='list of histograms to add')
    parser.add_argument('-p', '--phases', nargs='+',
                        help='list of phases to add')
    results = parser.parse_args(args)

    proj = G2Project(filename=filename)

    isPhase = False
    isPowderData = False
    isInstPrms = False
    instPrms = None

    # TODO how to associate phase with histogram?
    for arg in args[1:]:
        if arg == '--phases':
            isPhase = True
            isPowderData = False
            isInstPrms = False
        elif arg == '--powder':
            isPhase = False
            isPowderData = True
            isInstPrms = False
        # Instrument parameters must be specified before
        # any powder data files are passed
        elif arg == '--iparams':
            isPhase = False
            isPowderData = False
            isInstPrms = True
        elif isPhase:
            proj.add_phase(arg)
        elif isPowderData:
            proj.add_powder_histogram(arg, instPrms)
        elif isInstPrms:
            instPrms = arg
            isInstPrms = False
        else:
            print("Not sure what to do with: {}".format(arg))

    proj.save()


def dump(*args):
    """The dump subcommand"""
    import argparse
    parser = argparse.ArgumentParser(prog=' '.join([sys.argv[0], sys.argv[1]]))
    parser.add_argument('-g', '--histograms', action='store_true',
                        help='list histograms in files, overrides --raw')
    parser.add_argument('-p', '--phases', action='store_true',
                        help='list phases in files, overrides --raw')
    parser.add_argument('-r', '--raw', action='store_true',
                        help='dump raw file contents')
    parser.add_argument('files', nargs='*')
    results = parser.parse_args(args)

    if not results.histograms and not results.phases:
        results.raw = True
    if results.raw:
        import IPython.lib.pretty as pretty

    for fname in results.files:
        if results.raw:
            proj, nameList = LoadDictFromProjFile(fname)
            print("file:", fname)
            print("names:", nameList)
            for key, val in proj.items():
                print(key, ":")
                pretty.pprint(val)
        else:
            proj = G2Project(fname)
            if results.histograms:
                hists = proj.histograms()
                for h in hists:
                    print(fname, "hist", h.id, h.name)
            if results.phases:
                phase_list = proj.phases()
                for p in phase_list:
                    print(fname, "phase", p.id, p.name)


def IPyBrowse(*args):
    """Load a .gpx file and then open a IPython shell to browse it
    """
    filename = []
    for arg in args:
        fname = arg
        proj, nameList = LoadDictFromProjFile(fname)
        msg = "\nfile {} loaded into proj (dict) with names in nameList".format(fname)
        GSASIIpath.IPyBreak_base(msg)
        break


def refine(*args):
    """The refine subcommand"""
    proj = G2Project(args[0])
    if len(args) == 1:
        proj.refine()
    elif len(args) == 2:
        import json
        with open(args[1]) as refs:
            refs = json.load(refs)
        proj.do_refinements(refs['refinements'])
    else:
        print("Refine not sure what to do with args:", args)


def seqrefine(*args):
    """The seqrefine subcommand"""
    raise NotImplementedError("seqrefine is not yet implemented")


def export(*args):
    """The export subcommand"""
    # Export CIF or Structure or ...
    gpxfile, phase, exportfile = args
    proj = G2Project(gpxfile)
    phase = proj.phase(phase)
    phase.export_CIF(exportfile)


subcommands = {"create": create,
               "dump": dump,
               "refine": refine,
               "seqrefine": seqrefine,
               "export": export,
               "browse": IPyBrowse}


def main():
    '''The command line interface for GSASIIscriptable.

    Executes one of the following subcommands:

        * :func:`create`
        * :func:`dump`
        * :func:`refine`
        * :func:`seqrefine`
        * :func:`export`
        * :func:`browse`

    .. seealso::
        :func:`create`
        :func:`dump`
        :func:`refine`
        :func:`seqrefine`
        :func:`export`
        :func:`browse`
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('subcommand', choices=sorted(subcommands.keys()),
                        help='The subcommand to be executed')

    result = parser.parse_args(sys.argv[1:2])
    sub = result.subcommand
    subcommands[sub](*sys.argv[2:])

    # argv = sys.argv
    # if len(argv) > 1 and argv[1] in subcommands:
    #     subcommands[argv[1]](*argv[2:])
    # elif len(argv) == 1 or argv[1] in ('help', '--help', '-h'):
    #     # TODO print usage
    #     subcommand_names = ' | '.join(sorted(subcommands.keys()))
    #     print("USAGE: {} [ {} ] ...".format(argv[0], subcommand_names))
    # else:
    #     print("Unknown subcommand: {}".format(argv[1]))
    #     print("Available subcommands:")
    #     for name in sorted(subcommands.keys()):
    #         print("\t{}".format(name))
    #     sys.exit(-1)
    # sys.exit(0)

    # arg = sys.argv
    # print(arg)
    # if len(arg) > 1:
    #     GPXfile = arg[1]
    #     if not ospath.exists(GPXfile):
    #         print(u'ERROR - '+GPXfile+u" doesn't exist!")
    #         exit()
    #     Project,nameList = LoadDictFromProjFile(GPXfile)
    #     SaveDictToProjFile(Project,nameList,'testout.gpx')
    # else:
    #     print('ERROR - missing filename')
    #     exit()
    # print("Done")

if __name__ == '__main__':
    main()

    # from gpx_manipulatons.py
    # try:
    #     filename, authorname = sys.argv[1:3]
    #     proj, names = make_empty_project(authorname, filename)
    #     SaveDictToProjFile(proj, names, os.path.abspath(filename))
    # except ValueError:
    #     print("Usage: {} <filename> <author>".format(sys.argv[0]))
    #     sys.exit(-1)


    # from refinements.py
#      USAGE = """USAGE: {} datafile instparams phasefile projectname refinements

# datafile:    Input powder data
# intparams:   Corresponding instrument parameter file
# phasefile:   Phase to refine against data
# projectname: Project file to be created, should end in .gpx
# refinements: JSON file of refinements to be executed
# """
#     try:
#         datafile, instprm, phasefile, projectname, refinements = sys.argv[1:]
#     except ValueError:
#         print(USAGE.format(sys.argv[0]))
#         sys.exit(-1)

#     try:
#         with open(refinements) as f:
#             refinements = json.load(f)
#     except IOError:
#         print("No such refinements file: {}".format(refinements))

#     print("Creating project file \"{}\"...".format(projectname))
#     proj = G2Project(filename=projectname)
#     # Add the histogram
#     hist = proj.add_powder_histogram(datafile, instprm)
#     # Add the phase, and associate it with the histogram
#     proj.add_phase(phasefile, histograms=[hist.name])

#     proj.do_refinements(refinements['refinements'])
#     proj.save()


    # from gpx_dumper
    # import IPython.lib.pretty as pretty
    # proj, nameList = LoadDictFromProjFile(sys.argv[1])
    # print("names:", nameList)
    # for key, val in proj.items():
    #     print(key, ":", sep='')
    #     pretty.pprint(val)
