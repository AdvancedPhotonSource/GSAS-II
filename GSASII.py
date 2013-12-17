#!/usr/bin/env python
# -*- coding: utf-8 -*-
#GSASII
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSAS-II Main Module*
=====================

Main routines for the GSAS-II program
'''

import os
import sys
import math
import copy
import random as ran
import time
import copy
import glob
import imp
import inspect
import numpy as np
import scipy as sp
import wx
import matplotlib as mpl
try:
    import OpenGL as ogl
except ImportError:
    print('*******************************************************')
    print('PyOpenGL is missing from your python installation')
    print('     - we will try to install it')
    print('*******************************************************')
    def install_with_easyinstall(package):
        try: 
            print "trying a system-wide PyOpenGl install"
            easy_install.main(['-f',os.path.split(__file__)[0],package])
            return
        except:
            pass
        try: 
            print "trying a user level PyOpenGl install"
            easy_install.main(['-f',os.path.split(__file__)[0],'--user',package])
            return
        except:
            print "Install of '+package+' failed. Please report this information:"
            import traceback
            print traceback.format_exc()
            sys.exit()
    from setuptools.command import easy_install
    install_with_easyinstall('PyOpenGl')
    print('*******************************************************')         
    print('OpenGL has been installed. Please restart GSAS-II')
    print('*******************************************************')         
    sys.exit()

# load the GSAS routines
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import GSASIIplot as G2plt
import GSASIIpwd as G2pwd
import GSASIIpwdGUI as G2pdG
import GSASIIspc as G2spc
import GSASIIstrMain as G2stMn
import GSASIIstrIO as G2stIO
import GSASIImapvars as G2mv
import GSASIIobj as G2obj
import GSASIIlattice as G2lat

#wx inspector - use as needed
wxInspector = False

# print versions
print "Python module versions loaded:"
print "python:     ",sys.version[:5]
print "wxpython:   ",wx.__version__
print "matplotlib: ",mpl.__version__
print "numpy:      ",np.__version__
print "scipy:      ",sp.__version__
print "OpenGL:     ",ogl.__version__
try:
    import mkl
    print "Max threads ",mkl.get_max_threads()
except:
    pass
#    print "MKL module not present"
__version__ = '0.2.0'
G2gd.__version__ = __version__
print "This is GSAS-II version:     ",__version__,' revision '+str(GSASIIpath.GetVersionNumber())

# useful degree trig functions
sind = lambda x: math.sin(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
acosd = lambda x: 180.*math.acos(x)/math.pi
atan2d = lambda x,y: 180.*math.atan2(y,x)/math.pi

def create(parent):
    return GSASII(parent)

class GSASII(wx.Frame):
    '''Define the main GSAS-II frame and its associated menu items
    '''
    def _Add_FileMenuItems(self, parent):
        item = parent.Append(
            help='Open a GSAS-II project file (*.gpx)', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Open project...')
        self.Bind(wx.EVT_MENU, self.OnFileOpen, id=item.GetId())
        item = parent.Append(
            help='Save project under current name', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Save project')
        self.Bind(wx.EVT_MENU, self.OnFileSave, id=item.GetId())
        item = parent.Append(
            help='Save current project to new file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Save project as...')
        self.Bind(wx.EVT_MENU, self.OnFileSaveas, id=item.GetId())
        item = parent.Append(
            help='Create empty new project, saving current is optional', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&New project')
        self.Bind(wx.EVT_MENU, self.OnFileClose, id=item.GetId())
        item = parent.Append(
            help='Exit from GSAS-II', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='&Exit')
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=item.GetId())
        
    def _Add_DataMenuItems(self,parent):
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Read image data...')
        self.Bind(wx.EVT_MENU, self.OnImageRead, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Read Powder Pattern Peaks...')
        self.Bind(wx.EVT_MENU, self.OnReadPowderPeaks, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Sum powder data')
        self.Bind(wx.EVT_MENU, self.OnPwdrSum, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Sum image data')
        self.Bind(wx.EVT_MENU, self.OnImageSum, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Add phase')
        self.Bind(wx.EVT_MENU, self.OnAddPhase, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Delete phase')
        self.Bind(wx.EVT_MENU, self.OnDeletePhase, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Rename data') 
        self.Bind(wx.EVT_MENU, self.OnRenameData, id=item.GetId())
        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Delete data')
        self.Bind(wx.EVT_MENU, self.OnDataDelete, id=item.GetId())
                
    def _Add_CalculateMenuItems(self,parent):
        item = parent.Append(help='Make new PDFs from selected powder patterns', 
            id=wx.ID_ANY, kind=wx.ITEM_NORMAL,text='Make new PDFs')
        self.MakePDF.append(item)
#        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnMakePDFs, id=item.GetId())
        
        item = parent.Append(help='View least squares parameters', 
            id=wx.ID_ANY, kind=wx.ITEM_NORMAL,text='&View LS parms')
        self.Bind(wx.EVT_MENU, self.ShowLSParms, id=item.GetId())
        
        item = parent.Append(help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
            text='&Refine')
        if len(self.Refine): # extend state for new menus to match main (on mac)
            state = self.Refine[0].IsEnabled()
        else:
            state = False
        item.Enable(state)
        self.Refine.append(item)
        self.Bind(wx.EVT_MENU, self.OnRefine, id=item.GetId())
        
        item = parent.Append(help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
            text='Sequential refine')
        if len(self.SeqRefine): # extend state for new menus to match main (on mac)
            state = self.SeqRefine[0].IsEnabled()
        else:
            state = False
        item.Enable(state)
        self.SeqRefine.append(item)
        self.Bind(wx.EVT_MENU, self.OnSeqRefine, id=item.GetId())
        
    def _init_Imports(self):
        '''import all the G2phase*.py & G2sfact*.py & G2pwd*.py files that 
        are found in the path
        '''

        self.ImportPhaseReaderlist = []
        self._init_Import_routines('phase',self.ImportPhaseReaderlist,'Phase')
        self.ImportSfactReaderlist = []
        self._init_Import_routines('sfact',self.ImportSfactReaderlist,'Struct_Factor')
        self.ImportPowderReaderlist = []
        self._init_Import_routines('pwd',self.ImportPowderReaderlist,'Powder_Data')
        self.ImportMenuId = {}

    def _init_Import_routines(self,prefix,readerlist,errprefix):
        '''import all the import readers matching the prefix
        '''
        #path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # location of this file
        #pathlist = sys.path[:]
        #if path2GSAS2 not in pathlist: pathlist.append(path2GSAS2)
        #path2GSAS2 = os.path.join(
        #    os.path.dirname(os.path.realpath(__file__)), # location of this file
        #    'imports')
        pathlist = sys.path[:]
        #if path2GSAS2 not in pathlist: pathlist.append(path2GSAS2)

        filelist = []
        for path in pathlist:
            for filename in glob.iglob(os.path.join(
                path,
                "G2"+prefix+"*.py")):
                filelist.append(filename)    
                #print 'debug: found',filename
        filelist = sorted(list(set(filelist))) # remove duplicates
        for filename in filelist:
            path,rootname = os.path.split(filename)
            pkg = os.path.splitext(rootname)[0]
            try:
                fp = None
                fp, fppath,desc = imp.find_module(pkg,[path,])
                pkg = imp.load_module(pkg,fp,fppath,desc)
                for clss in inspect.getmembers(pkg): # find classes defined in package
                    if clss[0].startswith('_'): continue
                    if inspect.isclass(clss[1]):
                        # check if we have the required methods
                        for m in 'Reader','ExtensionValidator','ContentsValidator':
                            if not hasattr(clss[1],m): break
                            if not callable(getattr(clss[1],m)): break
                        else:
                            reader = clss[1]() # create an import instance
                            readerlist.append(reader)
            except AttributeError:
                print 'Import_'+errprefix+': Attribute Error'+str(filename)
                pass
            except ImportError:
                print 'Import_'+errprefix+': Error importing file'+str(filename)
                pass
            if fp: fp.close()

    def OnImportGeneric(self,reader,readerlist,label,multiple=False,usedRanIdList=[]):
        '''Used to import Phases, powder dataset or single
        crystal datasets (structure factor tables) using reader objects
        subclassed from :class:`GSASIIIO.ImportPhase`, :class:`GSASIIIO.ImportStructFactor`
        or :class:`GSASIIIO.ImportPowderData`. If a reader is specified, only
        that will be attempted, but if no reader is specified, every one
        that is potentially compatible (by file extension) will
        be tried on the selected file(s). 

        :param readerobject reader: This will be a reference to
          a particular object to be used to read a file or None,
          if every appropriate reader should be used.

        :param list readerlist: a list of reader objects appropriate for
          the current read attempt. At present, this will be either
          self.ImportPhaseReaderlist, self.ImportSfactReaderlist or
          self.ImportPowderReaderlist (defined in _init_Imports from
          the files found in the path), but in theory this list could
          be tailored. Used only when reader is None.

        :param str label: string to place on the open file dialog:
          Open `label` input file

        :param bool multiple: True if multiple files can be selected
          in the file dialog. False is default. At present True is used
          only for reading of powder data.

        :param list usedRanIdList: an optional list of random Ids that 
          have been used and should not be reused

        :returns: a list of reader objects (rd_list) that were able
          to read the specified file(s). This list may be empty. 
        '''
        self.lastimport = ''
        self.zipfile = None
        singlereader = True
        if reader is None:
            singlereader = False
            multiple = False
            #print "use all formats"
            choices = "any file (*.*)|*.*"
            choices += "|zip archive (.zip)|*.zip"
            extdict = {}
            # compile a list of allowed extensions
            for rd in readerlist:
                fmt = rd.formatName
                for extn in rd.extensionlist:
                    if not extdict.get(extn): extdict[extn] = []
                    extdict[extn] += [fmt,]
            for extn in sorted(extdict.keys(),cmp=lambda x,y: cmp(x.lower(), y.lower())):
                fmt = ''
                for f in extdict[extn]:
                    if fmt != "": fmt += ', '
                    fmt += f
                choices += "|" + fmt + " file (*" + extn + ")|*" + extn
        else:
            readerlist = [reader,]
            # compile a list of allowed extensions
            choices = reader.formatName + " file ("
            w = ""
            for extn in reader.extensionlist:
                if w != "": w += ";"
                w += "*" + extn
            choices += w + ")|" + w
            choices += "|zip archive (.zip)|*.zip"
            if not reader.strictExtension:
                choices += "|any file (*.*)|*.*"
        # get the file(s)
        if multiple:
            mode = style=wx.OPEN | wx.CHANGE_DIR | wx.MULTIPLE
        else:
            mode = style=wx.OPEN | wx.CHANGE_DIR
        dlg = wx.FileDialog(self, message="Choose "+label+" input file",
            defaultFile="",wildcard=choices, style=mode)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                if multiple:
                    filelist = dlg.GetPaths()
                    if len(filelist) == 0: return []
                else:
                    filename = dlg.GetPath()
                    filelist = [filename,]
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
            else: # cancel was pressed
                return []
        finally:
            dlg.Destroy()
        rd_list = []
        filelist1 = []
        for filename in filelist:
            # is this a zip file?
            if os.path.splitext(filename)[1].lower() == '.zip':
                extractedfiles = G2IO.ExtractFileFromZip(
                    filename,parent=self,
                    multipleselect=True)
                if extractedfiles is None: continue # error or Cancel 
                if extractedfiles != filename:
                    self.zipfile = filename # save zip name
                    filelist1 += extractedfiles
                    continue
            filelist1.append(filename)
        filelist = filelist1
        for filename in filelist:
            # is this a zip file?
            if os.path.splitext(filename)[1].lower() == '.zip':
                extractedfile = G2IO.ExtractFileFromZip(filename,parent=self)
                if extractedfile is None: continue # error or Cancel 
                if extractedfile != filename:
                    filename,self.zipfile = extractedfile,filename # now use the file that was created
            # set what formats are compatible with this file
            primaryReaders = []
            secondaryReaders = []
            for rd in readerlist:
                flag = rd.ExtensionValidator(filename)
                if flag is None: 
                    secondaryReaders.append(rd)
                elif flag:
                    primaryReaders.append(rd)
            if len(secondaryReaders) + len(primaryReaders) == 0:
                self.ErrorDialog('No Format','No matching format for file '+filename)
                return []

            fp = None
            msg = ''
            try:
                fp = open(filename,'Ur')
                if len(filelist) == 1:
                    # confirm we have the right file
                    rdmsg = 'File '+str(filename)+' begins:\n\n'
                    for i in range(3):
                        rdmsg += fp.readline()
                    rdmsg += '\n\nDo you want to read this file?'
                    if not all([ord(c) < 128 and ord(c) != 0 for c in rdmsg]): # show only if ASCII
                        rdmsg = 'File '+str(
                            filename)+' is a binary file. Do you want to read this file?'
                    result = wx.ID_NO
                    # it would be better to use something that
                    # would resize better, but this will do for now
                    dlg = wx.MessageDialog(
                        self, rdmsg,
                        'Is this the file you want?', 
                        wx.YES_NO | wx.ICON_QUESTION,
                        )
                    dlg.SetSize((700,300)) # does not resize on Mac
                    try:
                        result = dlg.ShowModal()
                    finally:
                        dlg.Destroy()
                    if result == wx.ID_NO: return []
                            
                self.lastimport = filename
                # try the file first with Readers that specify the
                # file's extension and later with ones that merely allow it
                errorReport = ''
                for rd in primaryReaders+secondaryReaders:
                    rd.ReInitialize() # purge anything from a previous read
                    fp.seek(0)  # rewind
                    rd.errors = "" # clear out any old errors
                    if not rd.ContentsValidator(fp): # rejected on cursory check
                        errorReport += "\n  "+rd.formatName + ' validator error'
                        if rd.errors: 
                            errorReport += ': '+rd.errors
                        continue 
                    repeat = True
                    rdbuffer = {} # create temporary storage for file reader
                    block = 0
                    while repeat: # loop if the reader asks for another pass on the file
                        block += 1
                        repeat = False
                        fp.seek(0)  # rewind
                        rd.objname = os.path.basename(filename)
                        flag = False
                        try:
                            flag = rd.Reader(filename,fp,self,
                                             buffer=rdbuffer,
                                             blocknum=block,
                                             usedRanIdList=usedRanIdList,
                                             )
                        except rd.ImportException as detail:
                            rd.errors += "\n  Read exception: "+str(detail)
                        except Exception as detail:
                            rd.errors += "\n  Unhandled read exception: "+str(detail)
                            rd.errors += "\n  Traceback info:\n"+str(traceback.format_exc())
                        if flag: # this read succeeded
                            rd_list.append(copy.deepcopy(rd)) # save the result before it is written over
                            if rd.repeat:
                                repeat = True
                            continue
                        errorReport += '\n'+rd.formatName + ' read error'
                        if rd.errors:
                            errorReport += ': '+rd.errors
                    if rd_list: # read succeeded, was there a warning or any errors? 
                        if rd.warnings:
                            self.ErrorDialog('Read Warning','The '+ rd.formatName+
                                             ' reader reported a warning message:\n\n'+
                                             rd.warnings)
                        break # success in reading, try no further
                else:
                    if singlereader:
                        self.ErrorDialog('Read Error','The '+ rd.formatName+
                                         ' reader was not able to read file '+filename+msg+
                                         '\n\nError message(s):\n'+errorReport
                                         )
                    else:
                        self.ErrorDialog('Read Error','No reader was able to read file '+filename+msg+
                                         '\n\nError messages:\n'+errorReport
                                         )
            except:
                import traceback
                print traceback.format_exc()
                self.ErrorDialog('Open Error','Unexpected error trying to open or read file '+filename)
            if fp: fp.close()
        return rd_list

    def _Add_ImportMenu_Phase(self,parent):
        '''configure the Import Phase menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Phase',
            submenu, help='Import phase data')
        for reader in self.ImportPhaseReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPhase, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
                              help='Import phase data, use file to try to determine format',
                              kind=wx.ITEM_NORMAL,
                              text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportPhase, id=item.GetId())

    def OnImportPhase(self,event):
        '''Called in response to an Import/Phase/... menu item
        to read phase information.
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.
        '''
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        
        # make a list of phase names, ranId's and the histograms used in those phases
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        usedHKLFhists = [] # used single-crystal histograms
        for p in usedHistograms:
            for h in usedHistograms[p]:
                if h.startswith('HKLF ') and h not in usedHKLFhists:
                    usedHKLFhists.append(h)
        rdlist = self.OnImportGeneric(reqrdr,
                                  self.ImportPhaseReaderlist,
                                  'phase',usedRanIdList=phaseRIdList)
        if len(rdlist) == 0: return
        # for now rdlist is only expected to have one element
        # but below will allow multiple phases to be imported
        # if ever the import routines ever implement multiple phase reads. 
        self.CheckNotebook()
        newPhaseList = []
        for rd in rdlist:
            dlg = wx.TextEntryDialog( # allow editing of phase name
                self, 'Enter the name for the new phase',
                'Edit phase name', rd.Phase['General']['Name'],
                style=wx.OK)
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                rd.Phase['General']['Name'] = dlg.GetValue()
            dlg.Destroy()
            # make new phase names unique
            rd.Phase['General']['Name'] = G2obj.MakeUniqueLabel(rd.Phase['General']['Name'],phaseNameList)
            PhaseName = rd.Phase['General']['Name']
            newPhaseList.append(PhaseName)
            print 'Read phase '+str(PhaseName)+' from file '+str(self.lastimport)
            if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
                sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
            else:
                sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
            psub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
            self.PatternTree.SetItemPyData(psub,rd.Phase)
            self.PatternTree.Expand(self.root) # make sure phases are seen
            self.PatternTree.Expand(sub) 
            self.PatternTree.Expand(psub)
            self.PatternTree.UnselectAll()
            self.PatternTree.SelectItem(psub) # show the page to complete the initialization (yuk!)
            wx.Yield() # make sure call of GSASII.OnPatternTreeSelChanged happens before we go on

            if rd.Constraints:
                sub = G2gd.GetPatternTreeItemId(self,self.root,'Constraints') # was created in CheckNotebook if needed
                Constraints = self.PatternTree.GetItemPyData(sub)                
                # TODO: make sure that NEWVAR names are unique here?
                for i in rd.Constraints:
                    if type(i) is dict:
                        #for j in i: print j,' --> ',i[j]
                        if '_Explain' not in Constraints: Constraints['_Explain'] = {}
                        Constraints['_Explain'].update(i)
                        continue
                    Constraints['Phase'].append(i)
        if not newPhaseList: return # somehow, no new phases
        # get a list of existing histograms
        PWDRlist = []
        HKLFlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        TextList = PWDRlist + HKLFlist
        if not len(TextList): return # no data loaded yet
        header = 'Select histogram(s) to add to new phase(s):'
        for phaseName in newPhaseList:
            header += '\n  '+str(phaseName)

        notOK = True
        while notOK:
            result = G2gd.ItemSelector(TextList,self,header,header='Add histogram(s)',multiple=True)
            if not result: return
            # check that selected single crystal histograms are not already in use!
            used = [TextList[i] for i in result if TextList[i] in usedHKLFhists]
            #for i in result:
            #    if TextList[i] in usedHKLFhists: used.append(TextList[i])
            if used:
                msg = 'The following single crystal histogram(s) are already in use'
                for i in used:
                    msg += '\n  '+str(i)
                msg += '\nAre you sure you want to add them to this phase? '
                msg += 'Associating a single crystal dataset to >1 histogram is usually an error, '
                msg += 'so No is suggested here.'
                if self.ErrorDialog('Likely error',msg,self,wtype=wx.YES_NO) == wx.ID_YES: notOK = False
            else:
                notOK = False
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.PatternTree.GetFirstChild(sub)
        while item: # loop over (new) phases
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if phaseName not in newPhaseList: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for i in result:
                histoName = TextList[i]
                if histoName in HKLFlist:
                    UseList[histoName] = {
                        'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                        'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
                        'Extinction':['Lorentzian','None', {'Tbar':0.1,'Cos2TM':0.955,
                                                            'Eg':[1.e-10,False],
                                                            'Es':[1.e-10,False],
                                                            'Ep':[1.e-10,False]},]
                        }
                    #redo UpdateHKLFdata(histoName) here:
                    Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                    refDict,reflData = self.PatternTree.GetItemPyData(Id)
                    G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                    for iref,ref in enumerate(reflData['RefList']):
                        H = list(ref[:3])
                        ref[4] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                        iabsnt,ref[3],Uniq,phi = G2spc.GenHKLf(H,SGData)
                elif histoName in PWDRlist:
                    UseList[histoName] = {
                        'Histogram':histoName,'Show':False,
                        'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{}],
                        'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                                [1.,1.,1.,0.,0.,0.],6*[False,]],
                        'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                    NShkl*[0.01,],NShkl*[False,]],
                        'HStrain':[NDij*[0.0,],NDij*[False,]],                          
                        'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}
                    Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                    refList = self.PatternTree.GetItemPyData(
                        G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'))
                    refList[generalData['Name']] = []
                else:
                    raise Exception('Unexpected histogram '+str(histoName))
        wx.EndBusyCursor()
        return # success
        
    def _Add_ImportMenu_Sfact(self,parent):
        '''configure the Import Structure Factor menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Structure Factor',
            submenu, help='Import Structure Factor data')
        for reader in self.ImportSfactReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,                
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportSfact, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
            help='Import Structure Factor, use file to try to determine format',
            kind=wx.ITEM_NORMAL,
            text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportSfact, id=item.GetId())

    def OnImportSfact(self,event):
        '''Called in response to an Import/Structure Factor/... menu item
        to read single crystal datasets. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.
        '''
        # get a list of existing histograms
        HKLFlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        rdlist = self.OnImportGeneric(reqrdr,self.ImportSfactReaderlist,
            'Structure Factor')
        if len(rdlist) == 0: return
        self.CheckNotebook()
        newHistList = []
        for rd in rdlist:
            HistName = rd.objname
            if len(rdlist) <= 2: 
                dlg = wx.TextEntryDialog( # allow editing of Structure Factor name
                    self, 'Enter the name for the new Structure Factor',
                    'Edit Structure Factor name', HistName,
                    style=wx.OK)
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    HistName = dlg.GetValue()
                dlg.Destroy()
            HistName = 'HKLF '+HistName
            # make new phase names unique
            HistName = G2obj.MakeUniqueLabel(HistName,HKLFlist)
            print 'Read structure factor table '+str(HistName)+' from file '+str(self.lastimport)
            if not rd.RefDict.get('FF'):
                rd.RefDict['FF'] = [{} for i in range(len(rd.RefDict['RefList']))]
            Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxint),
                }
            self.PatternTree.SetItemPyData(Id,[valuesdict,rd.RefDict])
            Sub = self.PatternTree.AppendItem(Id,text='Instrument Parameters')
            self.PatternTree.SetItemPyData(Sub,rd.Parameters)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='HKL Plot Controls'),
                rd.Controls)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Reflection List'),[])  #dummy entry for GUI use
            self.PatternTree.SelectItem(Id)
            self.PatternTree.Expand(Id)
            self.Sngl = Id
            newHistList.append(HistName)

        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\nsingle crystal dataset(s) to:'
        for Name in newHistList:
            header += '\n  '+str(Name)

        notOK = True
        result = G2gd.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.PatternTree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for histoName in newHistList:
                UseList[histoName] = {
                    'Histogram':histoName,'Show':False,'Scale':[1.0,True],
                    'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
                    'Extinction':['Lorentzian','None', {'Tbar':0.1,'Cos2TM':0.955,
                                                        'Eg':[1.e-10,False],
                                                        'Es':[1.e-10,False],
                                                        'Ep':[1.e-10,False]},]
                    }
                #redo UpdateHKLFdata(histoName) here:
                Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                refDict,reflData = self.PatternTree.GetItemPyData(Id)
                G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                for iref,ref in enumerate(reflData['RefList']):
                    H = list(ref[:3])
                    ref[4] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                    iabsnt,ref[3],Uniq,phi = G2spc.GenHKLf(H,SGData)
        wx.EndBusyCursor()
        
        return # success

    def _Add_ImportMenu_powder(self,parent):
        '''configure the Powder Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendMenu(wx.ID_ANY, 'Powder Data',
            submenu, help='Import Powder data')
        for reader in self.ImportPowderReaderlist:
            item = submenu.Append(wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='from '+reader.formatName+' file')
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,
            help='Import powder data, use file to try to determine format',
            kind=wx.ITEM_NORMAL,text='guess format from file')
        self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        submenu.AppendSeparator()
        item = submenu.Append(wx.ID_ANY,
            help='Create a powder data set entry that will be simulated',
            kind=wx.ITEM_NORMAL,text='Simulate a dataset')
        self.Bind(wx.EVT_MENU, self.OnDummyPowder, id=item.GetId())
           
    def ReadPowderInstprm(self,instfile):       #fix the write routine for [inst1,inst2] style
        '''Read a GSAS-II (new) instrument parameter file

        :param str instfile: name of instrument parameter file

        '''
        if os.path.splitext(instfile)[1].lower() != '.instprm': # invalid file
            return None            
        if not os.path.exists(instfile): # no such file
            return None
        File = open(instfile,'r')
        S = File.readline()
        if not S.startswith('#GSAS-II'): # not a valid file
            File.close()
            return None
        newItems = []
        newVals = []
        while S:
            if S[0] == '#':
                S = File.readline()
                continue
            [item,val] = S[:-1].split(':')
            newItems.append(item)
            try:
                newVals.append(float(val))
            except ValueError:
                newVals.append(val)                        
            S = File.readline()                
        File.close()
        # add a second MT dict here. TOF parms? (BHT)
        return G2IO.makeInstDict(newItems,newVals,len(newVals)*[False,]),{}
        
    def ReadPowderIparm(self,instfile,bank,databanks,rd):
        '''Read a GSAS (old) instrument parameter file

        :param str instfile: name of instrument parameter file

        :param int bank: the bank number read in the raw data file

        :param int databanks: the number of banks in the raw data file.
          If the number of banks in the data and instrument parameter files
          agree, then the sets of banks are assumed to match up and bank
          is used to select the instrument parameter file. If not, the user
          is asked to make a selection.

        :param obj rd: the raw data (histogram) data object. This
          sets rd.instbank.

        '''
        if not os.path.exists(instfile): # no such file
            return {}
        fp = 0
        try:
            fp = open(instfile,'Ur')
            Iparm = {}
            for S in fp:
                Iparm[S[:12]] = S[12:-1]
        except IOError:
            print('Error reading file:'+str(instfile))
        if fp:        
            fp.close()

        ibanks = int(Iparm.get('INS   BANK  ','1').strip())
        hType = Iparm['INS   HTYPE '].strip()
        if ibanks == 1: # there is only one bank here, return it
            rd.instbank = 1
            return Iparm
        if 'PNT' in hType:
            rd.instbank = bank
        elif ibanks != databanks:
            # number of banks in data and prm file not not agree, need a
            # choice from a human here
            choices = []
            for i in range(1,1+ibanks):
                choices.append('Bank '+str(i))
            bank = rd.BlockSelector(
                choices, self,
                title='Select an instrument parameter bank for '+
                os.path.split(rd.powderentry[0])[1]+' BANK '+str(bank)+
                '\nOr use Cancel to select from the default parameter sets',
                header='Block Selector')
        if bank is None: return {}
        # pull out requested bank # bank from the data, and change the bank to 1
        IparmS = {}
        for key in Iparm:
            if key[4:6] == "  ":
                IparmS[key] = Iparm[key]
            elif int(key[4:6].strip()) == bank:
                IparmS[key[:4]+' 1'+key[6:]] = Iparm[key]
        rd.instbank = bank
        return IparmS
                        
    def GetPowderIparm(self,rd, prevIparm, lastIparmfile, lastdatafile):
        '''Open and read an instrument parameter file for a data file
        Returns the list of parameters used in the data tree

        :param obj rd: the raw data (histogram) data object.

        :param str prevIparm: not used

        :param str lastIparmfile: Name of last instrument parameter
          file that was read, or a empty string. 

        :param str lastdatafile: Name of last data file that was read.

        :returns: a list of two dicts, the first containing instrument parameters
          and the second used for future TOF datasets (timemaps?)

        '''
        def SetPowderInstParms(Iparm, rd):
            '''extracts values from instrument parameters in rd.instdict
            or in array Iparm.
            Create and return the contents of the instrument parameter tree entry.
            '''
            DataType = Iparm['INS   HTYPE '].strip()[:3]  # take 1st 3 chars
            # override inst values with values read from data file
            if rd.instdict.get('type'):
                DataType = rd.instdict.get('type')
            data = [DataType,]
            instname = Iparm.get('INS  1INAME ')
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
                    wave1 = G2IO.sfloat(s[:10])
                    wave2 = G2IO.sfloat(s[10:20])
                v = (wave1,wave2,
                     G2IO.sfloat(s[20:30]),G2IO.sfloat(s[55:65]),G2IO.sfloat(s[40:50])) #get lam1, lam2, zero, pola & ratio
                if not v[1]:
                    names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','SH/L','Azimuth'] 
                    v = (v[0],v[2],v[4])
                    codes = [0,0,0,0]
                else:
                    names = ['Type','Lam1','Lam2','Zero','I(L2)/I(L1)','Polariz.','U','V','W','X','Y','SH/L','Azimuth']
                    codes = [0,0,0,0,0,0]
                data.extend(v)
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
                return [G2IO.makeInstDict(names,data,codes),{}]
            elif 'T' in DataType:
                names = ['Type','2-theta','difC','difA','Zero','alpha','beta-0','beta-1',
                    'beta-q','sig-0','sig-1','sig-q','X','Y','Azimuth']
                codes = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                azm = 0.
                if 'INS  1DETAZM' in Iparm:
                    azm = float(Iparm['INS  1DETAZM'])
                s = Iparm['INS  1BNKPAR'].split()
                data.extend([G2IO.sfloat(s[1]),])               #2-theta for bank
                s = Iparm['INS  1 ICONS'].split()
                data.extend([G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),G2IO.sfloat(s[2])])    #difC, difA, Zero
                s = Iparm['INS  1PRCF1 '].split()
                pfType = int(s[0])
                s = Iparm['INS  1PRCF11'].split()
                if abs(pfType) == 1:
                    data.extend([G2IO.sfloat(s[1]),G2IO.sfloat(s[2]),G2IO.sfloat(s[3])])
                    s = Iparm['INS  1PRCF12'].split()
                    data.extend([0.0,0.0,G2IO.sfloat(s[1]),0.0,0.0,0.0,azm])
                elif abs(pfType) in [3,4,5]:
                    data.extend([G2IO.sfloat(s[0]),G2IO.sfloat(s[1]),G2IO.sfloat(s[2])])
                    if abs(pfType) == 4:
                        data.extend([0.0,0.0,G2IO.sfloat(s[3]),0.0,0.0,0.0,azm])
                    else:
                        s = Iparm['INS  1PRCF12'].split()
                        data.extend([0.0,0.0,G2IO.sfloat(s[0]),0.0,0.0,0.0,azm])                       
                Inst1 = G2IO.makeInstDict(names,data,codes)
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
                        for i in range(8):
                            s = Iparm['INS  1IECOR'+str(i+1)].split()
                            Icovar += [float(S) for S in s]
                        Inst2['Icoeff'] = Icoeff
                        Inst2['Iesd'] = Iesd
                        Inst2['Icovar'] = Icovar
                return [Inst1,Inst2]

        # stuff we might need from the reader
        filename = rd.powderentry[0]
        bank = rd.powderentry[2]
        numbanks = rd.numbanks
        # is there an instrument parameter file defined for the current data set?
        if rd.instparm or (lastdatafile == filename and lastIparmfile):
            if rd.instparm:
                instfile = os.path.join(os.path.split(filename)[0],
                                    rd.instparm)
            else:
                # for multiple reads of one data file, reuse the inst parm file
                instfile = lastIparmfile
            if os.path.exists(instfile):
                #print 'debug: try read',instfile
                instParmList = self.ReadPowderInstprm(instfile)
                if instParmList is not None:
                    rd.instfile = instfile
                    rd.instmsg = 'GSAS-II file '+instfile
                    return instParmList
                Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                if Iparm:
                    #print 'debug: success'
                    rd.instfile = instfile
                    rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                    return SetPowderInstParms(Iparm,rd)
            else:
                self.ErrorDialog('Open Error','Error opening instrument parameter file '
                    +str(instfile)+' requested by file '+ filename)
        # is there an instrument parameter file matching the current file
        # with extension .inst or .prm? If so read it
        basename = os.path.splitext(filename)[0]
        for ext in '.instprm','.prm','.inst','.ins':
            instfile = basename + ext
            instParmList = self.ReadPowderInstprm(instfile)
            if instParmList is not None:
                rd.instfile = instfile
                rd.instmsg = 'GSAS-II file '+instfile
                return instParmList
            Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
            if Iparm:
                #print 'debug: success'
                rd.instfile = instfile
                rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                return SetPowderInstParms(Iparm,rd)
            else:
                #print 'debug: open/read failed',instfile
                pass # fail silently

        # did we read the data file from a zip? If so, look there for a
        # instrument parameter file
        if self.zipfile:
            for ext in '.instprm','.prm','.inst','.ins':
                instfile = G2IO.ExtractFileFromZip(
                    self.zipfile,
                    selection=os.path.split(basename + ext)[1],
                    parent=self)
                if instfile is not None and instfile != self.zipfile:
                    print 'debug:',instfile,'created from ',self.zipfile
                    instParmList = self.ReadPowderInstprm(instfile)
                    if instParmList is not None:
                        rd.instfile = instfile
                        rd.instmsg = 'GSAS-II file '+instfile
                        return instParmList
                    Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                    if Iparm:
                        rd.instfile = instfile
                        rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                        return SetPowderInstParms(Iparm,rd)
                    else:
                        #print 'debug: open/read for',instfile,'from',self.zipfile,'failed'
                        pass # fail silently

        while True: # loop until we get a file that works or we get a cancel
            instfile = ''
            dlg = wx.FileDialog(
                self,
                'Choose inst. param file for "'
                +rd.idstring
                +'" (or Cancel for default)',
                '.', '',
                'GSAS iparm file (*.prm,*.inst,*.ins)|*.prm;*.inst;*.ins|'
                'GSAS-II iparm file (*.instprm)|*.instprm|'
                'All files (*.*)|*.*', 
                wx.OPEN|wx.CHANGE_DIR)
            if os.path.exists(lastIparmfile):
                dlg.SetFilename(lastIparmfile)
            if dlg.ShowModal() == wx.ID_OK:
                instfile = dlg.GetPath()
            dlg.Destroy()
            if not instfile: break
            instParmList = self.ReadPowderInstprm(instfile)
            if instParmList is not None:
                rd.instfile = instfile
                rd.instmsg = 'GSAS-II file '+instfile
                return instParmList
            Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
            if Iparm:
                #print 'debug: success with',instfile
                rd.instfile = instfile
                rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                return SetPowderInstParms(Iparm,rd)
            else:
                self.ErrorDialog('Read Error',
                                 'Error opening/reading file '+str(instfile))
        
        # still no success: offer user choice of defaults
        while True: # loop until we get a choice
            choices = []
            head = 'Select from default instrument parameters for '+rd.idstring

            for l in rd.defaultIparm_lbl:
                choices.append('Defaults for '+l)
            res = rd.BlockSelector(
                choices,
                ParentFrame=self,
                title=head,
                header='Select default inst parms',
                useCancel=False)
            if res is None: continue
            rd.instfile = ''
            rd.instmsg = 'default: '+rd.defaultIparm_lbl[res]
            return SetPowderInstParms(rd.defaultIparms[res],rd)

    def OnImportPowder(self,event):
        '''Called in response to an Import/Powder Data/... menu item
        to read a powder diffraction data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        Also reads an instrument parameter file for each dataset.
        '''
        # get a list of existing histograms
        PWDRlist = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportPowderReaderlist,'Powder Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        Iparm = None
        lastIparmfile = ''
        lastdatafile = ''
        newHistList = []
        for rd in rdlist:
            # get instrument parameters for each dataset
            Iparm1,Iparm2 = self.GetPowderIparm(rd, Iparm, lastIparmfile, lastdatafile)
            if rd.repeat_instparm: 
                lastIparmfile = rd.instfile
            lastdatafile = rd.powderentry[0]
            HistName = rd.idstring
            HistName = 'PWDR '+HistName
            # make new phase names unique
            HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)
            print 'Read powder data '+str(HistName)+ \
                ' from file '+str(self.lastimport) + \
                ' with parameters from '+str(rd.instmsg)
            # data are read, now store them in the tree
            Id = self.PatternTree.AppendItem(parent=self.root,text=HistName)
            if 'T' in Iparm1['Type'][0]:
                if not rd.clockWd and rd.GSAS:
                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
                cw = np.diff(rd.powderdata[0])
                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
                rd.powderdata[1] = rd.powderdata[1][:-1]/cw
                rd.powderdata[2] = rd.powderdata[2][:-1]*cw**2  #1/var=w at this point
                if 'Itype' in Iparm2:
                    Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
                    Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
                    rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
                    YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
                    rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
                    var = 1./rd.powderdata[2][Ibeg:Ifin]
                    var += WYI*rd.powderdata[1]**2
                    var /= YI**2
                    rd.powderdata[2] = 1./var
                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])                                        
                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])                                        
                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])                                        
            Tmin = min(rd.powderdata[0])
            Tmax = max(rd.powderdata[0])
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxint),
                }
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.PatternTree.SetItemPyData(Id,[valuesdict,rd.powderdata])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Comments'),
                rd.comments)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Limits'),
                [(Tmin,Tmax),[Tmin,Tmax]])
            self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Background'),
                [['chebyschev',True,3,1.0,0.0,0.0],
                 {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Peak List')
                ,[])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Index Peak List'),
                [])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Unit Cells List'),
                [])
            self.PatternTree.SetItemPyData(
                self.PatternTree.AppendItem(Id,text='Reflection Lists'),
                {})
            self.PatternTree.Expand(Id)
            self.PatternTree.SelectItem(Id)
            newHistList.append(HistName)
            
        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = usedHistograms.keys() # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\npowder dataset(s) to:'
        for Name in newHistList:
            header += '\n  '+str(Name)

        notOK = True
        result = G2gd.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        item, cookie = self.PatternTree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            phaseName = self.PatternTree.GetItemText(item)
            data = self.PatternTree.GetItemPyData(item)
            item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for histoName in newHistList:
                UseList[histoName] = {
                    'Histogram':histoName,'Show':False,
                    'Scale':[1.0,False],'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{}],
                    'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                            [1.,1.,1.,0.,0.,0.],6*[False,]],
                    'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                                NShkl*[0.01,],NShkl*[False,]],
                    'HStrain':[NDij*[0.0,],NDij*[False,]],                          
                    'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}
                Id = G2gd.GetPatternTreeItemId(self,self.root,histoName)
                refList = self.PatternTree.GetItemPyData(
                    G2gd.GetPatternTreeItemId(self,Id,'Reflection Lists'))
                refList[generalData['Name']] = []
        return # success

    def OnDummyPowder(self,event):
        '''Called in response to Import/Powder Data/Simulate menu item
        to create a Dummy powder diffraction data set. 

        Reads an instrument parameter file and then gets input from the user
        '''
        rd = G2IO.ImportPowderData(
            extensionlist=tuple(),
            strictExtension=False,
            formatName = 'Simulate dataset',
            longFormatName = 'Compute a simulated pattern')
        rd.powderentry[0] = '' # no filename
        # #self.powderentry[1] = pos # bank offset (N/A here)
        rd.powderentry[2] = 1 # only one bank
        rd.comments.append('This is a dummy dataset for powder pattern simulation')
        self.CheckNotebook()
        Iparm = None
        lastIparmfile = ''
        lastdatafile = ''
        self.zipfile = None
        # get instrument parameters for 
        Iparm1,Iparm2 = self.GetPowderIparm(rd, Iparm, lastIparmfile, lastdatafile)
        if 'T' in Iparm1['Type'][0]:
            print('TOF simulation not supported yet')
            return False
        else:
            # need to get name, 2theta start, end, step
            rd.idstring = ' CW'
            if 'X' in Iparm1['Type'][0]:
                rd.idstring = 'CW x-ray simulation'
            else:
                rd.idstring = 'CW neutron simulation'
            # base initial range on wavelength
            wave = Iparm1.get('Lam')
            if wave:
                wave = wave[0]
            else:
                wave = Iparm1.get('Lam1')
                if wave:
                    wave = wave[0]
        N = 0
        while (N < 3): # insist on a dataset with a few points
            names = ('dataset name', 'start angle', 'end angle', 'step size')
            if not wave or wave < 1.0:
                inp = [rd.idstring, 10.,40.,0.005] # see names for what's what
            else:
                inp = [rd.idstring, 10.,80.,0.01] # see names for what's what
            dlg = G2gd.ScrolledMultiEditor(
                self,[inp] * len(inp),range(len(inp)),names,
                header='Enter simulation name and range',
                minvals=(None,0.001,0.001,0.0001),
                maxvals=(None,180.,180.,.1),
                sizevals=((225,-1),)
                )
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                if inp[1] > inp[2]:
                    end,start,step = inp[1:]
                else:                
                    start,end,step = inp[1:]
                step = abs(step)
            else:
                return False
            N = int((end-start)/step)+1
            x = np.linspace(start,end,N,True)
            N = len(x)
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
        # data are read, now store them in the tree
        Id = self.PatternTree.AppendItem(parent=self.root,
                                         text='PWDR '+inp[0])
        valuesdict = {
            'wtFactor':1.0,
            'Dummy':True,
            'ranId':ran.randint(0,sys.maxint),
            }
        self.PatternTree.SetItemPyData(Id,[valuesdict,rd.powderdata])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Comments'),
            rd.comments)
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        self.PatternId = G2gd.GetPatternTreeItemId(self,Id,'Limits')
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Background'),
            [['chebyschev',True,3,1.0,0.0,0.0],
             {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Instrument Parameters'),
            [Iparm1,Iparm2])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Sample Parameters'),
            rd.Sample)
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Peak List')
            ,[])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Index Peak List'),
            [])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Unit Cells List'),
            [])
        self.PatternTree.SetItemPyData(
            self.PatternTree.AppendItem(Id,text='Reflection Lists'),
            {})
        self.PatternTree.Expand(Id)
        self.PatternTree.SelectItem(Id)
        print('Added simulation powder data '+str(rd.idstring)+ 
              ' with parameters from '+str(rd.instmsg))
        return # success

    def _init_Exports(self,menu):
        '''Find exporter routines and add them into menus
        '''
        # set up the top-level menus
        projectmenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Entire project as',
            projectmenu, help='Export entire project')

        phasemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Phase as',
            phasemenu, help='Export phase or sometimes phases')

        powdermenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Powder data as',
            powdermenu, help='Export powder diffraction histogram(s)')

        singlemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Single crystal data as',
            singlemenu, help='Export single crystal histogram(s)')

        imagemenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Images as',
            imagemenu, help='Export powder image(s)')

        mapmenu = wx.Menu()
        item = menu.AppendMenu(
            wx.ID_ANY, 'Maps as',
            mapmenu, help='Export density map(s)')

        # pdfmenu = wx.Menu()
        # item = menu.AppendMenu(
        #     wx.ID_ANY, 'PDFs as',
        #     pdfmenu, help='Export pair distribution function(s)')

        # find all the exporter files
        pathlist = sys.path
        filelist = []
        for path in pathlist:
            for filename in glob.iglob(os.path.join(
                path,
                "G2export*.py")):
                filelist.append(filename)    
        filelist = sorted(list(set(filelist))) # remove duplicates
        exporterlist = []
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
                    if inspect.isclass(clss[1]):
                        # check if we have the required methods
                        for m in 'Exporter','loadParmDict':
                            if not hasattr(clss[1],m): break
                            if not callable(getattr(clss[1],m)): break
                        else:
                            exporter = clss[1](self) # create an export instance
                            exporterlist.append(exporter)
            except AttributeError:
                print 'Import_'+errprefix+': Attribute Error'+str(filename)
                pass
            except ImportError:
                print 'Import_'+errprefix+': Error importing file'+str(filename)
                pass
            if fp: fp.close()
        # Add submenu item(s) for each Exporter by its self-declared type (can be more than one)
        for obj in exporterlist:
            #print 'exporter',obj
            for typ in obj.exporttype:
                if typ == "project":
                    submenu = projectmenu
                elif typ == "phase":
                    submenu = phasemenu
                elif typ == "powder":
                    submenu = powdermenu
                elif typ == "single":
                    submenu = singlemenu
                elif typ == "image":
                    submenu = imagemenu
                elif typ == "map":
                    submenu = mapmenu
                # elif typ == "pdf":
                #     submenu = pdfmenu
                else:
                    print("Error, unknown type in "+str(obj))
                    break
                item = submenu.Append(
                    wx.ID_ANY,
                    help=obj.longFormatName,
                    kind=wx.ITEM_NORMAL,
                    text=obj.formatName)
                self.Bind(wx.EVT_MENU, obj.Exporter, id=item.GetId())
                self.ExportLookup[item.GetId()] = typ # lookup table for submenu item
        #code to debug an Exporter. much is hard-coded below, but code is reloaded before
        # each use allowing faster development
        # def DebugExport(event):
        #     print 'start reload'
        #     reload(G2IO)
        #     import G2export_ASCIImap as dev
        #     reload(dev)
        #     dev.ExportMapASCII(self).Exporter(event)
        # item = menu.Append(
        #     wx.ID_ANY,kind=wx.ITEM_NORMAL,
        #     help="debug exporter",text="test Export")
        # self.Bind(wx.EVT_MENU, DebugExport, id=item.GetId())
        # #self.ExportLookup[item.GetId()] = 'image'
        # self.ExportLookup[item.GetId()] = 'map'

    def _Add_ExportMenuItems(self,parent):
        item = parent.Append(
            help='Select PWDR item to enable',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export Powder Patterns...')
        self.ExportPattern.append(item)
        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPatterns, id=item.GetId())

        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export All Peak Lists...')
        self.ExportPeakList.append(item)
        item.Enable(True)
        self.Bind(wx.EVT_MENU, self.OnExportPeakList, id=item.GetId())

        item = parent.Append(
            help='',id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export HKLs...')
        self.ExportHKL.append(item)
        self.Bind(wx.EVT_MENU, self.OnExportHKL, id=item.GetId())

        item = parent.Append(
            help='Select PDF item to enable',
            id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,
            text='Export PDF...')
        self.ExportPDF.append(item)
        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPDF, id=item.GetId())

    def FillMainMenu(self,menubar):
        '''Define contents of the main GSAS-II menu for the (main) data tree window
        in the mac, used also for the data item windows as well.
        '''
        File = wx.Menu(title='')
        menubar.Append(menu=File, title='&File')
        self._Add_FileMenuItems(File)
        Data = wx.Menu(title='')
        menubar.Append(menu=Data, title='Data')
        self._Add_DataMenuItems(Data)
        Calculate = wx.Menu(title='')        
        menubar.Append(menu=Calculate, title='&Calculate')
        self._Add_CalculateMenuItems(Calculate)
        Import = wx.Menu(title='')        
        menubar.Append(menu=Import, title='Import')
        self._Add_ImportMenu_Phase(Import)
        self._Add_ImportMenu_powder(Import)
        self._Add_ImportMenu_Sfact(Import)
        #======================================================================
        # Code to help develop/debug an importer, much is hard-coded below
        # but module is reloaded before each use, allowing faster testing
        # def DebugImport(event):
        #     print 'start reload'
        #     import G2phase_ISO as dev
        #     reload(dev)
        #     rd = dev.ISODISTORTPhaseReader()
        #     self.ImportMenuId[event.GetId()] = rd
        #     self.OnImportPhase(event)
            # or ----------------------------------------------------------------------
            #self.OnImportGeneric(rd,[],'test of ISODISTORTPhaseReader')
            # special debug code
            # or ----------------------------------------------------------------------
            # filename = '/Users/toby/projects/branton/subgroup_cif.txt'
            # fp = open(filename,'Ur')
            # if not rd.ContentsValidator(fp):
            #     print 'not validated'
            #     # make a list of used phase ranId's
            # phaseRIdList = []
            # sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
            # if sub:
            #     item, cookie = self.PatternTree.GetFirstChild(sub)
            #     while item:
            #         phaseName = self.PatternTree.GetItemText(item)
            #         ranId = self.PatternTree.GetItemPyData(item).get('ranId')
            #         if ranId: phaseRIdList.append(ranId)
            #         item, cookie = self.PatternTree.GetNextChild(sub, cookie)
            # if rd.Reader(filename,fp,usedRanIdList=phaseRIdList):
            #     print 'read OK'
        # item = Import.Append(
        #     wx.ID_ANY,kind=wx.ITEM_NORMAL,
        #     help="debug importer",text="test importer")
        # self.Bind(wx.EVT_MENU, DebugImport, id=item.GetId())
        #======================================================================
        self.ExportMenu = wx.Menu(title='')
        menubar.Append(menu=self.ExportMenu, title='Export')
        self._init_Exports(self.ExportMenu)
        self._Add_ExportMenuItems(self.ExportMenu)
        HelpMenu=G2gd.MyHelp(self,helpType='Data tree',
            morehelpitems=[('&Tutorials','Tutorials')])
        menubar.Append(menu=HelpMenu,title='&Help')

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='GSASII', parent=parent,
            size=wx.Size(400, 250),style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II data tree')
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]))
        self._init_Imports()
        #initialize Menu item objects (these contain lists of menu items that are enabled or disabled)
        self.MakePDF = []
        self.Refine = []
        self.SeqRefine = []
        self.ExportPattern = []
        self.ExportPeakList = []
        self.ExportHKL = []
        self.ExportPDF = []
        self.ExportPhase = []
        self.ExportCIF = []
        #
        self.GSASIIMenu = wx.MenuBar()
        self.FillMainMenu(self.GSASIIMenu)
        self.SetMenuBar(self.GSASIIMenu)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.CreateStatusBar()
        self.mainPanel = wx.Panel(self,-1)
        
        wxID_PATTERNTREE = wx.NewId()
        self.PatternTree = wx.TreeCtrl(id=wxID_PATTERNTREE,
            parent=self.mainPanel, pos=wx.Point(0, 0),style=wx.TR_DEFAULT_STYLE )
        self.PatternTree.Bind(wx.EVT_TREE_SEL_CHANGED,
            self.OnPatternTreeSelChanged, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
            self.OnPatternTreeItemCollapsed, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
            self.OnPatternTreeItemExpanded, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_DELETE_ITEM,
            self.OnPatternTreeItemDelete, id=wxID_PATTERNTREE)
        self.PatternTree.Bind(wx.EVT_TREE_KEY_DOWN,
            self.OnPatternTreeKeyDown, id=wxID_PATTERNTREE)
        self.root = self.PatternTree.AddRoot('Loaded Data: ')
        
        plotFrame = wx.Frame(None,-1,'GSASII Plots',size=wx.Size(700,600), \
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(plotFrame)
        plotFrame.Show()
        
        self.dataDisplay = None
        
    def __init__(self, parent):
        self.ExportLookup = {}
        self._init_ctrls(parent)
        self.Image = wx.Image(
            os.path.join(GSASIIpath.path2GSAS2,'gsas2.ico'),
            wx.BITMAP_TYPE_ICO)
        if "wxMSW" in wx.PlatformInfo:
            img = self.Image.Scale(16, 16).ConvertToBitmap()
        elif "wxGTK" in wx.PlatformInfo:
            img = self.Image.Scale(22, 22).ConvertToBitmap()
        else:
            img = self.Image.ConvertToBitmap()
        self.SetIcon(wx.IconFromBitmap(img))
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        # various defaults
        self.oldFocus = None
        self.GSASprojectfile = ''
        self.dirname = os.path.expanduser('~')       #start in the users home directory by default; may be meaningless
        self.undofile = ''
        self.TreeItemDelete = False
        self.Offset = [0.0,0.0]
        self.delOffset = .02
        self.refOffset = -100.0
        self.refDelt = .01
        self.Weight = False
        self.IparmName = ''  # to be removed when SelectPowderData & GetInstrumentFile is
        self.IfPlot = False
        self.PatternId = 0
        self.PickId = 0
        self.PeakTable = []
        self.LimitsTable = []
        self.HKL = []
        self.Lines = []
        self.itemPicked = None
        self.dataFrame = None
        self.Interpolate = 'nearest'
        self.ContourColor = 'Paired'
        self.VcovColor = 'RdYlGn'
        self.RamaColor = 'Blues'
        self.Projection = 'equal area'
        self.logPlot = False
        self.qPlot = False
        self.Contour = False
        self.Legend = False
        self.SinglePlot = False
        self.SubBack = False
        self.plotView = 0
        self.Image = 0
        self.oldImagefile = ''
        self.ImageZ = []
        self.Integrate = 0
        self.imageDefault = {}
        self.Sngl = 0
        self.ifGetRing = False
        self.MaskKey = ''           #trigger for making image masks
        arg = sys.argv
        if len(arg) > 1:
            self.GSASprojectfile = os.path.splitext(arg[1])[0]+'.gpx'
            self.dirname = os.path.dirname(arg[1])
            if self.dirname: os.chdir(self.dirname)
            try:
                G2IO.ProjFileOpen(self)
                self.PatternTree.Expand(self.root)
                for item in self.Refine: item.Enable(True)
                for item in self.SeqRefine: item.Enable(True)
            except:
                print 'Error opening file',arg[1]

    def OnSize(self,event):
        'Called when the main window is resized. Not sure why'
        w,h = self.GetClientSizeTuple()
        self.mainPanel.SetSize(wx.Size(w,h))
        self.PatternTree.SetSize(wx.Size(w,h))
                        
    def OnPatternTreeSelChanged(self, event):
        '''Called when a data tree item is selected'''
        if self.TreeItemDelete:
            self.TreeItemDelete = False
        else:
            pltNum = self.G2plotNB.nb.GetSelection()
            if pltNum >= 0:                         #to avoid the startup with no plot!
                pltPage = self.G2plotNB.nb.GetPage(pltNum)
                pltPlot = pltPage.figure
            item = event.GetItem()
            G2gd.MovePatternTreeToGrid(self,item)
            if self.oldFocus:
                self.oldFocus.SetFocus()
        
    def OnPatternTreeItemCollapsed(self, event):
        'Called when a tree item is collapsed'
        event.Skip()

    def OnPatternTreeItemExpanded(self, event):
        'Called when a tree item is expanded'
        event.Skip()
        
    def OnPatternTreeItemDelete(self, event):
        'Called when a tree item is deleted -- not sure what this does'
        self.TreeItemDelete = True

    def OnPatternTreeItemActivated(self, event):
        'Called when a tree item is activated'
        event.Skip()
        
    def OnPatternTreeKeyDown(self,event):
        'Not sure what this does'
        key = event.GetKeyCode()
        item = self.PickId
        if type(item) is int: return # is this the toplevel in tree?
        if key == wx.WXK_UP:
            self.oldFocus = wx.Window.FindFocus()
            self.PatternTree.GetPrevSibling(item)
        elif key == wx.WXK_DOWN:
            self.oldFocus = wx.Window.FindFocus()
            self.PatternTree.GetNextSibling(item)
                
    def OnReadPowderPeaks(self,event):
        'Bound to menu Data/Read Powder Peaks -- still needed?'
        Cuka = 1.54052
        self.CheckNotebook()
        dlg = wx.FileDialog(self, 'Choose file with peak list', '.', '', 
            'peak files (*.txt)|*.txt|All files (*.*)|*.*',wx.OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.HKL = []
                self.powderfile = dlg.GetPath()
                comments,peaks = G2IO.GetPowderPeaks(self.powderfile)
                Id = self.PatternTree.AppendItem(parent=self.root,text='PKS '+os.path.basename(self.powderfile))
                data = ['PKS',Cuka,0.0]
                names = ['Type','Lam','Zero'] 
                codes = [0,0,0]
                inst = [G2IO.makeInstDict(names,data,codes),{}]
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),inst)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),comments)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),peaks)
                self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                self.PatternTree.Expand(Id)
                self.PatternTree.SelectItem(Id)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
        finally:
            dlg.Destroy()
                        
    def OnImageRead(self,event):
        'Called to read in an image in any known format'
        self.CheckNotebook()
        dlg = wx.FileDialog(
            self, 'Choose image files', '.', '',
            'Any image file (*.edf;*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img;*.G2img)|'
            '*.edf;*.tif;*.tiff;*.mar*;*.avg;*.sum;*.img;*.G2img;*.zip|'
            'European detector file (*.edf)|*.edf|'
            'Any detector tif (*.tif;*.tiff)|*.tif;*.tiff|'
            'MAR file (*.mar*)|*.mar*|'
            'GE Image (*.avg;*.sum)|*.avg;*.sum|'
            'ADSC Image (*.img)|*.img|'
            'GSAS-II Image (*.G2img)|*.G2img|'
            'Zip archive (*.zip)|*.zip|'
            'All files (*.*)|*.*',
            wx.OPEN | wx.MULTIPLE|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                imagefiles = dlg.GetPaths()
                imagefiles.sort()
                for imagefile in imagefiles:
                    # if a zip file, open and extract
                    if os.path.splitext(imagefile)[1].lower() == '.zip':
                        extractedfile = G2IO.ExtractFileFromZip(imagefile,parent=self)
                        if extractedfile is not None and extractedfile != imagefile:
                            imagefile = extractedfile
                    Comments,Data,Npix,Image = G2IO.GetImageData(self,imagefile)
                    if Comments:
                        Id = self.PatternTree.AppendItem(parent=self.root,text='IMG '+os.path.basename(imagefile))
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        Imax = np.amax(Image)
                        Imin = max(0.0,np.amin(Image))          #force positive
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                            Data['showLines'] = True
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['cutoff'] = 10
                            Data['pixLimit'] = 20
                            Data['edgemin'] = 100000000
                            Data['calibdmin'] = 0.5
                            Data['calibskip'] = 0
                            Data['ellipses'] = []
                            Data['calibrant'] = ''
                            Data['GonioAngles'] = [0.,0.,0.]
                            Data['DetDepth'] = 0.
                            Data['DetDepthRef'] = False
                        else:
                            Data['type'] = 'PWDR'
                            Data['color'] = 'Paired'
                            Data['tilt'] = 0.0
                            Data['rotation'] = 0.0
                            Data['showLines'] = False
                            Data['ring'] = []
                            Data['rings'] = []
                            Data['cutoff'] = 10
                            Data['pixLimit'] = 20
                            Data['calibdmin'] = 0.5
                            Data['calibskip'] = 0
                            Data['edgemin'] = 100000000
                            Data['ellipses'] = []
                            Data['GonioAngles'] = [0.,0.,0.]
                            Data['DetDepth'] = 0.
                            Data['DetDepthRef'] = False
                            Data['calibrant'] = ''
                            Data['IOtth'] = [2.0,5.0]
                            Data['LRazimuth'] = [135,225]
                            Data['azmthOff'] = 0.0
                            Data['outChannels'] = 2500
                            Data['outAzimuths'] = 1
                            Data['centerAzm'] = False
                            Data['fullIntegrate'] = False
                            Data['setRings'] = False
                            Data['background image'] = ['',1.0]                            
                        Data['setDefault'] = False
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Stress/Strain'),
                            {'Type':'True','d-zero':[],'Sample phi':0.0,'Sample z':0.0,'strain':np.zeros((3,3))})
                        self.PatternTree.SetItemPyData(Id,[Npix,imagefile])
                        self.PickId = Id
                        self.Image = Id
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
                self.PatternTree.SelectItem(G2gd.GetPatternTreeItemId(self,Id,'Image Controls'))             #show last one
        finally:
            path = dlg.GetDirectory()           # to get Mac/Linux to change directory!
            os.chdir(path)
            dlg.Destroy()

    def CheckNotebook(self):
        '''Make sure the data tree has the minimally expected controls.
        (BHT) correct?
        '''
        if not G2gd.GetPatternTreeItemId(self,self.root,'Notebook'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Notebook')
            self.PatternTree.SetItemPyData(sub,[''])
        if not G2gd.GetPatternTreeItemId(self,self.root,'Controls'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Controls')
            self.PatternTree.SetItemPyData(sub,copy.copy(G2obj.DefaultControls))
        if not G2gd.GetPatternTreeItemId(self,self.root,'Covariance'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Covariance')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Constraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Constraints')
            self.PatternTree.SetItemPyData(sub,{'Hist':[],'HAP':[],'Phase':[]})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Restraints'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Restraints')
            self.PatternTree.SetItemPyData(sub,{})
        if not G2gd.GetPatternTreeItemId(self,self.root,'Rigid bodies'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Rigid bodies')
            self.PatternTree.SetItemPyData(sub,{'Vector':{'AtInfo':{}},
                'Residue':{'AtInfo':{}},'RBIds':{'Vector':[],'Residue':[]}})
                
    class CopyDialog(wx.Dialog):
        '''Creates a dialog for copying control settings between
        data tree items'''
        def __init__(self,parent,title,text,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            ncols = len(data)/40+1
            dataGridSizer = wx.FlexGridSizer(rows=len(data),cols=ncols,hgap=2,vgap=2)
            for id,item in enumerate(self.data):
                ckbox = wx.CheckBox(panel,id,item[1])
                ckbox.Bind(wx.EVT_CHECKBOX,self.OnCopyChange)                    
                dataGridSizer.Add(ckbox,0,wx.LEFT,10)
            mainSizer.Add(dataGridSizer,0,wx.EXPAND)
            OkBtn = wx.Button(panel,-1,"Ok")
            OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((20,20),1)
            
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            panel.SetSizer(mainSizer)
            panel.Fit()
            self.Fit()
        
        def OnCopyChange(self,event):
            id = event.GetId()
            self.data[id][0] = self.FindWindowById(id).GetValue()        
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            
        def GetData(self):
            return self.data
        
    class SumDialog(wx.Dialog):
        'Allows user to supply scale factor(s) when summing data'
        def __init__(self,parent,title,text,dataType,data):
            wx.Dialog.__init__(self,parent,-1,title, 
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
            self.data = data
            panel = wx.Panel(self)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            topLabl = wx.StaticText(panel,-1,text)
            mainSizer.Add((10,10),1)
            mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
            mainSizer.Add((10,10),1)
            dataGridSizer = wx.FlexGridSizer(rows=len(data),cols=2,hgap=2,vgap=2)
            for id,item in enumerate(self.data[:-1]):
                name = wx.TextCtrl(panel,-1,item[1],size=wx.Size(200,20))
                name.SetEditable(False)
                scale = wx.TextCtrl(panel,id,'%.3f'%(item[0]),style=wx.TE_PROCESS_ENTER)
                scale.Bind(wx.EVT_TEXT_ENTER,self.OnScaleChange)
                scale.Bind(wx.EVT_KILL_FOCUS,self.OnScaleChange)
                dataGridSizer.Add(scale,0,wx.LEFT,10)
                dataGridSizer.Add(name,0,wx.RIGHT,10)
            if dataType:
                dataGridSizer.Add(wx.StaticText(panel,-1,'Sum result name: '+dataType),0, \
                    wx.LEFT|wx.TOP|wx.ALIGN_CENTER_VERTICAL,10)
                self.name = wx.TextCtrl(panel,-1,self.data[-1],size=wx.Size(200,20),style=wx.TE_PROCESS_ENTER)
                self.name.Bind(wx.EVT_TEXT_ENTER,self.OnNameChange)
                self.name.Bind(wx.EVT_KILL_FOCUS,self.OnNameChange)
                dataGridSizer.Add(self.name,0,wx.RIGHT|wx.TOP,10)
            mainSizer.Add(dataGridSizer,0,wx.EXPAND)
            OkBtn = wx.Button(panel,-1,"Ok")
            OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.BoxSizer(wx.HORIZONTAL)
            btnSizer.Add((20,20),1)
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((20,20),1)
            
            mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
            panel.SetSizer(mainSizer)
            panel.Fit()
            self.Fit()

        def OnScaleChange(self,event):
            id = event.GetId()
            value = self.FindWindowById(id).GetValue()
            try:
                self.data[id][0] = float(value)
                self.FindWindowById(id).SetValue('%.3f'%(self.data[id][0]))
            except ValueError:
                if value and '-' not in value[0]:
                    print 'bad input - numbers only'
                    self.FindWindowById(id).SetValue('0.000')
            
        def OnNameChange(self,event):
            self.data[-1] = self.name.GetValue() 
            
        def OnOk(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            
        def GetData(self):
            return self.data
                        
    def OnPwdrSum(self,event):
        'Sum together powder data(?)'
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = None
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetItemPyData(item)[1])    # (x,y,w,yc,yb,yd)
                    if not Inst:
                        Inst = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item, 'Instrument Parameters'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "PWDR" pattern')
                return
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each pattern in summation','PWDR',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    lenX = 0
                    Xminmax = [0,0]
                    Xsum = []
                    Ysum = []
                    Vsum = []
                    result = dlg.GetData()
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            x,y,w,yc,yb,yd = data   #numpy arrays!
                            v = 1./w
                            if lenX:
                                if lenX != len(x):
                                    self.ErrorDialog('Data length error','Data to be summed must have same number of points'+ \
                                        '\nExpected:'+str(lenX)+ \
                                        '\nFound:   '+str(len(x))+'\nfor '+name)
                                    return
                            else:
                                lenX = len(x)
                            if Xminmax[1]:
                                if Xminmax != [x[0],x[-1]]:
                                    self.ErrorDialog('Data range error','Data to be summed must span same range'+ \
                                        '\nExpected:'+str(Xminmax[0])+' '+str(Xminmax[1])+ \
                                        '\nFound:   '+str(x[0])+' '+str(x[-1])+'\nfor '+name)
                                    return
                                else:
                                    for j,yi in enumerate(y):
                                         Ysum[j] += scale*yi
                                         Vsum[j] += abs(scale)*v[j]
                            else:
                                Xminmax = [x[0],x[-1]]
                                YCsum = YBsum = YDsum = [0.0 for i in range(lenX)]
                                for j,yi in enumerate(y):
                                    Xsum.append(x[j])
                                    Ysum.append(scale*yi)
                                    Vsum.append(abs(scale*v[j]))
                    Wsum = 1./np.array(Vsum)
                    outname = 'PWDR '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                                self.PatternTree.Delete(Id)
                        finally:
                            dlg2.Destroy()
                    Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        Sample = G2pdG.SetDefaultSample()
                        valuesdict = {
                            'wtFactor':1.0,
                            'Dummy':False,
                            'ranId':ran.randint(0,sys.maxint),
                            }
                        self.PatternTree.SetItemPyData(Id,[valuesdict,[np.array(Xsum),np.array(Ysum),np.array(Wsum),
                            np.array(YCsum),np.array(YBsum),np.array(YDsum)]])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)                    
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Background'),[['chebyschev',True,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Instrument Parameters'),Inst)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Sample Parameters'),Sample)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Index Peak List'),[])
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Unit Cells List'),[])             
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Reflection Lists'),{})             
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
            finally:
                dlg.Destroy()

    def OnImageSum(self,event):
        'Sum together image data(?)'
        TextList = []
        DataList = []
        SumList = []
        Names = []
        Inst = []
        SumItemList = []
        Comments = ['Sum equals: \n']
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append([0.0,name])
                    DataList.append(self.PatternTree.GetItemPyData(item))        #Size,Image
                    Data = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,item,'Image Controls'))
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "IMG" pattern')
                return
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data','Enter scale for each image in summation','IMG',TextList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    imSize = 0
                    result = dlg.GetData()
                    First = True
                    Found = False
                    for i,item in enumerate(result[:-1]):
                        scale,name = item
                        data = DataList[i]
                        if scale:
                            Found = True                                
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            Npix,imagefile = data
                            image = G2IO.GetImageData(self,imagefile,imageOnly=True)
                            if First:
                                newImage = np.zeros_like(image)
                                First = False
                            if imSize:
                                if imSize != Npix:
                                    self.ErrorDialog('Image size error','Images to be summed must be same size'+ \
                                        '\nExpected:'+str(imSize)+ \
                                        '\nFound:   '+str(Npix)+'\nfor '+name)
                                    return
                                newImage = newImage+scale*image
                            else:
                                imSize = Npix
                                newImage = newImage+scale*image
                            del(image)
                    if not Found:
                        self.ErrorDialog('Image sum error','No nonzero image multipliers found')
                        return
                        
                    newImage = np.asfarray(newImage,dtype=np.float32)                        
                    outname = 'IMG '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = G2gd.GetPatternTreeItemId(self,self.root,name)
                        finally:
                            dlg2.Destroy()
                    else:
                        Id = self.PatternTree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        dlg = wx.FileDialog(self, 'Choose sum image filename', '.', '', 
                            'G2img files (*.G2img)|*.G2img', 
                            wx.SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
                        if dlg.ShowModal() == wx.ID_OK:
                            newimagefile = dlg.GetPath()
                            newimagefile = G2IO.FileDlgFixExt(dlg,newimagefile)
                            G2IO.PutG2Image(newimagefile,Comments,Data,Npix,newImage)
                            Imax = np.amax(newImage)
                            Imin = np.amin(newImage)
                            newImage = []
                            self.PatternTree.SetItemPyData(Id,[imSize,newimagefile])
                            self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Comments'),Comments)
                        del(newImage)
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['cutoff'] = 10
                        Data['pixLimit'] = 20
                        Data['ellipses'] = []
                        Data['calibrant'] = ''
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Image Controls'),Data)                                            
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]]}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Masks'),Masks)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='Stress/Strain'),{})
                        self.PatternTree.SelectItem(Id)
                        self.PatternTree.Expand(Id)
                        self.PickId = G2gd.GetPatternTreeItemId(self,self.root,outname)
                        self.Image = self.PickId
            finally:
                dlg.Destroy()
                      
    def OnAddPhase(self,event):
        'Add a new, empty phase to the tree. Called by Data/Add Phase menu'
        if not G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = self.PatternTree.AppendItem(parent=self.root,text='Phases')
        else:
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            PhaseName = dlg.GetValue()
        dlg.Destroy()
        sub = self.PatternTree.AppendItem(parent=sub,text=PhaseName)
        E,SGData = G2spc.SpcGroup('P 1')
        self.PatternTree.SetItemPyData(sub,G2IO.SetNewPhase(Name=PhaseName,SGData=SGData))
        
    def OnDeletePhase(self,event):
        'Delete a phase from the tree. Called by Data/Delete Phase menu'
        #Hmm, also need to delete this phase from Reflection Lists for each PWDR histogram
        if self.dataFrame:
            self.dataFrame.Clear() 
        TextList = []
        DelList = []
        DelItemList = []
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            return
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                TextList.append(self.PatternTree.GetItemText(item))
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)                
            dlg = wx.MultiChoiceDialog(self, 'Which phase to delete?', 'Delete phase', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append([i,TextList[i]])
                    item, cookie = self.PatternTree.GetFirstChild(sub)
                    i = 0
                    while item:
                        if [i,self.PatternTree.GetItemText(item)] in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(sub, cookie)
                        i += 1
                    for item in DelItemList:
                        name = self.PatternTree.GetItemText(item)
                        self.PatternTree.Delete(item)
                        self.G2plotNB.Delete(name)
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        name = self.PatternTree.GetItemText(item)
                        if 'PWDR' in name:
                            Id = G2gd.GetPatternTreeItemId(self,item, 'Reflection Lists')
                            refList = self.PatternTree.GetItemPyData(Id)
                            for i,item in DelList:
                                del(refList[item])
                            self.PatternTree.SetItemPyData(Id,refList)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            finally:
                dlg.Destroy()
                
    def OnRenameData(self,event):
        'Renames an existing phase. Called by Data/Rename Phase menu'
        name = self.PatternTree.GetItemText(self.PickId)      
        if 'PWDR' in name or 'HKLF' in name or 'IMG' in name:
            dataType = name[:name.index(' ')+1]                 #includes the ' '
            dlg = wx.TextEntryDialog(self,'Data name: '+dataType,'Change data name',
                defaultValue=name[name.index(' ')+1:])
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    self.PatternTree.SetItemText(self.PickId,dataType+dlg.GetValue())
            finally:
                dlg.Destroy()
        
    def GetFileList(self,fileType,skip=None):        #potentially useful?
        'Appears unused. Note routine of same name in GSASIIpwdGUI'
        fileList = []
        Source = ''
        id, cookie = self.PatternTree.GetFirstChild(self.root)
        while id:
            name = self.PatternTree.GetItemText(id)
            if fileType in name:
                if id == skip:
                    Source = name
                else:
                    fileList.append([False,name,id])
            id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
            
    def OnDataDelete(self, event):
        '''Delete one or more histograms from data tree. Called by the
        Data/DeleteData menu
        '''
        TextList = ['All Data']
        DelList = []
        DelItemList = []
        ifPWDR = False
        ifIMG = False
        ifHKLF = False
        ifPDF = False
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name not in ['Notebook','Controls','Covariance','Constraints',
                    'Restraints','Phases','Rigid bodies']:
                    if 'PWDR' in name: ifPWDR = True
                    if 'IMG' in name: ifIMG = True
                    if 'HKLF' in name: ifHKLF = True
                    if 'PDF' in name: ifPDF = True
                    TextList.append(name)
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if ifPWDR: TextList.insert(1,'All PWDR')
            if ifIMG: TextList.insert(1,'All IMG')
            if ifHKLF: TextList.insert(1,'All HKLF')
            if ifPDF: TextList.insert(1,'All PDF')                
            dlg = wx.MultiChoiceDialog(self, 'Which data to delete?', 'Delete data', TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: DelList.append(TextList[i])
                    if 'All Data' in DelList:
                        DelList = [item for item in TextList if item[:3] != 'All']
                    elif 'All PWDR' in DelList:
                        DelList = [item for item in TextList if item[:4] == 'PWDR']
                    elif 'All IMG' in DelList:
                        DelList = [item for item in TextList if item[:3] == 'IMG']
                    elif 'All HKLF' in DelList:
                        DelList = [item for item in TextList if item[:4] == 'HKLF']
                    elif 'All PDF' in DelList:
                        DelList = [item for item in TextList if item[:3] == 'PDF']
                    item, cookie = self.PatternTree.GetFirstChild(self.root)
                    while item:
                        if self.PatternTree.GetItemText(item) in DelList: DelItemList.append(item)
                        item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
                    for item in DelItemList:
                        self.PatternTree.Delete(item)
                    self.PickId = 0
                    wx.CallAfter(G2plt.PlotPatterns,self,True)                        #so plot gets updated
            finally:
                dlg.Destroy()

    def OnFileOpen(self, event, filename=None):
        '''Reads in a GSAS-II .gpx project file in response to the
        File/Open Project menu button
        '''
        result = wx.ID_OK
        Id = 0
        if self.PatternTree.GetChildrenCount(self.root,False):
            if self.dataFrame:
                self.dataFrame.Clear() 
            dlg = wx.MessageDialog(
                self,
                'Do you want to overwrite the current project? '
                'Any unsaved changes will be lost. Press OK to continue.',
                'Overwrite?',  wx.OK | wx.CANCEL)
            try:
                result = dlg.ShowModal()
                if result == wx.ID_OK:
                    self.PatternTree.DeleteChildren(self.root)
                    self.GSASprojectfile = ''
                    if self.HKL: self.HKL = []
                    if self.G2plotNB.plotList:
                        self.G2plotNB.clear()
            finally:
                dlg.Destroy()
        if result != wx.ID_OK: return

        if not filename:
            if self.dataDisplay: self.dataDisplay.Destroy()
            dlg = wx.FileDialog(self, 'Choose GSAS-II project file', '.', '', 
                'GSAS-II project file (*.gpx)|*.gpx',wx.OPEN|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() != wx.ID_OK: return
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.dirname = dlg.GetDirectory()
            finally:
                dlg.Destroy()
        else:
            self.GSASprojectfile = os.path.splitext(filename)[0]+'.gpx'
            self.dirname = os.path.split(filename)[0]

        G2IO.ProjFileOpen(self)
        self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
        self.PatternTree.Expand(self.root)
        self.HKL = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item and not Id:
            name = self.PatternTree.GetItemText(item)
            if name[:4] in ['PWDR','HKLF','IMG ','PDF ']:
                Id = item
            elif name == 'Controls':
                data = self.PatternTree.GetItemPyData(item)
                if data:
                    for item in self.Refine: item.Enable(True)
                    for item in self.SeqRefine: item.Enable(True)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if Id:
            self.PatternTree.SelectItem(Id)
        self.CheckNotebook()
        os.chdir(self.dirname)           # to get Mac/Linux to change directory!

    def OnFileClose(self, event):
        '''Clears the data tree in response to the
        File/New Project menu button. User is given option to save
        the project.
        '''
        if self.dataFrame:
            self.dataFrame.Clear()
            self.dataFrame.SetLabel('GSAS-II data display') 
        dlg = wx.MessageDialog(self, 'Save current project?', ' ', wx.YES | wx.NO | wx.CANCEL)
        try:
            result = dlg.ShowModal()
            if result == wx.ID_OK:
                self.OnFileSaveMenu(event)
            if result != wx.ID_CANCEL:
                self.GSASprojectfile = ''
                self.PatternTree.SetItemText(self.root,'Loaded Data: ')
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
        finally:
            dlg.Destroy()

    def OnFileSave(self, event):
        '''Save the current project in response to the
        File/Save Project menu button
        '''
        
        if self.GSASprojectfile: 
            self.PatternTree.SetItemText(self.root,'Loaded Data: '+self.GSASprojectfile)
            G2IO.ProjFileSave(self)
        else:
            self.OnFileSaveas(event)

    def OnFileSaveas(self, event):
        '''Save the current project in response to the
        File/Save as menu button
        '''
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file name', '.', '', 
            'GSAS-II project file (*.gpx)|*.gpx',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.PatternTree.SetItemText(self.root,'Saving project as'+self.GSASprojectfile)
                self.SetTitle("GSAS-II data tree: "+os.path.split(self.GSASprojectfile)[1])
                G2IO.ProjFileSave(self)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
        finally:
            dlg.Destroy()

    def ExitMain(self, event):
        '''Called if the main window is closed'''
        if self.undofile:
            os.remove(self.undofile)
        sys.exit()
        
    def OnFileExit(self, event):
        '''Called in response to the File/Quit menu button'''
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def OnExportPatterns(self,event):
        names = ['All']
        exports = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            name = self.PatternTree.GetItemText(item)
            if 'PWDR' in name:
                names.append(name)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if names:
            dlg = wx.MultiChoiceDialog(self,'Select','Powder patterns to export',names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                if sel[0] == 0:
                    exports = names[1:]
                else:
                    for x in sel:
                        exports.append(names[x])
            dlg.Destroy()
        if exports:
            dlg = wx.FileDialog(self, 'Choose output powder file name', '.', '', 
                'GSAS fxye file (*.fxye)|*.fxye|xye file (*.xye)|*.xye',
                wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    powderfile = dlg.GetPath()
                    powderfile = G2IO.FileDlgFixExt(dlg,powderfile)
                    if 'fxye' in powderfile:
                        G2IO.powderFxyeSave(self,exports,powderfile)
                    else:       #just xye
                        G2IO.powderXyeSave(self,exports,powderfile)
            finally:
                dlg.Destroy()
        
    def OnExportPeakList(self,event):
        dlg = wx.FileDialog(self, 'Choose output peak list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item:
                    name = self.PatternTree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.PatternTree.GetFirstChild(item)
                        while item2:
                            name2 = self.PatternTree.GetItemText(item2)
                            if name2 == 'Peak List':
                                peaks = self.PatternTree.GetItemPyData(item2)
                                file.write("%s \n" % (name+' Peak List'))                
                                for peak in peaks:
                                    file.write("%10.5f %12.2f %10.3f %10.3f \n" % \
                                        (peak[0],peak[2],peak[4],peak[6]))
                            item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportHKL(self,event):
        dlg = wx.FileDialog(self, 'Choose output reflection list file name', '.', '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item:
                    name = self.PatternTree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.PatternTree.GetFirstChild(item)
                        while item2:
                            name2 = self.PatternTree.GetItemText(item2)
                            if name2 == 'Reflection Lists':
                                data = self.PatternTree.GetItemPyData(item2)
                                phases = data.keys()
                                for phase in phases:
                                    peaks = data[phase]
                                    file.write("%s %s %s \n" % (name,phase,' Reflection List'))
                                    file.write('%s \n'%(' h  k  l  m  d-space 2-theta wid F**2'))                
                                    for peak in peaks:
                                        FWHM = G2pwd.getgamFW(peak[7],peak[6])/50.      #to get delta-2-theta in deg.
                                        file.write(" %3d %3d %3d %3d %10.5f %10.5f %10.5f %10.3f \n" % \
                                            (int(peak[0]),int(peak[1]),int(peak[2]),int(peak[3]),peak[4],peak[5],FWHM,peak[8]))
                            item2, cookie2 = self.PatternTree.GetNextChild(item, cookie2)                            
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportPDF(self,event):
        #need S(Q) and G(R) to be saved here - probably best from selection?
        names = ['All']
        exports = []
        item, cookie = self.PatternTree.GetFirstChild(self.root)
        while item:
            name = self.PatternTree.GetItemText(item)
            if 'PDF' in name:
                names.append(name)
            item, cookie = self.PatternTree.GetNextChild(self.root, cookie)
        if names:
            dlg = wx.MultiChoiceDialog(self,'Select','PDF patterns to export',names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                if sel[0] == 0:
                    exports = names[1:]
                else:
                    for x in sel:
                        exports.append(names[x])
            dlg.Destroy()
        if exports:
            G2IO.PDFSave(self,exports)
        
    def OnMakePDFs(self,event):
        '''Calculates PDFs
        '''
        tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
        TextList = ['All PWDR']
        PDFlist = []
        Names = []
        if self.PatternTree.GetCount():
            id, cookie = self.PatternTree.GetFirstChild(self.root)
            while id:
                name = self.PatternTree.GetItemText(id)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append(name)
                id, cookie = self.PatternTree.GetNextChild(self.root, cookie)
            if len(TextList) == 1:
                self.ErrorDialog('Nothing to make PDFs for','There must be at least one "PWDR" pattern')
                return
            dlg = wx.MultiChoiceDialog(self,'Make PDF controls','Make PDF controls for:',TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    for i in result: PDFlist.append(TextList[i])
                    if 0 in result:
                        PDFlist = [item for item in TextList if item[:4] == 'PWDR']                        
                    for item in PDFlist:
                        PWDRname = item[4:]
                        Id = self.PatternTree.AppendItem(parent=self.root,text='PDF '+PWDRname)
                        Data = {
                            'Sample':{'Name':item,'Mult':1.0,'Add':0.0},
                            'Sample Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container':{'Name':'','Mult':-1.0,'Add':0.0},
                            'Container Bkg.':{'Name':'','Mult':-1.0,'Add':0.0},'ElList':{},
                            'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':10.0,
                            'DetType':'Image plate','ObliqCoeff':0.2,'Ruland':0.025,'QScaleLim':[0,100],
                            'Lorch':True,}
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='PDF Controls'),Data)
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='I(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='S(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='F(Q)'+PWDRname),[])        
                        self.PatternTree.SetItemPyData(self.PatternTree.AppendItem(Id,text='G(R)'+PWDRname),[])        
                for item in self.ExportPDF: item.Enable(True)
            finally:
                dlg.Destroy()
                
    def GetPWDRdatafromTree(self,PWDRname):
        ''' Returns powder data from GSASII tree

        :param str PWDRname: a powder histogram name as obtained from
          :meth:`GSASIIstruct.GetHistogramNames`

        :returns: PWDRdata = powder data dictionary with
          Powder data arrays, Limits, Instrument Parameters,
          Sample Parameters            
        '''
        PWDRdata = {}
        try:
            PWDRdata.update(self.PatternTree.GetItemPyData(PWDRname)[0])            #wtFactor + ?
        except ValueError:
            PWDRdata['wtFactor'] = 1.0
        PWDRdata['Data'] = self.PatternTree.GetItemPyData(PWDRname)[1]          #powder data arrays
        PWDRdata['Limits'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Limits'))
        PWDRdata['Background'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Background'))
        PWDRdata['Instrument Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Instrument Parameters'))
        PWDRdata['Sample Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Sample Parameters'))
        PWDRdata['Reflection Lists'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,PWDRname,'Reflection Lists'))
        if 'ranId' not in PWDRdata:  # patch, add a random Id
            PWDRdata['ranId'] = ran.randint(0,sys.maxint)
        if 'ranId' not in PWDRdata['Sample Parameters']:  # I hope this becomes obsolete at some point
            PWDRdata['Sample Parameters']['ranId'] = PWDRdata['ranId']
        return PWDRdata

    def GetHKLFdatafromTree(self,HKLFname):
        ''' Returns single crystal data from GSASII tree

        :param str HKLFname: a single crystal histogram name as obtained
          from
          :meth:`GSASIIstruct.GetHistogramNames`

        :returns: HKLFdata = single crystal data list of reflections

        '''
        HKLFdata = {}
        HKLFdata.update(self.PatternTree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        try:
#            HKLFdata.update(self.PatternTree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        except ValueError:
#            HKLFdata['wtFactor'] = 1.0
        HKLFdata['Data'] = self.PatternTree.GetItemPyData(HKLFname)[1]
        HKLFdata['Instrument Parameters'] = self.PatternTree.GetItemPyData(G2gd.GetPatternTreeItemId(self,HKLFname,'Instrument Parameters'))
        return HKLFdata
        
    def GetPhaseData(self):
        '''Returns a dict with defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseData` also exists to
        get same info from GPX file.
        '''
        phaseData = {}
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            print 'no phases found in GetPhaseData'
            sub = None
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phaseName = self.PatternTree.GetItemText(item)
                phaseData[phaseName] =  self.PatternTree.GetItemPyData(item)
                if 'ranId' not in phaseData[phaseName]:
                    phaseData[phaseName]['ranId'] = ran.randint(0,sys.maxint)          
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseData

    def GetPhaseInfofromTree(self):
        '''Get the phase names and their rId values,
        also the histograms used in each phase. 

        :returns: (phaseRIdList, usedHistograms) where

          * 
          ,phaseRIdList,usedHistograms
        '''
        phaseRIdList = []
        usedHistograms = {}
        sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phaseName = self.PatternTree.GetItemText(item)
                ranId = self.PatternTree.GetItemPyData(item).get('ranId')
                if ranId: phaseRIdList.append(ranId)
                data = self.PatternTree.GetItemPyData(item)
                UseList = data['Histograms']
                usedHistograms[phaseName] = UseList.keys()
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseRIdList,usedHistograms

    def GetPhaseNames(self):
        '''Returns a list of defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseNames` also exists to
        get same info from GPX file.
        '''
        phaseNames = []
        if G2gd.GetPatternTreeItemId(self,self.root,'Phases'):
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Phases')
        else:
            print 'no phases found in GetPhaseNames'
            sub = None
        if sub:
            item, cookie = self.PatternTree.GetFirstChild(sub)
            while item:
                phase = self.PatternTree.GetItemText(item)
                phaseNames.append(phase)
                item, cookie = self.PatternTree.GetNextChild(sub, cookie)
        return phaseNames
    
    def GetHistogramNames(self,hType):
        """ Returns a list of histogram names found in the GSASII data tree
        Note routine :func:`GSASIIstrIO.GetHistogramNames` also exists to
        get same info from GPX file.
        
        :param str hType: list of histogram types
        :return: list of histogram names
        
        """
        HistogramNames = []
        if self.PatternTree.GetCount():
            item, cookie = self.PatternTree.GetFirstChild(self.root)
            while item:
                name = self.PatternTree.GetItemText(item)
                if name[:4] in hType:
                    HistogramNames.append(name)        
                item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                

        return HistogramNames

                    
    def GetUsedHistogramsAndPhasesfromTree(self):
        ''' Returns all histograms that are found in any phase
        and any phase that uses a histogram.
        This also assigns numbers to used phases and histograms by the
        order they appear in the file. 
        Note routine :func:`GSASIIstrIO.GetUsedHistogramsAndPhasesfromTree` also exists to
        get same info from GPX file.

        :returns: (Histograms,Phases)

            * Histograms = dictionary of histograms as {name:data,...}
            * Phases = dictionary of phases that use histograms
        '''
        Histograms = {}
        Phases = {}
        phaseNames = self.GetPhaseNames()
        phaseData = self.GetPhaseData()
        histoList = self.GetHistogramNames(['PWDR','HKLF'])

        for phase in phaseData:
            Phase = phaseData[phase]
            pId = phaseNames.index(phase)
            Phase['pId'] = pId
            if Phase['Histograms']:
                if phase not in Phases:
                    Phases[phase] = Phase
                for hist in Phase['Histograms']:
                    if 'Use' not in Phase['Histograms'][hist]:      #patch: add Use flag as True
                        Phase['Histograms'][hist]['Use'] = True         
                    if hist not in Histograms and Phase['Histograms'][hist]['Use']:
                        item = G2gd.GetPatternTreeItemId(self,self.root,hist)
                        if item:
                            if 'PWDR' in hist[:4]: 
                                Histograms[hist] = self.GetPWDRdatafromTree(item)
                            elif 'HKLF' in hist[:4]:
                                Histograms[hist] = self.GetHKLFdatafromTree(item)
                            hId = histoList.index(hist)
                            Histograms[hist]['hId'] = hId
                        else: # would happen if a referenced histogram were renamed or deleted
                            print('For phase "'+str(phase)+
                                  '" unresolved reference to histogram "'+str(hist)+'"')
        #G2obj.IndexAllIds(Histograms=Histograms,Phases=Phases)
        G2obj.IndexAllIds(Histograms=Histograms,Phases=phaseData)
        return Histograms,Phases
        
    def MakeLSParmDict(self):
        '''Load all parameters used for computation from the tree into a
        dict of paired values [value, refine flag]. Note that this is
        different than the parmDict used in the refinement, which only has
        values.

        :returns: (parmDict,varyList) where:

         * parmDict is a dict with values and refinement flags
           for each parameter and
         * varyList is a list of variables (refined parameters).
        '''
        parmDict = {}
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        for phase in Phases:
            if 'pId' not in Phases[phase]:
                self.ErrorDialog('View parameter error','You must run least squares at least once')
                raise Exception,'No pId for phase '+str(phase)
        rigidbodyDict = self.PatternTree.GetItemPyData(   
            G2gd.GetPatternTreeItemId(self,self.root,'Rigid bodies'))
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable = G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)        
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,Print=False)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,Print=False)
        varyList = rbVary+phaseVary+hapVary+histVary
        parmDict.update(rbDict)
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        for parm in parmDict:
            if parm.split(':')[-1] in ['Azimuth','Gonio. radius','Lam1','Lam2',
                'Omega','Chi','Phi','nDebye','nPeaks']:
                parmDict[parm] = [parmDict[parm],'-']
            elif parm.split(':')[-2] in ['Ax','Ay','Az','SHmodel','SHord']:
                parmDict[parm] = [parmDict[parm],'-']
            elif parm in varyList:
                parmDict[parm] = [parmDict[parm],'T']
            else:
                parmDict[parm] = [parmDict[parm],'F']
        return parmDict,varyList

    def ShowLSParms(self,event):
        '''Displays a window showing all parameters in the refinement.
        Called from the Calculate/View LS Parms menu.
        '''
        parmDict,varyList = self.MakeLSParmDict()
        parmValDict = {}
        for i in parmDict:
            parmValDict[i] = parmDict[i][0]
            
        reqVaryList = tuple(varyList) # save requested variables
        try:
            # process constraints
            sub = G2gd.GetPatternTreeItemId(self,self.root,'Constraints') 
            Constraints = self.PatternTree.GetItemPyData(sub)
            constList = []
            for item in Constraints:
                if item.startswith('_'): continue
                constList += Constraints[item]
            G2mv.InitVars()
            constrDict,fixedList,ignored = G2stIO.ProcessConstraints(constList)
            groups,parmlist = G2mv.GroupConstraints(constrDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList)
            G2mv.Map2Dict(parmValDict,varyList)
        except:
            pass
        dlg = G2gd.ShowLSParms(self,'Least Squares Parameters',parmValDict,varyList,reqVaryList)
        dlg.ShowModal()
        dlg.Destroy()
        
    def OnRefine(self,event):
        '''Perform a refinement.
        Called from the Calculate/Refine menu.
        '''        
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.CheckConstraints(self.GSASprojectfile)
        if errmsg:
            print('Error in constraints:\n'+errmsg+
                  '\nRefinement not possible')
            self.ErrorDialog('Constraint Error',
                             'Error in constraints:\n'+errmsg+
                             '\nRefinement not possible')
            return
        if warnmsg:
            print('Conflict between refinment flag settings and constraints:\n'+
                  warnmsg+'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                             'Conflict between refinment flag settings and constraints:\n'+
                             warnmsg+
                             '\nRefinement not possible')
            return
        #works - but it'd be better if it could restore plots
        dlg = wx.ProgressDialog('Residual','All data Rw =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        dlg.SetSize(Size)
        Rw = 100.00
        try:
            Rw = G2stMn.Refine(self.GSASprojectfile,dlg)
        finally:
            dlg.Destroy()
        oldId =  self.PatternTree.GetSelection()
        oldName = self.PatternTree.GetItemText(oldId)
        parentId = self.PatternTree.GetItemParent(oldId)
        parentName = ''
        if parentId:
            parentName = self.PatternTree.GetItemText(parentId)
        dlg = wx.MessageDialog(self,'Load new result?','Refinement results, Rw =%.3f'%(Rw),wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id = 0
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
                G2IO.ProjFileOpen(self)
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item and not Id:
                    name = self.PatternTree.GetItemText(item)
                    if name[:4] in ['PWDR','HKLF']:
                        Id = item
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                if parentName:
                    parentId = G2gd.GetPatternTreeItemId(self, self.root, parentName)
                    if parentId:
                        itemId = G2gd.GetPatternTreeItemId(self, parentId, oldName)
                    else:
                        itemId = G2gd.GetPatternTreeItemId(self, self.root, oldName)
                    self.PatternTree.SelectItem(itemId)
                elif Id:
                    self.PatternTree.SelectItem(Id)
        finally:
            dlg.Destroy()

    def OnSeqRefine(self,event):
        '''Perform a sequential refinement.
        Called from the Calculate/Sequential refine menu.
        '''        
        Id = G2gd.GetPatternTreeItemId(self,self.root,'Sequential results')
        if not Id:
            Id = self.PatternTree.AppendItem(self.root,text='Sequential results')
            self.PatternTree.SetItemPyData(Id,{})            
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.CheckConstraints(self.GSASprojectfile)
        if errmsg:
            print('Error in constraints:\n'+errmsg+
                  '\nRefinement not possible')
            self.ErrorDialog('Constraint Error',
                             'Error in constraints:\n'+errmsg+
                             '\nRefinement not possible')
            return
        if warnmsg:
            print('Conflict between refinment flag settings and constraints:\n'+
                  warnmsg+'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                             'Conflict between refinment flag settings and constraints:\n'+
                             warnmsg+'\nRefinement not possible')
            return
        dlg = wx.ProgressDialog('Residual for histogram 0','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        screenSize = wx.ClientDisplayRect()
        Size = dlg.GetSize()
        Size = (int(Size[0]*1.2),Size[1]) # increase size a bit along x
        dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
        dlg.SetSize(Size)
        try:
            G2stMn.SeqRefine(self.GSASprojectfile,dlg)
        finally:
            dlg.Destroy()        
        dlg = wx.MessageDialog(self,'Load new result?','Refinement results',wx.OK|wx.CANCEL)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Id = 0
                self.PatternTree.DeleteChildren(self.root)
                if self.HKL: self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
                G2IO.ProjFileOpen(self)
                item, cookie = self.PatternTree.GetFirstChild(self.root)
                while item and not Id:
                    name = self.PatternTree.GetItemText(item)
                    if name[:4] in ['PWDR','HKLF']:
                        Id = item
                    item, cookie = self.PatternTree.GetNextChild(self.root, cookie)                
                if Id:
                    self.PatternTree.SelectItem(Id)
        finally:
            dlg.Destroy()
        
    def ErrorDialog(self,title,message,parent=None, wtype=wx.OK):
        'Display an error message'
        result = None
        if parent is None:
            dlg = wx.MessageDialog(self, message, title,  wtype)
        else:
            dlg = wx.MessageDialog(parent, message, title,  wtype)
            dlg.CenterOnParent() # not working on Mac
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        return result

class GSASIImain(wx.App):
    '''Defines a wxApp for GSAS-II

    Creates a wx frame (self.main) which contains the display of the
    data tree.
    '''
    def OnInit(self):
        '''Called automatically when the app is created.'''
        self.main = GSASII(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True
    def MacOpenFile(self, filename):
        '''Called on Mac every time a file is dropped on the app when it is running,
        treat this like a File/Open project menu action.
        Should be ignored on other platforms
        '''
        self.main.OnFileOpen(None,filename)

def main():
    '''Start up the GSAS-II application'''
    #application = GSASIImain() # don't redirect output, someday we
    # may want to do this if we can 
    application = GSASIImain(0)
    if wxInspector:
        import wx.lib.inspection as wxeye
        wxeye.InspectionTool().Show()

    #application.main.OnRefine(None)
    application.MainLoop()
    
if __name__ == '__main__':
    main() # start the GUI
