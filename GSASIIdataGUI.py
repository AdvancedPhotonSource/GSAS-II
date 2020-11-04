# -*- coding: utf-8 -*-
#GSASIIdataGUI - Main GUI routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIdataGUI: Main GSAS-II GUI*
------------------------------------

Module that defines GUI routines and classes for the main GUI Frame (window)
and the main routines that define the GSAS-II tree panel and much of the
data editing panel. 
'''
from __future__ import division, print_function
import platform
import time
import math
import random as ran
import copy
import sys
import os
import inspect
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    try:
        import _pickle as cPickle
    except:
        print('Warning: failed to import the optimized Py3 pickle (_pickle)')
        import pickle as cPickle
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
try:
    import OpenGL as ogl
except ImportError:
    pass
import scipy as sp
import scipy.optimize as so
try:
    import wx
    import wx.grid as wg
    #import wx.wizard as wz
    #import wx.aui
    import wx.lib.scrolledpanel as wxscroll
except ImportError:
    pass
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIImath as G2mth
import GSASIIIO as G2IO
import GSASIIfiles as G2fil
import GSASIIstrIO as G2stIO
import GSASIIlattice as G2lat
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIIconstrGUI as G2cnstG
import GSASIIrestrGUI as G2restG
import GSASIIobj as G2obj
import GSASIIexprGUI as G2exG
import GSASIIlog as log
import GSASIIctrlGUI as G2G
import GSASIIElem as G2elem
import GSASIIpwd as G2pwd
import GSASIIstrMain as G2stMn
import defaultIparms as dI
import GSASIIfpaGUI as G2fpa

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)

# Define short names for convenience
WACV = wx.ALIGN_CENTER_VERTICAL
VERY_LIGHT_GREY = wx.Colour(240,240,240)
DULL_YELLOW = (230,230,190)

# define Ids for wx menu items
commonTrans = {'abc':np.eye(3),'a-cb':np.array([[1.,0.,0.],[0.,0.,-1.],[0.,1.,0.]]),
    'ba-c':np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,-1.]]),'-cba':np.array([[0.,0.,-1.],[0.,1.,0.],[1.,0.,0.]]),
    'bca':np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]]),'cab':np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]]),
    'R->H':np.array([[1.,-1.,0.],[0.,1.,-1.],[1.,1.,1.]]),'H->R':np.array([[2./3,1./3,1./3],[-1./3,1./3,1./3],[-1./3,-2./3,1./3]]),
    'P->A':np.array([[-1.,0.,0.],[0.,-1.,1.],[0.,1.,1.]]),'R->O':np.array([[-1.,0.,0.],[0.,-1.,0.],[0.,0.,1.]]),
    'P->B':np.array([[-1.,0.,1.],[0.,-1.,0.],[1.,0.,1.]]),'B->P':np.array([[-.5,0.,.5],[0.,-1.,0.],[.5,0.,.5]]),
    'P->C':np.array([[1.,1.,0.],[1.,-1.,0.],[0.,0.,-1.]]),'C->P':np.array([[.5,.5,0.],[.5,-.5,0.],[0.,0.,-1.]]),
    'P->F':np.array([[-1.,1.,1.],[1.,-1.,1.],[1.,1.,-1.]]),'F->P':np.array([[0.,.5,.5],[.5,0.,.5],[.5,.5,0.]]),   
    'P->I':np.array([[0.,1.,1.],[1.,0.,1.],[1.,1.,0.]]),'I->P':np.array([[-.5,.5,.5],[.5,-.5,.5],[.5,.5,-.5]]),    
    'A->P':np.array([[-1.,0.,0.],[0.,-.5,.5],[0.,.5,.5]]),'O->R':np.array([[-1.,0.,0.],[0.,-1.,0.],[0.,0.,1.]]), 
    'abc*':np.eye(3), }
commonNames = ['abc','bca','cab','a-cb','ba-c','-cba','P->A','A->P','P->B','B->P','P->C','C->P',
    'P->I','I->P','P->F','F->P','H->R','R->H','R->O','O->R','abc*','setting 1->2']          #don't put any new ones after the setting one!

def SetDefaultDData(dType,histoName,NShkl=0,NDij=0):
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
        
def GetDisplay(pos):
    '''Gets display number (0=main display) for window position (pos). If pos outside all displays
    returns None
    '''
    displays = np.array([list(wx.Display(i).GetGeometry()) for i in range(wx.Display.GetCount())])
    for ip,display in enumerate(displays):
        display[2:3] += display[0:1]
        if (display[0] < pos[0] < display[2]) and (display[1] < pos[1] < display[3]):
            return ip
    return None

################################################################################
#### class definitions used for main GUI
################################################################################
               
class MergeDialog(wx.Dialog):
    ''' HKL transformation & merge dialog
    
    :param wx.Frame parent: reference to parent frame (or None)
    :param data: HKLF data
    
    '''        
               
    def __init__(self,parent,data):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,'Setup HKLF merge', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.data = data
        self.Super = data[1]['Super']
        if self.Super:
            self.Trans = np.eye(4)
        else:
            self.Trans = np.eye(3)
        self.Cent = 'noncentrosymmetric'
        self.Laue = '1'
        self.Class = 'triclinic'
        self.Common = 'abc'
        self.Draw()
        
    def Draw(self):
                
        def OnCent(event):
            Obj = event.GetEventObject()
            self.Cent = Obj.GetValue()
            self.Laue = ''
            wx.CallAfter(self.Draw)
            
        def OnLaue(event):
            Obj = event.GetEventObject()
            self.Laue = Obj.GetValue()
            wx.CallAfter(self.Draw)
            
        def OnClass(event):
            Obj = event.GetEventObject()
            self.Class = Obj.GetValue()
            self.Laue = ''
            wx.CallAfter(self.Draw)
            
        def OnCommon(event):
            Obj = event.GetEventObject()
            self.Common = Obj.GetValue()
            self.Trans = commonTrans[self.Common]
            wx.CallAfter(self.Draw)
       
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        MatSizer = wx.BoxSizer(wx.HORIZONTAL)
        transSizer = wx.BoxSizer(wx.VERTICAL)
        transSizer.Add(wx.StaticText(self.panel,label=" HKL Transformation matrix: M*H = H'"))
        if self.Super:
            Trmat = wx.FlexGridSizer(4,4,0,0)
        else:
            commonSizer = wx.BoxSizer(wx.HORIZONTAL)
            commonSizer.Add(wx.StaticText(self.panel,label=' Common transformations: '),0,WACV)
            common = wx.ComboBox(self.panel,value=self.Common,choices=commonNames[:-2], #not the last two!
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            common.Bind(wx.EVT_COMBOBOX,OnCommon)
            commonSizer.Add(common,0,WACV)
            transSizer.Add(commonSizer)
            Trmat = wx.FlexGridSizer(3,3,0,0)
        for iy,line in enumerate(self.Trans):
            for ix,val in enumerate(line):
                item = G2G.ValidatedTxtCtrl(self.panel,self.Trans[iy],ix,nDig=(10,3),size=(65,25))
                Trmat.Add(item)
        transSizer.Add(Trmat)
        MatSizer.Add((10,0),0)
        MatSizer.Add(transSizer)
        mainSizer.Add(MatSizer)
        laueClass = ['triclinic','monoclinic','orthorhombic','trigonal(H)','tetragonal','hexagonal','cubic']
        centroLaue = {'triclinic':['-1',],'monoclinic':['2/m','1 1 2/m','2/m 1 1',],
            'orthorhombic':['m m m',],'trigonal(H)':['-3','-3 m 1','-3 1 m',],    \
            'tetragonal':['4/m','4/m m m',],'hexagonal':['6/m','6/m m m',],'cubic':['m 3','m 3 m']}
        noncentroLaue = {'triclinic':['1',],'monoclinic':['2','2 1 1','1 1 2','m','m 1 1','1 1 m',],
            'orthorhombic':['2 2 2','m m 2','m 2 m','2 m m',],
            'trigonal(H)':['3','3 1 2','3 2 1','3 m 1','3 1 m',],
            'tetragonal':['4','-4','4 2 2','4 m m','-4 2 m','-4 m 2',], \
            'hexagonal':['6','-6','6 2 2','6 m m','-6 m 2','-6 2 m',],'cubic':['2 3','4 3 2','-4 3 m']}
        centChoice = ['noncentrosymmetric','centrosymmetric']
        mainSizer.Add(wx.StaticText(self.panel,label=' Select Laue class for new lattice:'),0)
        Class = wx.ComboBox(self.panel,value=self.Class,choices=laueClass,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Class.Bind(wx.EVT_COMBOBOX,OnClass)
        mainSizer.Add(Class,0)
        mainSizer.Add(wx.StaticText(self.panel,label=' Target Laue symmetry:'),0)
        Cent = wx.ComboBox(self.panel,value=self.Cent,choices=centChoice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Cent.Bind(wx.EVT_COMBOBOX,OnCent)
        mergeSizer = wx.BoxSizer(wx.HORIZONTAL)
        mergeSizer.Add(Cent,0,WACV)
        mergeSizer.Add((10,0),0)
        Choice = centroLaue[self.Class]
        if 'non' in self.Cent:
            Choice = noncentroLaue[self.Class]
        Laue = wx.ComboBox(self.panel,value=self.Laue,choices=Choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        Laue.Bind(wx.EVT_COMBOBOX,OnLaue)
        mergeSizer.Add(Laue,0,WACV)
        mainSizer.Add(mergeSizer)

        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(self.panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        if self.Laue:
            btnSizer.Add(OkBtn)
            btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)
        
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        
    def GetSelection(self):
        return self.Trans,self.Cent,self.Laue

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)
        
def GUIpatches():
    'Misc fixes that only needs to be done when running a GUI'
    try:  # patch for LANG environment var problem on occasional OSX machines
        import locale
        locale.getdefaultlocale()
    except ValueError:
        print('Fixing location (see https://github.com/matplotlib/matplotlib/issues/5420.)')
        os.environ['LC_ALL'] = 'en_US.UTF-8'
        locale.getdefaultlocale()
    try:
        import OpenGL
        OpenGL  # avoids unused package error
    except ImportError:
        print('*******************************************************')
        print('PyOpenGL is missing from your python installation')
        print('     - we will try to install it')
        print('*******************************************************')
        def install_with_easyinstall(package):
            try: 
                print ("trying a system-wide PyOpenGl install")
                easy_install.main(['-f',os.path.split(__file__)[0],package])
                return
            except:
                pass
            try: 
                print ("trying a user level PyOpenGl install")
                easy_install.main(['-f',os.path.split(__file__)[0],'--user',package])
                return
            except:
                print (u"Install of '+package+u' failed. Please report this information:")
                import traceback
                print (traceback.format_exc())
                sys.exit()
        from setuptools.command import easy_install
        install_with_easyinstall('PyOpenGl')
        print('*******************************************************')         
        print('OpenGL has been installed. Restarting GSAS-II')
        print('*******************************************************')         
        loc = os.path.dirname(__file__)
        import subprocess
        subprocess.Popen([sys.executable,os.path.join(loc,'GSASII.py')])
        sys.exit()
    # PATCH: for Mavericks (OS X 10.9.x), wx produces an annoying warning about LucidaGrandeUI.
    # In case stderr has been suppressed there, redirect python error output to stdout. Nobody
    # else should care much about this. 
    sys.stderr = sys.stdout

def ShowVersions():
    '''Show the versions all of required Python packages, etc.
    '''
    import numpy as np
    import scipy as sp
    import wx
    import matplotlib as mpl
    import OpenGL as ogl
    import GSASIIpath
    # print (versions)
    print ("Python module versions loaded:")
    print ("  Python:     %s from %s"%(sys.version.split()[0],sys.executable))
    print ("  wx:         %s"%wx.__version__)
    print ("  matplotlib: %s"%mpl.__version__)
    print ("  numpy:      %s"%np.__version__)
    print ("  scipy:      %s"%sp.__version__)
    print ("  OpenGL:     %s"%ogl.__version__)
    Image = None
    version = '?'
    try:
        from PIL import Image
    except ImportError:
        try:
            import Image
        except ImportError:
            pass
    if Image is None:
        print ("Image module not present; Note that PIL (Python Imaging Library) or pillow is needed for some image operations")
    else:
        # version # can be in various places, try standard dunderscore first
        for ver in '__version__','VERSION','PILLOW_VERSION':
            if hasattr(Image,ver):
                try:
                    version = eval('Image.'+ver)
                    break
                except:
                    pass
        print ("  Image:      %s (PIL or Pillow)"%version)
    print ("  Platform:   %s %s %s"%(sys.platform,platform.architecture()[0],platform.machine()))
    try:
        import mkl
        print ("  Max threads:%s"%mkl.get_max_threads())
    except:
        pass
    rev = GSASIIpath.svnGetRev()
    if rev is None: 
        "no SVN"
    else:
        rev = "SVN version {}".format(rev)
    print ("Latest GSAS-II revision (from .py files): {} ({})\n".format(
        GSASIIpath.GetVersionNumber(),rev))
    # patch 11/2020: warn if GSASII path has not been updated past v4576.
    # For unknown reasons on Mac with gsas2full, there have been checksum
    # errors in the .so files that prevented svn from completing updates.
    # If GSASIIpath.svnChecksumPatch is not present, then the fix for that
    # has not been retrieved, so warn. Keep for a year or so. 
    try:
        GSASIIpath.svnChecksumPatch
    except:
        print('Warning GSAS-II incompletely updated. Please contact toby@anl.gov')
    # end patch

def warnNumpyVersion(application):
    dlg = wx.MessageDialog(application.main,
                'version '+ np.__version__ +
                ' of numpy produces .gpx files that are not compatible with older versions: please upgrade or downgrade numpy to avoid this version',
                'numpy Warning',wx.OK)
    try:
        dlg.CenterOnParent()
        dlg.ShowModal()
    finally:
        dlg.Destroy()
        
###############################################################################
#### GUI creation
###############################################################################
def GSASIImain(application):
    '''Start up the GSAS-II GUI'''                        
    ShowVersions()
    GUIpatches()
    knownVersions = ['2.7','3.6','3.7','3.8']
    if platform.python_version()[:3] not in knownVersions: 
        dlg = wx.MessageDialog(None, 
                'GSAS-II requires Python 2.7.x or 3.6+\n Yours is '+sys.version.split()[0],
                'Python version error',  wx.OK)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()
        sys.exit()

    if platform.python_version()[:3] == '2.7':
        msg = '''The end-of-life for python 2.7 was January 1, 2020. 
We strongly recommend reinstalling GSAS-II from a new installation kit as we may not be able to offer support for operation of GSAS-II in python 2.7. See instructions for details.
'''
        download = ''
        cmds = []
        instructions = 'https://subversion.xray.aps.anl.gov/trac/pyGSAS'
        if sys.platform == "win32":
            download = 'https://subversion.xray.aps.anl.gov/admin_pyGSAS/downloads/gsas2full-Latest-Windows-x86_64.exe'
            instructions = 'https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/SingleStepWindowsIllustrated'
        elif sys.platform == "darwin":
            cmds = ['echo starting download, please wait...',
                    '''echo 'curl "https://subversion.xray.aps.anl.gov/admin_pyGSAS/downloads/gsas2full-Latest-MacOSX-x86_64.sh" > /tmp/g2.sh; bash /tmp/g2.sh' ''',
                    'curl "https://subversion.xray.aps.anl.gov/admin_pyGSAS/downloads/gsas2full-Latest-MacOSX-x86_64.sh" > /tmp/g2.sh; bash /tmp/g2.sh'
                    ]
            instructions = 'https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/MacSingleStepInstallerFigs'
        elif sys.platform.startswith("linux"):
            download = 'https://subversion.xray.aps.anl.gov/admin_pyGSAS/downloads/gsas2full-Latest-Linux-x86_64.sh'
            instructions = 'https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/LinuxSingleStepInstaller'
        else:
            print(u'Unknown platform: '+sys.platform)
        if platform.architecture()[0] != '64bit' and sys.platform == "win32":
            msg += '''\nYou are currently using 32-bit Python. Please check if you are running 32-bit windows or 64-bit windows (use Start/Settings/System/About & look for "System type".
            We recommend using the 64-bit installer if you have 64-bit windows.'''
            download = ''
        elif platform.architecture()[0] != '64bit' and sys.platform.startswith("linux"):
            msg += '''\nYou are using 32-bit Python. We now only package for 64-bit linux.
            If you are running 32-bit linux you will need to install Python yourself.
            See instructions at https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallLinux'''
            instructions = 'https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/InstallLinux'
        dlg = wx.Dialog(None,wx.ID_ANY,'End-Of-Life warning for Python 2.7',
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        txt = wx.StaticText(dlg,wx.ID_ANY,G2G.StripIndents(msg))
        mainSizer.Add(txt)
        txt.Wrap(600)
        dlg.SetSizer(mainSizer)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btnsizer.Add((1,1),1,wx.EXPAND,1)
        OKbtn = wx.Button(dlg, wx.ID_OK,'Continue')
        OKbtn.SetDefault()
        OKbtn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
        btnsizer.Add(OKbtn)

        btn = wx.Button(dlg, wx.ID_ANY,'Show Instructions')
        def openInstructions(event):
            G2G.ShowWebPage(instructions,None)
        btn.Bind(wx.EVT_BUTTON, openInstructions)
        btnsizer.Add(btn)
        if download:
            btn = wx.Button(dlg, wx.ID_ANY,'Start Download')
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_YES))
            btnsizer.Add(btn)
        elif cmds:
            btn = wx.Button(dlg, wx.ID_ANY,'Start Install')
            btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
            btnsizer.Add(btn)

        #btn = wx.Button(dlg, wx.ID_CANCEL)
        #btnsizer.AddButton(btn)
        btnsizer.Add((1,1),1,wx.EXPAND,1)
        #btnsizer.Realize()
        mainSizer.Add((-1,5),1,wx.EXPAND,1)
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER,0)
        mainSizer.Add((-1,10))
        res = 0
        try:
            res = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if res == wx.ID_YES:
            G2G.ShowWebPage(download,None)
            G2G.ShowWebPage(instructions,None)
            wx.Sleep(1)
            dlg = wx.MessageDialog(None,G2G.StripIndents(
                    '''Download has been started in your browser; installation instructions will also be shown in a web page\n\nPress OK to exit GSAS-II, Cancel to continue.'''),
                    'start install',wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:
                sys.exit()
        elif res == wx.ID_CANCEL:
            dlg = wx.MessageDialog(None,G2G.StripIndents(
                    '''Press OK to continue. Instructions will be shown in a web page. 
                    Download and installation will start in the terminal window after you press OK. Respond to questions there.'''),
                    'start install',wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:
                G2G.ShowWebPage(instructions,None)
                GSASIIpath.runScript(cmds, wait=True)    
                sys.exit()
    application.main = GSASII(None)  # application.main is the main wx.Frame (G2frame in most places)
    application.SetTopWindow(application.main)
    # save the current package versions
    application.main.PackageVersions = G2fil.get_python_versions([wx, mpl, np, sp, ogl])
    if np.__version__ == '1.16.0':
        wx.CallLater(100,warnNumpyVersion,application)
    if GSASIIpath.GetConfigValue('wxInspector'):
        import wx.lib.inspection as wxeye
        wxeye.InspectionTool().Show()

    try:
        application.SetAppDisplayName('GSAS-II')
    except:
        pass
    #application.GetTopWindow().SendSizeEvent()
    application.GetTopWindow().Show(True)

################################################################################
#### Create main frame (window) for GUI
################################################################################
class GSASII(wx.Frame):
    '''Define the main GSAS-II frame and its associated menu items.

    :param parent: reference to parent application

    '''                
    def MenuBinding(self,event):
        '''Called when a menu is clicked upon; looks up the binding in table
        '''
        log.InvokeMenuCommand(event.GetId(),self,event)
            
#    def Bind(self,eventtype,handler,*args,**kwargs):
#        '''Override the Bind function so that we can wrap calls that will be logged.
#        
#        N.B. This is a bit kludgy. Menu bindings with an id are wrapped and
#        menu bindings with an object and no id are not. 
#        '''
#        if eventtype == wx.EVT_MENU and 'id' in kwargs:
#            menulabels = log.SaveMenuCommand(kwargs['id'],self,handler)
#            if menulabels:
#                wx.Frame.Bind(self,eventtype,self.MenuBinding,*args,**kwargs)
#                return
#        wx.Frame.Bind(self,eventtype,handler,*args,**kwargs)      
    
    def _Add_FileMenuItems(self, parent):
        '''Add items to File menu
        '''
        item = parent.Append(wx.ID_ANY,'&Open project...\tCtrl+O','Open a GSAS-II project file (*.gpx)')            
        self.Bind(wx.EVT_MENU, self.OnFileOpen, id=item.GetId())
        if sys.platform == "darwin": 
            item = parent.Append(wx.ID_ANY,'&Open in new window...','Open a GSAS-II project file (*.gpx) in a separate process')
            self.Bind(wx.EVT_MENU, self.OnNewGSASII, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Reopen recent...\tCtrl+E','Reopen a previously used GSAS-II project file (*.gpx)')
        self.Bind(wx.EVT_MENU, self.OnFileReopen, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'&Save project\tCtrl+S','Save project under current name')
        self.Bind(wx.EVT_MENU, self.OnFileSave, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Save project as...','Save current project to new file')
        self.Bind(wx.EVT_MENU, self.OnFileSaveas, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'&New project','Create empty new project, saving current is optional')
        self.Bind(wx.EVT_MENU, self.OnFileClose, id=item.GetId())
        item = parent.Append(wx.ID_PREFERENCES,"&Preferences",'')
        self.Bind(wx.EVT_MENU, self.OnPreferences, item)
        if GSASIIpath.whichsvn():
            item = parent.Append(wx.ID_ANY,'Edit proxy...','Edit proxy internet information (used for updates)')
            self.Bind(wx.EVT_MENU, self.EditProxyInfo, id=item.GetId())
        if GSASIIpath.GetConfigValue('debug'):
            def OnIPython(event):
                GSASIIpath.IPyBreak()
            item = parent.Append(wx.ID_ANY,"IPython Console",'')
            self.Bind(wx.EVT_MENU, OnIPython, item)
            def OnwxInspect(event):
                import wx.lib.inspection as wxeye
                wxeye.InspectionTool().Show()
            item = parent.Append(wx.ID_ANY,"wx inspection tool",'')
            self.Bind(wx.EVT_MENU, OnwxInspect, item)
            
        item = parent.Append(wx.ID_EXIT,'Exit\tALT+F4','Exit from GSAS-II')
        self.Bind(wx.EVT_MENU, self.ExitMain, id=item.GetId())
        
    def _Add_DataMenuItems(self,parent):
        '''Add items to Data menu
        '''
        # item = parent.Append(
        #     help='',id=wx.ID_ANY,
        #     kind=wx.ITEM_NORMAL,
        #     text='Read image data...')
        # self.Bind(wx.EVT_MENU, self.OnImageRead, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Read Powder Pattern Peaks...','')
        self.Bind(wx.EVT_MENU, self.OnReadPowderPeaks, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Sum or Average powder data','')
        self.Bind(wx.EVT_MENU, self.OnPwdrSum, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Sum image data','')
        self.Bind(wx.EVT_MENU, self.OnImageSum, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Add new phase','')
        self.Bind(wx.EVT_MENU, self.OnAddPhase, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Delete phase entries','')
        self.Bind(wx.EVT_MENU, self.OnDeletePhase, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Rename data entry',
                'Rename the selected data tree item (PWDR, HKLF or IMG)')
        self.Bind(wx.EVT_MENU, self.OnRenameData, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Delete data entries',
                'Delete selected data items from data tree')
        self.Bind(wx.EVT_MENU, self.OnDataDelete, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'Delete plots','Delete selected plots')
        self.Bind(wx.EVT_MENU, self.OnPlotDelete, id=item.GetId())
        expandmenu = wx.Menu()
        item = parent.AppendSubMenu(expandmenu,'Expand tree items',  
            'Expand items of type in GSAS-II data tree')
        for s in 'all','IMG','PWDR','PDF','HKLF','SASD','REFD':
            if s == 'all':
                help = 'Expand all items in GSAS-II data tree'
            else:
                help = 'Expand '+s+' type items in GSAS-II data tree'
            item = expandmenu.Append(wx.ID_ANY,s,help)
            self.Bind(wx.EVT_MENU,self.ExpandAll,id=item.GetId())
        movemenu = wx.Menu()
        item = parent.AppendSubMenu(movemenu,'Move tree items',  
            'Move items of type items to end of GSAS-II data tree')
        for s in 'IMG','PWDR','PDF','HKLF','SASD','REFD','Phase':
            help = 'Move '+s+' type items to end of GSAS-II data tree'
            item = movemenu.Append(wx.ID_ANY,s,help)
            self.Bind(wx.EVT_MENU,self.MoveTreeItems,id=item.GetId())

    def _Add_CalculateMenuItems(self,parent):
        item = parent.Append(wx.ID_ANY,'Setup PDFs','Create PDF tree entries for selected powder patterns')
        self.MakePDF.append(item)
        self.Bind(wx.EVT_MENU, self.OnMakePDFs, id=item.GetId())
        
        item = parent.Append(wx.ID_ANY,'&View LS parms\tCTRL+L','View least squares parameters')
        self.Bind(wx.EVT_MENU, self.OnShowLSParms, id=item.GetId())
        
        item = parent.Append(wx.ID_ANY,'&Refine\tCTRL+R','Perform a refinement') 
        if len(self.Refine): # extend state for new menus to match main (on mac)
            state = self.Refine[0].IsEnabled()
        else:
            state = False
        item.Enable(state)
        self.Refine.append(item)
        self.Bind(wx.EVT_MENU, self.OnRefine, id=item.GetId())

        item = parent.Append(wx.ID_ANY,'&Run Fprime','X-ray resonant scattering')
        self.Bind(wx.EVT_MENU, self.OnRunFprime, id=item.GetId())
        item = parent.Append(wx.ID_ANY,'&Run Absorb','x-ray absorption')
        self.Bind(wx.EVT_MENU, self.OnRunAbsorb, id=item.GetId())

#        if GSASIIpath.GetConfigValue('debug'): # allow exceptions for debugging
#            item = parent.Append(help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
#                text='tree test')
#            self.Bind(wx.EVT_MENU, self.TreeTest, id=item.GetId())

    def _init_Imports(self):
        '''import all the G2phase*.py & G2sfact*.py & G2pwd*.py files that 
        are found in the path
        '''

        self.ImportPhaseReaderlist = G2fil.LoadImportRoutines('phase','Phase')
        self.ImportSfactReaderlist = G2fil.LoadImportRoutines('sfact','Struct_Factor')
        self.ImportPowderReaderlist = G2fil.LoadImportRoutines('pwd','Powder_Data')
        self.ImportSmallAngleReaderlist = G2fil.LoadImportRoutines('sad','SmallAngle_Data')
        self.ImportReflectometryReaderlist = G2fil.LoadImportRoutines('rfd','Reflectometry_Data')
        self.ImportPDFReaderlist = G2fil.LoadImportRoutines('pdf','PDF_Data')
        self.ImportImageReaderlist = G2fil.LoadImportRoutines('img','Images')
        self.ImportMenuId = {}

    def testSeqRefineMode(self):
        '''Returns the list of histograms included in a sequential refinement or
        an empty list if a standard (non-sequential) refinement.
        Also sets Menu item status depending on mode
        '''
        cId = GetGPXtreeItemId(self,self.root, 'Controls')
        if cId:
            controls = self.GPXtree.GetItemPyData(cId)
            seqSetting = controls.get('Seq Data',[])
        else:
            seqSetting = None
            
        if seqSetting:
            for item in self.Refine:
                item.SetItemLabel('Se&quential refine\tCtrl+R')    #might fail on old wx
            seqMode = True
        else:
            for item in self.Refine:
                item.SetItemLabel('&Refine\tCtrl+R')    #might fail on old wx
            seqMode = False
        for menu,Id in self.ExportSeq:
            menu.Enable(Id,seqMode)
        for menu,Id in self.ExportNonSeq:
            menu.Enable(Id,not seqMode)
        return seqSetting

        
    def PreviewFile(self,filename):
        'utility to confirm we have the right file'
        fp = open(filename,'r')
        rdmsg = u'File '+ filename +u' begins:\n\n'
        try:
            rdmsg += fp.read(80)
            rdmsg += '\n\nDo you want to read this file?'
        except UnicodeDecodeError:
            rdmsg = None
        fp.close()
        if rdmsg is None or not all([ord(c) < 128 and ord(c) != 0 for c in rdmsg]): # show only if ASCII
            rdmsg = u'File '+ filename +u' is a binary file. Do you want to read this file?'
        # it would be better to use something that
        # would resize better, but this will do for now
        dlg = wx.MessageDialog(
            self, rdmsg,
            'Is this the file you want?', 
            wx.YES_NO | wx.ICON_QUESTION,
            )
        dlg.SetSize((700,300)) # does not resize on Mac
        result = wx.ID_NO
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_NO: return True
        return False
    
    def OnImportGeneric(self,reader,readerlist,label,multiple=False,
        usedRanIdList=[],Preview=True,load2Tree=False):
        '''Used for all imports, including Phases, datasets, images...

        Called from :meth:`GSASII.OnImportPhase`, :meth:`GSASII.OnImportImage`,
        :meth:`GSASII.OnImportSfact`, :meth:`GSASII.OnImportPowder`,
        :meth:`GSASII.OnImportSmallAngle` and :meth:'GSASII.OnImportReflectometry`

        Uses reader_objects subclassed from :class:`GSASIIobj.ImportPhase`,
        :class:`GSASIIobj.ImportStructFactor`,
        :class:`GSASIIobj.ImportPowderData`,
        :class:`GSASIIobj.ImportSmallAngleData`
        :class:`GSASIIobj.ImportReflectometryData` or
        :class:`GSASIIobj.ImportImage`.
        If a specific reader is specified, only that method will be called,
        but if no reader is specified, every one that is potentially
        compatible (by file extension) will be tried on the file(s)
        selected in the Open File dialog.

        :param reader_object reader: This will be a reference to
          a particular object to be used to read a file or None,
          if every appropriate reader should be used.

        :param list readerlist: a list of reader objects appropriate for
          the current read attempt. At present, this will be either
          self.ImportPhaseReaderlist, self.ImportSfactReaderlist
          self.ImportPowderReaderlist or self.ImportImageReaderlist
          (defined in _init_Imports from the files found in the path),
          but in theory this list could be tailored.
          Used only when reader is None.

        :param str label: string to place on the open file dialog:
          Open `label` input file

        :param bool multiple: True if multiple files can be selected
          in the file dialog. False is default. At present True is used
          only for reading of powder data.

        :param list usedRanIdList: an optional list of random Ids that 
          have been used and should not be reused

        :param bool Preview: indicates if a preview of the file should
          be shown. Default is True, but set to False for image files
          which are all binary.

        :param bool load2Tree: indicates if the file should be loaded
          into the data tree immediately (used for images only). True
          only when called from :meth:`OnImportImage`; causes return
          value to change to a list of True values rather than
          reader objects.

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
            for extn in sorted(extdict.keys(),key=lambda k: k.lower()):
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
            mode = wx.FD_OPEN|wx.FD_MULTIPLE
        else:
            mode = wx.FD_OPEN
        if len(readerlist) > 1: 
            typ = ' (type to be guessed)'
        else:
            typ = '( type '+readerlist[0].formatName+')'
        filelist = G2G.GetImportFile(self,
                    message="Choose "+label+" input file"+typ,
                    defaultFile="",wildcard=choices,style=mode)
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
        Start = True    #1st time read - clear selections below
        for filename in filelist:
            # is this a zip file?
            if os.path.splitext(filename)[1].lower() == '.zip':
                extractedfile = G2IO.ExtractFileFromZip(filename,parent=self)
                if extractedfile is None: continue # error or Cancel 
                if extractedfile != filename:
                    filename,self.zipfile = extractedfile,filename # now use the file that was created
            # determine which formats are compatible with this file
            primaryReaders = []
            secondaryReaders = []
            for rd in readerlist:
                flag = rd.ExtensionValidator(filename)
                if flag is None: 
                    secondaryReaders.append(rd)
                elif flag:
                    primaryReaders.append(rd)
            if len(secondaryReaders) + len(primaryReaders) == 0 and reader:
                self.ErrorDialog('Not supported','The selected reader cannot read file '+filename)
                return []
            elif len(secondaryReaders) + len(primaryReaders) == 0:
                self.ErrorDialog('No Format','No matching format for file '+filename)
                return []

            fp = None
            msg = ''
            if len(filelist) == 1 and Preview:
                if self.PreviewFile(filename): return []
            self.lastimport = filename # this is probably not what I want to do -- it saves only the
            # last name in a series. See rd.readfilename for a better name.

            # try the file first with Readers that specify the
            # file's extension and later with ones that merely allow it
            errorReport = ''
            for rd in primaryReaders+secondaryReaders:
                if Start:   #clear old bank selections to allow new ones to be selected by user
                    rd.selections = []
                    rd.dnames = []
                rd.ReInitialize() # purge anything from a previous read
                rd.errors = "" # clear out any old errors
                if not rd.ContentsValidator(filename): # rejected on cursory check
                    errorReport += "\n  "+rd.formatName + ' validator error'
                    if rd.errors: 
                        errorReport += ': '+rd.errors
                    continue
                if len(rd.selections)>1 and Start:
                    dlg = G2G.G2MultiChoiceDialog(self,'Dataset Selector','Select data to read from the list below',rd.dnames)
                    if dlg.ShowModal() == wx.ID_OK:
                        rd.selections = dlg.GetSelections()
                    Start = False
                    dlg.Destroy()
                repeat = True
                rdbuffer = {} # create temporary storage for file reader
                block = 0
                while repeat: # loop if the reader asks for another pass on the file
                    block += 1
                    repeat = False
                    rd.objname = os.path.basename(filename)
                    flag = False
                    if GSASIIpath.GetConfigValue('debug'): # allow exceptions for debugging
                        flag = rd.Reader(filename,self,buffer=rdbuffer,blocknum=block,
                            usedRanIdList=usedRanIdList,)
                    else:
                        try:
                            flag = rd.Reader(filename,self,buffer=rdbuffer,
                                blocknum=block,usedRanIdList=usedRanIdList,)
                        except rd.ImportException as detail:
                            rd.errors += "\n  Read exception: "+str(detail)
                        except Exception as detail:
                            import traceback
                            rd.errors += "\n  Unhandled read exception: "+str(detail)
                            rd.errors += "\n  Traceback info:\n"+str(traceback.format_exc())
                    if flag: # this read succeeded
                        if rd.SciPy:        #was default read by scipy; needs 1 time fixes
                            G2IO.EditImageParms(self,rd.Data,rd.Comments,rd.Image,filename)
                            rd.SciPy = False
                        rd.readfilename = filename
                        if load2Tree:   #images only
                            if rd.repeatcount == 1 and not rd.repeat: # skip image number if only one in set
                                rd.Data['ImageTag'] = None
                            else:
                                rd.Data['ImageTag'] = rd.repeatcount
                            rd.Data['formatName'] = rd.formatName
                            if rd.sumfile:
                                rd.readfilename = rd.sumfile
                            # Load generic metadata, as configured
                            G2fil.GetColumnMetadata(rd)
                            G2IO.LoadImage2Tree(rd.readfilename,self,rd.Comments,rd.Data,rd.Npix,rd.Image)
                            rd_list.append(True) # save a stub the result before it is written over
                            del rd.Image
                        else:                                                   
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
                            ' reader reported a warning message:\n\n'+rd.warnings)
                    break # success in reading, try no further
            else:
                if singlereader:
                    msg += '\n'+rd.warnings
                    print(u'The '+ rd.formatName+u' reader was not able to read file '+filename+msg)
                    try:
                        print(u'\n\nError message(s):\n\t'+errorReport)
                    except:
                        pass
                    self.ErrorDialog('Read Error','The '+ rd.formatName+
                        ' reader was not able to read file '+filename+msg)
                else:
                    print('No reader was able to read file '+filename+msg)
                    try:
                        print('\n\nError message(s):\n\t'+errorReport)
                    except:
                        pass
                    self.ErrorDialog('Read Error','No reader was able to read file '+filename+msg)
            if fp: fp.close()
        return rd_list

    def _Add_ImportMenu_Phase(self,parent):
        '''configure the Import Phase menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'Phase','Import phase data')
        for reader in self.ImportPhaseReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPhase, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,'guess format from file','Import phase data, use file to try to determine format')
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
        phaseNameList = list(usedHistograms.keys()) # phase names in use
        usedHKLFhists = [] # used single-crystal histograms
        for p in usedHistograms:
            for h in usedHistograms[p]:
                if h.startswith('HKLF ') and h not in usedHKLFhists:
                    usedHKLFhists.append(h)
                    
                    
        rdlist = self.OnImportGeneric(reqrdr,self.ImportPhaseReaderlist,
            'phase',usedRanIdList=phaseRIdList)
        if len(rdlist) == 0: return
        # for now rdlist is only expected to have one element
        # but below will allow multiple phases to be imported
        # if ever the import routines ever implement multiple phase reads. 
        self.CheckNotebook()
        newPhaseList = []
        for rd in rdlist:
            PhaseName = ''
            dlg = wx.TextEntryDialog(self, 'Enter the name for the new phase',
                'Edit phase name', rd.Phase['General']['Name'],style=wx.OK)
            while PhaseName == '':
                dlg.CenterOnParent()
                if dlg.ShowModal() == wx.ID_OK:
                    PhaseName = dlg.GetValue().strip()
                else:
                    dlg.Destroy()
                    return
            dlg.Destroy()
            # make new phase names unique
            rd.Phase['General']['Name'] = G2obj.MakeUniqueLabel(PhaseName,phaseNameList)
            if rd.Phase['General']['SGData']['SpGrp'] in G2spc.spg2origins:
                wx.MessageDialog(self,' The atom positions for this structure must be in the 2nd setting for the space group '+
                    rd.Phase['General']['SGData']['SpGrp'], 
                    'Is this 2nd setting?',  wx.OK).ShowModal()
                
            PhaseName = rd.Phase['General']['Name'][:]
            newPhaseList.append(PhaseName)
            print(u'Read phase {} from file {}'.format(PhaseName,self.lastimport))
            if not GetGPXtreeItemId(self,self.root,'Phases'):
                sub = self.GPXtree.AppendItem(parent=self.root,text='Phases')
            else:
                sub = GetGPXtreeItemId(self,self.root,'Phases')
            psub = self.GPXtree.AppendItem(parent=sub,text=PhaseName)
            self.GPXtree.SetItemPyData(psub,rd.Phase)
            wx.CallAfter(self.GPXtree.SelectItem,psub) # should call SelectDataTreeItem
            try:
                rd.MPhase['General']['Name'] = G2obj.MakeUniqueLabel(PhaseName+' mag',phaseNameList)
                PhaseName = rd.MPhase['General']['Name'][:]
                newPhaseList.append(PhaseName)
                psub = self.GPXtree.AppendItem(parent=sub,text=PhaseName)
                self.GPXtree.SetItemPyData(psub,rd.MPhase)
                wx.CallAfter(self.GPXtree.SelectItem,psub) # should call SelectDataTreeItem
            except (AttributeError,TypeError):
                pass
            self.GPXtree.Expand(self.root) # make sure phases are seen
            self.GPXtree.Expand(sub) 
            self.GPXtree.Expand(psub)
            self.PickIdText = None
            
            # add constraints imported with phase to tree
            #    at present, constraints are generated only in ISODISTORT_proc in the
            #    CIF import
            if rd.Constraints:
                sub = GetGPXtreeItemId(self,self.root,'Constraints') # was created in CheckNotebook if needed
                Constraints = self.GPXtree.GetItemPyData(sub)                
                for i in rd.Constraints:
                    if type(i) is dict:
                        if '_Explain' not in Constraints: Constraints['_Explain'] = {}
                        Constraints['_Explain'].update(i)
                    else:
                        Constraints['Phase'].append(i)
        if not newPhaseList: return # somehow, no new phases
        # get a list of existing histograms
        PWDRlist = []
        HKLFlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        TextList = PWDRlist + HKLFlist
        if not TextList:
            return          #no histograms
        header = 'Select histogram(s) to add to new phase(s):'
        for phaseName in newPhaseList:
            header += '\n  '+phaseName

        notOK = True
        while notOK:
            result = G2G.ItemSelector(TextList,self,header,header='Add histogram(s)',multiple=True)
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
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.GPXtree.GetFirstChild(sub)
        while item: # loop over (new) phases
            phaseName = self.GPXtree.GetItemText(item)
            data = self.GPXtree.GetItemPyData(item)
            item, cookie = self.GPXtree.GetNextChild(sub, cookie)
            if phaseName not in newPhaseList: continue
            generalData = data['General']
            SGData = generalData['SGData']
            Super = generalData.get('Super',0)
            SuperVec = []
            if Super:
                SuperVec = np.array(generalData['SuperVec'][0])
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for i in result:
                histoName = TextList[i]
                if histoName in HKLFlist:
                    #redo UpdateHKLFdata(histoName) here:
                    Id = GetGPXtreeItemId(self,self.root,histoName)
                    refDict,reflData = self.GPXtree.GetItemPyData(Id)
                    G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                    Super = reflData.get('Super',0)
                    for iref,ref in enumerate(reflData['RefList']):
                        hkl = ref[:3]
                        if Super:
                            H = list(hkl+SuperVec*ref[3])
                        else:
                            H = hkl
                        ref[4+Super] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                        iabsnt = G2spc.GenHKLf(H,SGData)[0]
                        if iabsnt:  #flag space gp. absences
                            if Super: 
                                if not ref[2+Super]:
                                    ref[3+Super] = 0
                                else:
                                    ref[3+Super] = 1    #twin id
                            else:
                                ref[3] = 0
                    UseList[histoName] = SetDefaultDData(reflData['Type'],histoName)
                elif histoName in PWDRlist:
                    Id = GetGPXtreeItemId(self,self.root,histoName)
                    refList = self.GPXtree.GetItemPyData(
                        GetGPXtreeItemId(self,Id,'Reflection Lists'))
                    refList[generalData['Name']] = {}
                    UseList[histoName] = SetDefaultDData('PWDR',histoName,NShkl=NShkl,NDij=NDij)
                else:
                    raise Exception('Unexpected histogram '+histoName)
        wx.EndBusyCursor()
        self.EnableRefineCommand()
        
        return # success
        
    def _Add_ImportMenu_Image(self,parent):
        '''configure the Import Image menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu, 'Image','Import image file')
        for reader in self.ImportImageReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportImage, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,'guess format from file','Import image data, use file to try to determine format')
        self.Bind(wx.EVT_MENU, self.OnImportImage, id=item.GetId())
        
    def OnImportImage(self,event):
        '''Called in response to an Import/Image/... menu item
        to read an image from a file. Like all the other imports,
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        A reader object is filled each time an image is read. 
        '''
        self.CheckNotebook()
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        rdlist = self.OnImportGeneric(reqrdr,self.ImportImageReaderlist,
            'image',multiple=True,Preview=False,load2Tree=True)
        if rdlist: 
            self.GPXtree.SelectItem(GetGPXtreeItemId(self,self.Image,'Image Controls'))             #show last image to have beeen read
                    
    def _Add_ImportMenu_Sfact(self,parent):
        '''configure the Import Structure Factor menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'Structure Factor','Import Structure Factor data')
        for reader in self.ImportSfactReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)               
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportSfact, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,'guess format from file','Import Structure Factor, use file to try to determine format')
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
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())
        rdlist = self.OnImportGeneric(reqrdr,self.ImportSfactReaderlist,
            'Structure Factor',multiple=True)
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
            HistName = 'HKLF '+G2obj.StripUnicode(HistName,'_')
            # make new histogram names unique
            if len(rd.Banks):
                for Bank in rd.Banks:
                    valuesdict = {'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),}
                    HistName = G2obj.MakeUniqueLabel(HistName,HKLFlist)
                    print (u'Read structure factor table '+HistName+u' from file '+self.lastimport)
                    Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
                    if not Bank['RefDict'].get('FF'):
                        Bank['RefDict']['FF'] = {}
                    self.GPXtree.SetItemPyData(Id,[valuesdict,Bank['RefDict']])
                    Sub = self.GPXtree.AppendItem(Id,text='Instrument Parameters')
                    self.GPXtree.SetItemPyData(Sub,copy.copy(rd.Parameters))
                    self.GPXtree.SetItemPyData(
                        self.GPXtree.AppendItem(Id,text='Reflection List'),{})  #dummy entry for GUI use
                    newHistList.append(HistName)
            else:
                valuesdict = {'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),}
                HistName = G2obj.MakeUniqueLabel(HistName,HKLFlist)
                print (u'Read structure factor table '+HistName+u' from file '+self.lastimport)
                if not rd.RefDict.get('FF'):
                    rd.RefDict['FF'] = {}
                Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
                self.GPXtree.SetItemPyData(Id,[valuesdict,rd.RefDict])
                Sub = self.GPXtree.AppendItem(Id,text='Instrument Parameters')
                self.GPXtree.SetItemPyData(Sub,rd.Parameters)
                self.GPXtree.SetItemPyData(
                    self.GPXtree.AppendItem(Id,text='Reflection List'),{})  #dummy entry for GUI use
                newHistList.append(HistName)
                
            self.GPXtree.SelectItem(Id)
            self.GPXtree.Expand(Id)
            self.Sngl = True

        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = list(usedHistograms.keys()) # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\nsingle crystal dataset(s) to:'
        for Name in newHistList:
            header += '\n  '+str(Name)
        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        wx.BeginBusyCursor()
        item, cookie = self.GPXtree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            data = self.GPXtree.GetItemPyData(item)
            item, cookie = self.GPXtree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            Super = generalData.get('Super',0)
            SuperVec = []
            if Super:
                SuperVec = np.array(generalData['SuperVec'][0])
            UseList = data['Histograms']
            for histoName in newHistList:
                #redo UpdateHKLFdata(histoName) here:
                Id = GetGPXtreeItemId(self,self.root,histoName)
                refDict,reflData = self.GPXtree.GetItemPyData(Id)
                UseList[histoName] = SetDefaultDData(reflData['Type'],histoName)
                G,g = G2lat.cell2Gmat(generalData['Cell'][1:7])
                if 'TwMax' in reflData:     #nonmerohedral twins present
                    UseList[histoName]['Twins'] = []
                    for iT in range(reflData['TwMax'][0]+1):
                        if iT in reflData['TwMax'][1]:
                            UseList[histoName]['Twins'].append([False,0.0])
                        else:
                            UseList[histoName]['Twins'].append([np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,reflData['TwMax'][0]]])
                else:   #no nonmerohedral twins
                    UseList[histoName]['Twins'] = [[np.array([[1,0,0],[0,1,0],[0,0,1]]),[1.0,False,0]],]
                for iref,ref in enumerate(reflData['RefList']):
                    hkl = ref[:3]
                    if Super:
                        H = list(hkl+SuperVec*ref[3])
                    else:
                        H = hkl
                    ref[4+Super] = np.sqrt(1./G2lat.calc_rDsq2(H,G))
                    iabsnt,mul,Uniq,phi = G2spc.GenHKLf(H,SGData)
                    if iabsnt:  #flag space gp. absences
                        if Super: 
                            if not ref[2+Super]:
                                ref[3+Super] = 0
                            else:
                                ref[3+Super] = 1    #twin id?
                        else:
                            ref[3] = 0
        wx.EndBusyCursor()
        self.EnableRefineCommand()        
        return # success

    def _Add_ImportMenu_powder(self,parent):
        '''configure the Powder Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'Powder Data','Import Powder data')
        for reader in self.ImportPowderReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,'guess format from file','Import powder data, use file to try to determine format')
        self.Bind(wx.EVT_MENU, self.OnImportPowder, id=item.GetId())
        submenu.AppendSeparator()
        item = submenu.Append(wx.ID_ANY,'Simulate a dataset','Create a powder data set entry that will be simulated')
        self.Bind(wx.EVT_MENU, self.OnDummyPowder, id=item.GetId())
        item = submenu.Append(wx.ID_ANY,'Auto Import','Import data files as found')
        def OnAutoImport(event):
            G2G.AutoLoadFiles(self,FileTyp='pwd')
        self.Bind(wx.EVT_MENU, OnAutoImport, id=item.GetId())
        
        item = submenu.Append(wx.ID_ANY,'Fit instr. profile from fundamental parms...','')
        self.Bind(wx.EVT_MENU, self.OnPowderFPA, id=item.GetId())
        
    def OpenPowderInstprm(self,instfile):
        '''Read a GSAS-II (new) instrument parameter file

        :param str instfile: name of instrument parameter file

        '''
        File = open(instfile,'r')
        lines = File.readlines()
        File.close()
        return lines        
           
    def ReadPowderInstprm(self,instLines,bank,databanks,rd):
        '''Read lines from a GSAS-II (new) instrument parameter file
        similar to G2pwdGUI.OnLoad
        If instprm file has multiple banks each with header #Bank n: ..., this 
        finds matching bank no. to load - problem with nonmatches?

        Note that this routine performs a similar role to :func:`GSASIIfiles.ReadPowderInstprm`,
        but this will call a GUI routine for selection when needed. TODO: refactor to combine

        :param list instLines: strings from GSAS-II parameter file; can be concatenated with ';'
        :param int  bank: bank number to check when instprm file has '#BANK n:...' strings
            when bank = n then use parameters; otherwise skip that set. Ignored if BANK n: 
            not present. NB: this kind of instprm file made by a Save all profile command in Instrument Parameters 
        :return dict: Inst  instrument parameter dict if OK, or
                str: Error message if failed    
        '''
        if 'GSAS-II' not in instLines[0]: # not a valid file
            return 'Not a valid GSAS-II instprm file'
        newItems = []
        newVals = []
        Found = False
        il = 0
        if bank is None: # no bank was specified in the input file, is more than one present in file?
            banklist = set([])
            for S in instLines:
                if S[0] == '#' and 'Bank' in S:
                    banklist.add(int(S.split(':')[0].split()[1]))
            if len(banklist) > 1: # yes, the user must make a selection
                choices = [str(i) for i in banklist]
                bank = int(G2G.ItemSelector(choices,self,multiple=False))
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
                                return 'Bank %d not found in instprm file'%(bank)
                            S = instLines[il]
                        continue
                else:   #a non #Bank file
                    il += 1
                    S = instLines[il]
            Found = True
            if '"""' in S:
                delim = '"""' 
            elif "'''" in S:
                delim = "'''"
            else:
                S = S.replace(' ','')
                SS = S.strip().split(';')
                for s in SS:
                    [item,val] = s.split(':',1)
                    newItems.append(item)
                    try:
                        newVals.append(float(val))
                    except ValueError:
                        newVals.append(val)
                il += 1
                continue
            # read multiline values, delimited by ''' or """
            item,val = S.strip().split(':',1)
            val = val.replace(delim,'').rstrip()
            val += '\n'
            while True:
                il += 1
                if il >= len(instLines): break
                S = instLines[il]
                if delim in S:
                    val += S.replace(delim,'').rstrip()
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
        return [G2fil.makeInstDict(newItems,newVals,len(newVals)*[False,]),{}]
        
    def ReadPowderIparm(self,instfile,bank,databanks,rd):
        '''Read a GSAS (old) instrument parameter file

        :param str instfile: name of instrument parameter file
        :param int bank: the bank number read in the raw data file
        :param int databanks: the number of banks in the raw data file.
          If the number of banks in the data and instrument parameter files
          agree, then the sets of banks are assumed to match up and bank
          is used to select the instrument parameter file. If not and not TOF, 
          the user is asked to make a selection.
        :param obj rd: the raw data (histogram) data object. This
          sets rd.instbank.

        '''
        if not os.path.exists(instfile): # no such file
            return {}
        fp = 0
        try:
            fp = open(instfile,'r')
            Iparm = {}
            for S in fp:
                if '#' in S[0]:
                    continue
                Iparm[S[:12]] = S[12:-1]
        except IOError:
            print(u'Error reading file: {}'.format(instfile))
        if fp:        
            fp.close()

        ibanks = int(Iparm.get('INS   BANK  ','1').strip())
        if ibanks == 1: # there is only one bank here, return it
            rd.instbank = 1
            rd.powderentry[2] = 1
            return Iparm
        if 'PNT' in Iparm['INS   HTYPE ']:      #allow mismatch between banks in data  iparm file for TOF
            rd.instbank = bank
        elif ibanks != databanks or bank is None:
            choices = []
            for i in range(1,1+ibanks):
                choices.append('Bank '+str(i))
            bank = 1 + G2G.BlockSelector(
                choices, self,
                title=u'Select an instrument parameter bank for '+
                os.path.split(rd.powderentry[0])[1]+u' BANK '+str(bank)+
                u'\nOr use Cancel to select from the default parameter sets',
                header='Block Selector')
        if bank is None: return {}
        # pull out requested bank # bank from the data, and change the bank to 1
        IparmS = {}
        for key in Iparm:
            if 'INS' in key[:3]:    #skip around rubbish lines in some old iparm files
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
          and the second used for TOF lookup tables for profile coeff.
        '''
                
        def GetDefaultParms(self,rd):
            '''Solicits from user a default set of parameters & returns Inst parm dict
            param: self: refers to the GSASII main class
            param: rd: importer data structure
            returns: dict: Instrument parameter dictionary
            '''       
            sind = lambda x: math.sin(x*math.pi/180.)
            tand = lambda x: math.tan(x*math.pi/180.)
            while True: # loop until we get a choice
                choices = []
                head = 'Select from default instrument parameters for '+rd.idstring
    
                for l in dI.defaultIparm_lbl:
                    choices.append('Defaults for '+l)
                res = G2G.BlockSelector(choices,ParentFrame=self,title=head,
                    header='Select default inst parms',useCancel=True)
                if res is None: return None
                rd.instfile = ''
                if 'lab data' in choices[res]:
                    rd.Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
                        'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
                else:
                    rd.Sample.update({'Type':'Debye-Scherrer','Absorption':[0.,False],'DisplaceX':[0.,False],
                        'DisplaceY':[0.,False]})
                if 'Generic' in choices[res]:
                    dlg = G2G.MultiDataDialog(self,title='Generic TOF detector bank',
                        prompts=['Total FP','2-theta',],values=[25.0,150.,],
                            limits=[[6.,200.],[5.,175.],],formats=['%6.2f','%6.1f',])
                    if dlg.ShowModal() == wx.ID_OK: #strictly empirical approx.
                        FP,tth = dlg.GetValues()
                        difC = 505.632*FP*sind(tth/2.)
                        sig1 = 50.+2.5e-6*(difC/tand(tth/2.))**2
                        bet1 = .00226+7.76e+11/difC**4
                        rd.instmsg = 'default: '+dI.defaultIparm_lbl[res]
                        Inst = self.ReadPowderInstprm(dI.defaultIparms[res],bank,numbanks,rd)
                        Inst[0]['difC'] = [difC,difC,0]
                        Inst[0]['sig-1'] = [sig1,sig1,0]
                        Inst[0]['beta-1'] = [bet1,bet1,0]
                        return Inst    #this is [Inst1,Inst2] a pair of dicts
                    dlg.Destroy()
                else:
                    rd.instmsg = 'default: '+dI.defaultIparm_lbl[res]
                    inst1,inst2 = self.ReadPowderInstprm(dI.defaultIparms[res],bank,numbanks,rd)
                    if rd.instdict.get('wave'):
                        inst1['Lam'][0] = rd.instdict.get('wave')
                        inst1['Lam'][1] = rd.instdict.get('wave')
                    return [inst1,inst2]

        # stuff we might need from the reader
        filename = rd.powderentry[0]
        bank = rd.powderentry[2]
        numbanks = rd.numbanks
        #1st priority: is there an instrument parameter file matching the current file
        # with extension .instprm, .prm, .inst, or .ins? If so read it
        basename = os.path.splitext(filename)[0]
        for ext in '.prm','.inst','.ins','.instprm':
            if self.zipfile:
                instfile = G2IO.ExtractFileFromZip(self.zipfile,
                    selection=os.path.split(basename + ext)[1],parent=self)
                if instfile == None:
                    continue
            else:
                instfile = basename + ext
            if not os.path.exists(instfile):
                continue
            if 'instprm' in instfile:
                Lines = self.OpenPowderInstprm(instfile)
                instParmList = self.ReadPowderInstprm(Lines,bank,numbanks,rd)    #this is [Inst1,Inst2] a pair of dicts
                if 'list' in str(type(instParmList)):
                    rd.instfile = instfile
                    rd.instmsg = 'GSAS-II file '+instfile
                    return instParmList
                else:
                    #print 'debug: open/read failed',instfile
                    pass # fail silently
            else:
                Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                if Iparm:
                    #print 'debug: success'
                    rd.instfile = instfile
                    rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                    return G2fil.SetPowderInstParms(Iparm,rd)
                else:
                    #print 'debug: open/read failed',instfile
                    pass # fail silently

        #2nd priority: is there an instrument parameter file defined for the current data set?
        # or if this is a read on a set of set of files, use the last one again
        #if rd.instparm as found in data file header or (lastdatafile == filename and lastIparmfile):
        if rd.instparm or lastIparmfile:
            if rd.instparm:
                instfile = os.path.join(os.path.split(filename)[0],rd.instparm)
            else:
                # for multiple reads of one data file, reuse the inst parm file
                instfile = lastIparmfile
#            if self.zipfile:
#                instfile = G2IO.ExtractFileFromZip(self.zipfile,
#                    selection=os.path.split(instfile)[1],parent=self)
            if instfile != None and os.path.exists(instfile):
                #print 'debug: try read',instfile
                if 'instprm' in instfile:   #GSAS-II file must have .instprm as extension
                    Lines = self.OpenPowderInstprm(instfile)
                    if Lines is not None:
                        instParmList = self.ReadPowderInstprm(Lines,bank,numbanks,rd)   #this is [Inst1,Inst2] a pair of dicts
                else:   #old GSAS style iparm file - could be named anything!
                    Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                    if Iparm:
                        #print 'debug: success'
                        rd.instfile = instfile
                        rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                        instParmList = G2fil.SetPowderInstParms(Iparm,rd)     #this is [Inst1,Inst2] a pair of dicts
                if 'list' in str(type(instParmList)):   #record stuff & return stuff
                    rd.instfile = instfile
                    rd.instmsg = 'GSAS-II file '+instfile
                    return instParmList
                else:   #bad iparms - try default
                    rd.instmsg = instParmList   #an error message
                    return GetDefaultParms(self,rd)
            else:
                self.ErrorDialog('Open Error',u'Error opening instrument parameter file '
                    +u'{} requested by file '.format(instfile,filename))
        #Finally - ask user for Instrument parametrs file - seems it can't be in a zip file
        while True: # loop until we get a file that works or we get a cancel
            instfile = ''
            pth = G2G.GetImportPath(self)
            if not pth: pth = '.'
            extOrd = [0,1]
            if GSASIIpath.GetConfigValue('Instprm_default',False):
                extOrd = [1,0]
            extList = ['GSAS iparm file (*.prm,*.inst,*.ins)|*.prm;*.inst;*.ins|','GSAS-II iparm file (*.instprm)|*.instprm|']
            dlg = wx.FileDialog(self,
                u'Choose inst. param file for "'+rd.idstring+u'" (or Cancel for default)',
                pth, '',extList[extOrd[0]]+extList[extOrd[1]]+'All files (*.*)|*.*', wx.FD_OPEN)
            if os.path.exists(lastIparmfile):
                dlg.SetFilename(lastIparmfile)
            if dlg.ShowModal() == wx.ID_OK:
                instfile = dlg.GetPath()
            dlg.Destroy()
            if not instfile: 
                return GetDefaultParms(self,rd) #on Cancel/break
            if 'instprm' in instfile:
                Lines = self.OpenPowderInstprm(instfile)
                if Lines is not None:
                    instParmList = self.ReadPowderInstprm(Lines,bank,numbanks,rd)    #this is [Inst1,Inst2] a pair of dicts
                if 'list' in str(type(instParmList)):
                    rd.instfile = instfile
                    rd.instmsg = 'GSAS-II file '+instfile
                    return instParmList
                else:
                    rd.instmsg = instParmList   #an error message
                    return GetDefaultParms(self,rd)
            else:
                Iparm = self.ReadPowderIparm(instfile,bank,numbanks,rd)
                if Iparm:
                    #print 'debug: success with',instfile
                    rd.instfile = instfile
                    rd.instmsg = instfile + ' bank ' + str(rd.instbank)
                    return G2fil.SetPowderInstParms(Iparm,rd)
                else:
                    self.ErrorDialog('Read Error',
                                     u'Error opening/reading file {}'.format(instfile))
    def EnableRefineCommand(self):
        haveData = False
        # check for phases connected to histograms
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if sub: 
            item, cookie = self.GPXtree.GetFirstChild(sub)
            while item: # loop over phases
                data = self.GPXtree.GetItemPyData(item)
                item, cookie = self.GPXtree.GetNextChild(sub, cookie)
                UseList = data['Histograms']
                if UseList: haveData = True
            if haveData:
                self.dataWindow.DataMenu.Enable(G2G.wxID_DATADELETE,True)
                for item in self.Refine: item.Enable(True)
        else:
            self.dataWindow.DataMenu.Enable(G2G.wxID_DATADELETE,False)
            for item in self.Refine: item.Enable(False)

        
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
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
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
#        lastVals = []
        self.EnablePlot = False
        Iparms = {}
        for rd in rdlist:
            if 'Instrument Parameters' in rd.pwdparms:
                Iparm1,Iparm2 = rd.pwdparms['Instrument Parameters']
            elif Iparms and not lastIparmfile:
                Iparm1,Iparm2 = Iparms
            else:
                # get instrument parameters for each dataset, unless already set
#                if lastIparmfile:  # is this histogram like previous?
#                    if lastVals != (rd.powderdata[0].min(),rd.powderdata[0].max(),len(rd.powderdata[0])):
#                        lastIparmfile = ''
                Iparms = self.GetPowderIparm(rd, Iparm, lastIparmfile, lastdatafile)
                if not Iparms:  #may have bailed out
                    Id = 0
                    continue
                Iparm1,Iparm2 = Iparms
                if rd.repeat_instparm: 
                    lastIparmfile = rd.instfile
                else:
                    Iparms = {}
#                    lastVals = (rd.powderdata[0].min(),rd.powderdata[0].max(),len(rd.powderdata[0]))
                # override any keys in read instrument parameters with ones set in import
                for key in Iparm1:  
                    if key in rd.instdict:
                        Iparm1[key] = rd.instdict[key]
            lastdatafile = rd.powderentry[0]
            if 'phoenix' in wx.version():
                HistName = 'PWDR '+rd.idstring
            else:
                HistName = 'PWDR '+G2obj.StripUnicode(rd.idstring,'_')
            # make new histogram names unique
            if HistName in PWDRlist:
                dlg = wx.MessageDialog(self,'Skip %s?'%(HistName),'Duplicate data name',wx.YES_NO)
                try:
                    if dlg.ShowModal() == wx.ID_YES:
                        Id = 0
                        continue
                finally:
                    dlg.Destroy()
            HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)
            try:
                print('Read powder data '+HistName+ 
                    ' from file '+G2obj.StripUnicode(rd.readfilename) +
                    ' (format: '+ rd.formatName + 
                    '). Inst parameters from '+G2obj.StripUnicode(rd.instmsg))
            except:
                print('Read powder data') 
            # data are read, now store them in the tree
            Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
            if 'T' in Iparm1['Type'][0]:
                if not rd.clockWd and rd.GSAS:
                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
                cw = np.diff(rd.powderdata[0])
                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
                if rd.GSAS:     #NB: old GSAS wanted intensities*CW even if normalized!
                    npts = min(len(rd.powderdata[0]),len(rd.powderdata[1]),len(cw))
                    rd.powderdata[1] = rd.powderdata[1][:npts]/cw[:npts]
                    rd.powderdata[2] = rd.powderdata[2][:npts]*cw[:npts]**2  #1/var=w at this point
                else:       #NB: from topas/fullprof type files
                    rd.powderdata[1] = rd.powderdata[1][:-1]
                    rd.powderdata[2] = rd.powderdata[2][:-1]
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
                rd.powderdata[1] = np.where(np.isinf(rd.powderdata[1]),0.,rd.powderdata[1])
                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])
                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])
                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])
            Ymin = np.min(rd.powderdata[1])                 
            Ymax = np.max(rd.powderdata[1])                 
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxsize),
                'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-.1*Ymax,'refDelt':0.1*Ymax,
                'Yminmax':[Ymin,Ymax]
                }
            # apply user-supplied corrections to powder data
            if 'CorrectionCode' in Iparm1:
                print('Applying corrections from instprm file')
                corr = Iparm1['CorrectionCode'][0]
                try:
                    exec(corr)
                    print('done')
                except Exception as err:
                    print(u'error: {}'.format(err))
                    print('with commands -------------------')
                    print(corr)
                    print('---------------------------------')
                finally:
                    del Iparm1['CorrectionCode']
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.GPXtree.SetItemPyData(Id,[valuesdict,rd.powderdata])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Comments'),
                rd.comments)
            Tmin = min(rd.powderdata[0])
            Tmax = max(rd.powderdata[0])
            Tmin1 = Tmin
            if 'NT' in Iparm1['Type'][0] and G2lat.Pos2dsp(Iparm1,Tmin) < 0.4:                
                Tmin1 = G2lat.Dsp2pos(Iparm1,0.4)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Limits'),
                rd.pwdparms.get('Limits',[(Tmin,Tmax),[Tmin1,Tmax]])
                )
            self.PatternId = GetGPXtreeItemId(self,Id,'Limits')
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Background'),
                rd.pwdparms.get('Background',
                    [['chebyschev-1',True,3,1.0,0.0,0.0],{'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                    )
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Peak List')
                ,{'peaks':[],'sigDict':{}})
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Index Peak List'),
                [[],[]])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Unit Cells List'),
                [])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Reflection Lists'),
                {})
            # if any Control values have been set, move them into tree
            Controls = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,self.root, 'Controls'))
            Controls.update(rd.Controls)
            newHistList.append(HistName)
            rd.repeat_instparm = False  #clear the iparm reuse flag
        else:
            self.EnablePlot = True
            if Id:
                self.GPXtree.Expand(Id)
                self.GPXtree.SelectItem(Id)

        if not newHistList: return # somehow, no new histograms
        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = list(usedHistograms.keys()) # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to link\nto the newly-read data:'
        for Name in newHistList:
            header += '\n  '+str(Name)

        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        item, cookie = self.GPXtree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            data = self.GPXtree.GetItemPyData(item)
            item, cookie = self.GPXtree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            for histoName in newHistList:
                UseList[histoName] = SetDefaultDData('PWDR',histoName,NShkl=NShkl,NDij=NDij)
                Id = GetGPXtreeItemId(self,self.root,histoName)
                refList = self.GPXtree.GetItemPyData(
                    GetGPXtreeItemId(self,Id,'Reflection Lists'))
                refList[generalData['Name']] = []
        self.EnableRefineCommand()
        return # success

    def OnDummyPowder(self,event):
        '''Called in response to Import/Powder Data/Simulate menu item
        to create a Dummy powder diffraction data set. 

        Reads an instrument parameter file and then gets input from the user
        '''
        # get a list of existing histograms
        PWDRlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # Initialize a base class reader
        rd = G2obj.ImportPowderData(
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
        lastdatafile = ''
        self.zipfile = None
        # get instrument parameters for it
        Iparm = self.GetPowderIparm(rd, Iparm, '', lastdatafile)
        if Iparm is None:
            return
        Iparm1, Iparm2 = Iparm
        if 'T' in Iparm1['Type'][0]:
            rd.idstring = ' TOF neutron simulation'
            simType = 'TOF'
        else:
            # need to get name, 2theta start, end, step
            rd.idstring = ' CW'
            simType = 'CW'
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
            if 'TOF' in rd.idstring:
                names = ('dataset name', 'start TOF(ms)', 'end TOF(ms)', 'DT/T')
                inp = [rd.idstring, 10.,80.,0.0005] # see names for what's what
                dlg = G2G.ScrolledMultiEditor(
                    self,[inp] * len(inp),range(len(inp)),names,
                    header='Enter simulation name and range',
                    minvals=(None,.5,1.0,0.0001),
                    maxvals=(None,200.,200.,.001),
                    sizevals=((225,-1),)
                    )
            else:
                names = ('dataset name', 'start angle', 'end angle', 'step size')
                if not wave or wave < 1.0:
                    inp = [rd.idstring, 10.,40.,0.005] # see names for what's what
                else:
                    inp = [rd.idstring, 10.,80.,0.01] # see names for what's what
                dlg = G2G.ScrolledMultiEditor(
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
            if 'TOF' in rd.idstring:
                N = (np.log(end)-np.log(start))/step
                x = np.exp((np.arange(0,N))*step+np.log(start*1000.))
                N = len(x)
            else:            
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
        HistName = inp[0]
        HistName = 'PWDR '+HistName
        HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)  # make new histogram names unique
        Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
        Ymin = np.min(rd.powderdata[1])
        Ymax = np.max(rd.powderdata[1])
        valuesdict = {
            'wtFactor':1.0,
            'Dummy':True,'simType':simType,
            'ranId':ran.randint(0,sys.maxsize),
            'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-.1*Ymax,'refDelt':0.1*Ymax,
            'Yminmax':[Ymin,Ymax]
            }
        self.GPXtree.SetItemPyData(Id,[valuesdict,rd.powderdata])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Comments'),
            rd.comments)
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        self.PatternId = GetGPXtreeItemId(self,Id,'Limits')
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Background'),
            [['chebyschev-1',True,3,1.0,0.0,0.0],
             {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Instrument Parameters'),
            [Iparm1,Iparm2])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Sample Parameters'),
            rd.Sample)
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Peak List')
            ,{'peaks':[],'sigDict':{}})
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Index Peak List'),
            [[],[]])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Unit Cells List'),
            [])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Reflection Lists'),
            {})
        self.GPXtree.Expand(Id)
        self.GPXtree.SelectItem(Id)
        print(u'Added simulation powder data {}'.format(HistName)+
              ' with parameters from {}'.format(rd.instmsg))

        # make a list of phase names
        phaseRIdList,usedHistograms = self.GetPhaseInfofromTree()
        phaseNameList = list(usedHistograms.keys()) # phase names in use
        if not phaseNameList: return # no phases yet, nothing to do
        header = 'Select phase(s) to add the new\npowder simulation (dummy) dataset to:'
        result = G2G.ItemSelector(phaseNameList,self,header,header='Add to phase(s)',multiple=True)
        if not result: return
        # connect new phases to histograms
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if not sub:
            raise Exception('ERROR -- why are there no phases here?')
        item, cookie = self.GPXtree.GetFirstChild(sub)
        iph = -1
        while item: # loop over (new) phases
            iph += 1
            data = self.GPXtree.GetItemPyData(item)
            item, cookie = self.GPXtree.GetNextChild(sub, cookie)
            if iph not in result: continue
            generalData = data['General']
            SGData = generalData['SGData']
            UseList = data['Histograms']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            UseList[HistName] = SetDefaultDData('PWDR',HistName,NShkl=NShkl,NDij=NDij)
            Id = GetGPXtreeItemId(self,self.root,HistName)
            refList = self.GPXtree.GetItemPyData(
                GetGPXtreeItemId(self,Id,'Reflection Lists'))
            refList[generalData['Name']] = []
        cId = GetGPXtreeItemId(self,self.root, 'Controls')
        Controls = self.GPXtree.GetItemPyData(cId)
        Controls['max cyc'] = 0
        self.EnableRefineCommand()
        return # success

    def AddSimulatedPowder(self,ttArr,intArr,HistName,Lam1,Lam2):
        '''Create a PWDR entry for a computed powder pattern
        '''
        # get a list of existing histograms
        PWDRlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('PWDR ') and name not in PWDRlist:
                    PWDRlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # Initialize a base class reader
        rd = G2obj.ImportPowderData(
            extensionlist=tuple(),
            strictExtension=False,
            formatName = 'FPA Simulated dataset',
            longFormatName = 'Fundamental Parameters simulated pattern')
        rd.powderentry[0] = '' # no filename
        # #self.powderentry[1] = pos # bank offset (N/A here)
        rd.powderentry[2] = 1 # only one bank
        rd.comments.append('This is a powder pattern simulated with Fundamental Parameters')
        self.CheckNotebook()
        #self.zipfile = None
        # get instrument parameters for it
        rd.Sample.update({'Type':'Bragg-Brentano','Shift':[0.,False],'Transparency':[0.,False],
            'SurfRoughA':[0.,False],'SurfRoughB':[0.,False]})
        Iparm1, Iparm2 = G2fil.ReadPowderInstprm(dI.defaultIparms[0],1,1,rd)
        rd.idstring = ' CW'
        simType = 'CW'
        # set wavelength
        if Lam2:
            Iparm1['Lam1'][0] = Lam1
            Iparm1['Lam2'][0] = Lam2
            Iparm1['Lam1'][1] = Lam1
            Iparm1['Lam2'][1] = Lam2
        else:
            Iparm1['Lam'] = Iparm1['Lam1']
            del Iparm1['Lam1'],Iparm1['Lam2']
            Iparm1['Lam'][0] = Lam1
            Iparm1['Lam'][1] = Lam1

        rd.powderdata = [
            np.array(ttArr), # x-axis values
            np.array(intArr), # powder pattern intensities
            np.ones_like(ttArr), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(intArr), # calc. intensities (zero)
            np.zeros_like(ttArr), # calc. background (zero)
            np.zeros_like(ttArr), # obs-calc profiles
            ]
        Tmin = rd.powderdata[0][0]
        Tmax = rd.powderdata[0][-1]
        # data are read, now store them in the tree
        HistName = 'PWDR '+HistName
        HistName = G2obj.MakeUniqueLabel(HistName,PWDRlist)  # make new histogram names unique
        Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
        Ymin = np.min(rd.powderdata[1])
        Ymax = np.max(rd.powderdata[1])
        valuesdict = {
            'wtFactor':1.0,
            'Dummy':True,'simType':simType,
            'ranId':ran.randint(0,sys.maxsize),
            'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-.1*Ymax,'refDelt':0.1*Ymax,
            'Yminmax':[Ymin,Ymax]
            }
        self.GPXtree.SetItemPyData(Id,[valuesdict,rd.powderdata])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Comments'),
            rd.comments)
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        self.PatternId = GetGPXtreeItemId(self,Id,'Limits')
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Background'),
            [['chebyschev-1',True,3,1.0,0.0,0.0],
             {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Instrument Parameters'),
            [Iparm1,Iparm2])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Sample Parameters'),
            rd.Sample)
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Peak List')
            ,{'peaks':[],'sigDict':{}})
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Index Peak List'),
            [[],[]])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Unit Cells List'),
            [])
        self.GPXtree.SetItemPyData(
            self.GPXtree.AppendItem(Id,text='Reflection Lists'),
            {})
        self.GPXtree.Expand(Id)
        self.GPXtree.SelectItem(Id)
        print(u'Added simulation powder data {}'.format(HistName))
        return Id
    
    def OnPreferences(self,event):
        'Edit the GSAS-II configuration variables'
        dlg = G2G.SelectConfigSetting(self)
        dlg.ShowModal() == wx.ID_OK
        dlg.Destroy()

    def EditProxyInfo(self,event):
        '''Edit the proxy information used by subversion
        '''
        h,p,e = host,port,etc = GSASIIpath.getsvnProxy()
        labels = ['Proxy address','proxy port']
        values = [host,port]
        i = 1
        for item in etc:
            i += 1
            labels.append('extra svn arg #'+str(i))
            values.append(item)
        msg = '''This dialog allows customization of the subversion (svn) 
        command. If a proxy server is needed, the address/host and port 
        can be added supplied here. This will generate command-line options 

        --config-option servers:global:http-proxy-host=*host* 
        --config-option servers:global:http-proxy-port=*port*

        Additional subversion command line options can be supplied here
        by pressing the '+' button. As examples of options that might be of 
        value, use two extra lines to add:

        --config-dir
        DIR

        to specify an alternate configuration location. 

        Or, use four extra lines to add

        --config-option
        servers:global:http-proxy-username=*account*
        --config-option
        servers:global:http-proxy-password=*password*

        to specify a proxy user name and password.

        Note that strings marked *value* are items that will be configured 
        by the user. See http://svnbook.red-bean.com for more information on
        subversion. 
        '''
        dlg = G2G.MultiStringDialog(self,'Enter proxy values',
                            labels,values,size=300,addRows=True,hlp=msg)
        if dlg.Show():
            values = dlg.GetValues()
            h,p = values[:2]
            e = values[2:]
        dlg.Destroy()
        if h != host or p != port or etc != e:
            localproxy = proxyinfo = os.path.join(
                os.path.expanduser('~/.G2local/'),
                "proxyinfo.txt")
            if not os.path.exists(proxyinfo):
                proxyinfo = os.path.join(GSASIIpath.path2GSAS2,"proxyinfo.txt")
            GSASIIpath.setsvnProxy(h,p,e)
            if not h.strip() and not e:
                if os.path.exists(localproxy): os.remove(localproxy)
                if os.path.exists(proxyinfo): os.remove(proxyinfo)
                return
            try:
                fp = open(proxyinfo,'w')
            except:
                fp = open(localproxy,'w')
                proxyinfo = localproxy
            try:
                fp.write(h.strip()+'\n')
                fp.write(p.strip()+'\n')
                for i in e:
                    if i.strip():
                        fp.write(i.strip()+'\n')
                fp.close()
            except Exception as err:
                print('Error writing file {}:\n{}'.format(proxyinfo,err))
            print('File {} written'.format(proxyinfo))
                
    def _Add_ImportMenu_smallangle(self,parent):
        '''configure the Small Angle Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'Small Angle Data','Import small angle data')
        for reader in self.ImportSmallAngleReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportSmallAngle, id=item.GetId())
        # item = submenu.Append(wx.ID_ANY,
        #     help='Import small angle data, use file to try to determine format',
        #     kind=wx.ITEM_NORMAL,text='guess format from file')
        # self.Bind(wx.EVT_MENU, self.OnImportSmallAngle, id=item.GetId())

    def OnImportSmallAngle(self,event):
        '''Called in response to an Import/Small Angle Data/... menu item
        to read a small angle diffraction data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        '''
        
        def GetSASDIparm(reader):
            parm = reader.instdict
            Iparm = {'Type':[parm['type'],parm['type'],0],'Lam':[parm['wave'],
                parm['wave'],0],'Azimuth':[0.,0.,0]}           
            return Iparm,{}
            
        # get a list of existing histograms
        SASDlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('SASD ') and name not in SASDlist:
                    SASDlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportSmallAngleReaderlist,'Small Angle Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        newHistList = []
        self.EnablePlot = False
        for rd in rdlist:
            HistName = rd.idstring
            HistName = 'SASD '+HistName
            # make new histogram names unique
            HistName = G2obj.MakeUniqueLabel(HistName,SASDlist)
            print ('Read small angle data '+HistName+ \
                ' from file '+self.lastimport)
            # data are read, now store them in the tree
            Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
            Iparm1,Iparm2 = GetSASDIparm(rd)
#            if 'T' in Iparm1['Type'][0]:
#                if not rd.clockWd and rd.GSAS:
#                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
#                cw = np.diff(rd.powderdata[0])
#                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
#                rd.powderdata[1] = rd.powderdata[1][:-1]/cw
#                rd.powderdata[2] = rd.powderdata[2][:-1]*cw**2  #1/var=w at this point
#                if 'Itype' in Iparm2:
#                    Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
#                    Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
#                    rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
#                    YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
#                    rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
#                    var = 1./rd.powderdata[2][Ibeg:Ifin]
#                    var += WYI*rd.powderdata[1]**2
#                    var /= YI**2
#                    rd.powderdata[2] = 1./var
#                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])                                        
            Tmin = min(rd.smallangledata[0])
            Tmax = max(rd.smallangledata[0])
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxsize),
                'Offset':[0.0,0.0],
                }
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.GPXtree.SetItemPyData(Id,[valuesdict,rd.smallangledata])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Comments'),
                rd.comments)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Limits'),
                [(Tmin,Tmax),[Tmin,Tmax]])
            self.PatternId = GetGPXtreeItemId(self,Id,'Limits')
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Substances'),G2pdG.SetDefaultSubstances())
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Models'),G2pdG.SetDefaultSASDModel())
            newHistList.append(HistName)
        else:
            self.EnablePlot = True
            self.GPXtree.Expand(Id)
            self.GPXtree.SelectItem(Id)
            
        if not newHistList: return # somehow, no new histograms
        return # success
        
    def _Add_ImportMenu_reflectometry(self,parent):
        '''configure the reflectometry Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'Reflectometry Data','Import reflectometry data')
        for reader in self.ImportReflectometryReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportReflectometry, id=item.GetId())
        # item = submenu.Append(wx.ID_ANY,
        #     help='Import reflectometry data, use file to try to determine format',
        #     kind=wx.ITEM_NORMAL,text='guess format from file')
        # self.Bind(wx.EVT_MENU, self.OnImportReflectometry, id=item.GetId())
        
    def OnImportReflectometry(self,event):
        '''Called in response to an Import/Reflectometry Data/... menu item
        to read a reflectometry data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        '''
        
        def GetREFDIparm(reader):
            parm = reader.instdict
            Iparm = {'Type':[parm['type'],parm['type'],0],'Lam':[parm['wave'],
                parm['wave'],0],'Azimuth':[0.,0.,0]}           
            return Iparm,{}
            
        # get a list of existing histograms
        REFDlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('REFD ') and name not in REFDlist:
                    REFDlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportReflectometryReaderlist,'Reflectometry Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        newHistList = []
        self.EnablePlot = False
        for rd in rdlist:
            HistName = rd.idstring
            HistName = 'REFD '+HistName
            # make new histogram names unique
            HistName = G2obj.MakeUniqueLabel(HistName,REFDlist)
            print ('Read reflectometry data '+HistName+ \
                ' from file '+self.lastimport)
            # data are read, now store them in the tree
            Id = self.GPXtree.AppendItem(parent=self.root,text=HistName)
            Iparm1,Iparm2 = GetREFDIparm(rd)
#            if 'T' in Iparm1['Type'][0]:
#                if not rd.clockWd and rd.GSAS:
#                    rd.powderdata[0] *= 100.        #put back the CW centideg correction
#                cw = np.diff(rd.powderdata[0])
#                rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
#                rd.powderdata[1] = rd.powderdata[1][:-1]/cw
#                rd.powderdata[2] = rd.powderdata[2][:-1]*cw**2  #1/var=w at this point
#                if 'Itype' in Iparm2:
#                    Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
#                    Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
#                    rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
#                    YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
#                    rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
#                    var = 1./rd.powderdata[2][Ibeg:Ifin]
#                    var += WYI*rd.powderdata[1]**2
#                    var /= YI**2
#                    rd.powderdata[2] = 1./var
#                rd.powderdata[3] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[4] = np.zeros_like(rd.powderdata[0])                                        
#                rd.powderdata[5] = np.zeros_like(rd.powderdata[0])                                        
            Tmin = min(rd.reflectometrydata[0])
            Tmax = max(rd.reflectometrydata[0])
            ifDQ = np.any(rd.reflectometrydata[5])
            valuesdict = {
                'wtFactor':1.0,
                'Dummy':False,
                'ranId':ran.randint(0,sys.maxsize),
                'Offset':[0.0,0.0],
                'ifDQ':ifDQ
                }
            rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
            self.GPXtree.SetItemPyData(Id,[valuesdict,rd.reflectometrydata])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Comments'),
                rd.comments)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Limits'),
                [(Tmin,Tmax),[Tmin,Tmax]])
            self.PatternId = GetGPXtreeItemId(self,Id,'Limits')
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Instrument Parameters'),
                [Iparm1,Iparm2])
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Substances'),G2pdG.SetDefaultSubstances())
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Sample Parameters'),
                rd.Sample)
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='Models'),G2pdG.SetDefaultREFDModel())
            newHistList.append(HistName)
        else:
            self.EnablePlot = True
            self.GPXtree.Expand(Id)
            self.GPXtree.SelectItem(Id)
            
        if not newHistList: return # somehow, no new histograms
        return # success

    def _Add_ImportMenu_PDF(self,parent):
        '''configure the PDF Data menus accord to the readers found in _init_Imports
        '''
        submenu = wx.Menu()
        item = parent.AppendSubMenu(submenu,'PDF G(R) Data','Import PDF G(R) data')
        for reader in self.ImportPDFReaderlist:
            item = submenu.Append(wx.ID_ANY,u'from '+reader.formatName+u' file',reader.longFormatName)
            self.ImportMenuId[item.GetId()] = reader
            self.Bind(wx.EVT_MENU, self.OnImportPDF, id=item.GetId())
        submenu.AppendSeparator()
        item = submenu.Append(wx.ID_ANY,'Auto Import','Import PDF files as found')
        def OnAutoImport(event):
            G2G.AutoLoadFiles(self,FileTyp='gr')
        self.Bind(wx.EVT_MENU, OnAutoImport, id=item.GetId())
        # item = submenu.Append(wx.ID_ANY,
        #     help='Import reflectometry data, use file to try to determine format',
        #     kind=wx.ITEM_NORMAL,text='guess format from file')
        # self.Bind(wx.EVT_MENU, self.OnImportReflectometry, id=item.GetId())
        
    def OnImportPDF(self,event):
        '''Called in response to an Import/PDF G(R) Data/... menu item
        to read a PDF G(R) data set. 
        dict self.ImportMenuId is used to look up the specific
        reader item associated with the menu item, which will be
        None for the last menu item, which is the "guess" option
        where all appropriate formats will be tried.

        '''
        
        # get a list of existing histograms
        PDFlist = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name.startswith('PDF ') and name not in PDFlist:
                    PDFlist.append(name)
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        # look up which format was requested
        reqrdr = self.ImportMenuId.get(event.GetId())  
        rdlist = self.OnImportGeneric(
            reqrdr,self.ImportPDFReaderlist,'PDF G(R) Data',multiple=True)
        if len(rdlist) == 0: return
        self.CheckNotebook()
        newHistList = []
        self.EnablePlot = False
        for rd in rdlist:
            HistName = rd.idstring
            HistName = 'PDF '+HistName
            # make new histogram names unique
            HistName = G2obj.MakeUniqueLabel(HistName,PDFlist)
            print ('Read PDF G(R) data '+HistName+ \
                ' from file '+self.lastimport)
            # data are read, now store them in the tree
            Id = self.GPXtree.AppendItem(self.root,text=HistName)
            Ymin = np.min(rd.pdfdata[1])                 
            Ymax = np.max(rd.pdfdata[1])                 
            valuesdict = {
                'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),
                'Offset':[0.0,0.0],'delOffset':0.02*Ymax,
                'Yminmax':[Ymin,Ymax],
                }
            self.GPXtree.SetItemPyData(
                self.GPXtree.AppendItem(Id,text='PDF Controls'),
                    {'G(R)':[valuesdict,rd.pdfdata,HistName],
                         'diffGRname':'','diffMult':1.0,'Rmax':Ymax,})
            self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='PDF Peaks'),
                {'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]})
        else:
            self.EnablePlot = True
            self.GPXtree.Expand(Id)
            self.GPXtree.SelectItem(Id)
            
        if not newHistList: return # somehow, no new histograms
        return # success
    
    def AddToNotebook(self,text):
        Id =  GetGPXtreeItemId(self,self.root,'Notebook')
        data = self.GPXtree.GetItemPyData(Id)
        data.append('Notebook entry @ %s: %s'%(time.ctime(),text))    
                   
###############################################################################
#Command logging
###############################################################################

    def OnMacroRecordStatus(self,event,setvalue=None):
        '''Called when the record macro menu item is used which toggles the
        value. Alternately a value to be set can be provided. Note that this
        routine is made more complex because on the Mac there are lots of menu
        items (listed in self.MacroStatusList) and this loops over all of them.
        '''
        nextvalue = log.ShowLogStatus() != True
        if setvalue is not None:
            nextvalue = setvalue 
        if nextvalue:
            log.LogOn()
            set2 = True
        else:
            log.LogOff()
            set2 = False
        for menuitem in self.MacroStatusList:
            menuitem.Check(set2)

    def _init_Macro(self):
        '''Define the items in the macro menu.
        '''
        menu = self.MacroMenu
        item = menu.Append(
                help='Start or stop recording of menu actions, etc.', id=wx.ID_ANY,
                kind=wx.ITEM_CHECK,text='Record actions')
        self.MacroStatusList.append(item)
        item.Check(log.ShowLogStatus())
        self.Bind(wx.EVT_MENU, self.OnMacroRecordStatus, item)

        # this may only be of value for development work
        item = menu.Append(
            help='Show logged commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Show log')
        def OnShowLog(event):
            print (70*'=')
            print ('List of logged actions')
            for i,line in enumerate(log.G2logList):
                if line: print ('%d %s'%(i,line))
            print (70*'=')
        self.Bind(wx.EVT_MENU, OnShowLog, item)

        item = menu.Append(
            help='Clear logged commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Clear log')
        def OnClearLog(event): log.G2logList=[None]
        self.Bind(wx.EVT_MENU, OnClearLog, item)
        
        item = menu.Append(
            help='Save logged commands to file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Save log')
        def OnSaveLog(event):
            defnam = os.path.splitext(os.path.split(self.GSASprojectfile)[1])[0]+'.gcmd'
            dlg = wx.FileDialog(self,
                'Choose an file to save past actions', '.', defnam, 
                'GSAS-II cmd output (*.gcmd)|*.gcmd',
                wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    # make sure extension is correct
                    filename = os.path.splitext(filename)[0]+'.gcmd'
                else:
                    filename = None
            finally:
                dlg.Destroy()
            if filename:
                fp = open(filename,'wb')
                fp.write(str(len(log.G2logList)-1)+'\n')
                for item in log.G2logList:
                    if item: cPickle.dump(item,fp)
                fp.close()
        self.Bind(wx.EVT_MENU, OnSaveLog, item)

        item = menu.Append(
            help='Load logged commands from file', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Load log')
        def OnLoadLog(event):
            # this appends. Perhaps we should ask to clear? 
            defnam = os.path.splitext(
                os.path.split(self.GSASprojectfile)[1])[0]+'.gcmd'
            dlg = wx.FileDialog(self,
                'Choose an file to read saved actions', '.', defnam, 
                'GSAS-II cmd output (*.gcmd)|*.gcmd',
                wx.FD_OPEN)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    # make sure extension is correct
                    filename = os.path.splitext(filename)[0]+'.gcmd'
                else:
                    filename = None
            finally:
                dlg.Destroy()
            if filename and os.path.exists(filename):
                fp = open(filename,'rb')
                lines = fp.readline()
                for i in range(int(lines)):
                    log.G2logList.append(cPickle.load(fp))
                fp.close()
        self.Bind(wx.EVT_MENU, OnLoadLog, item)

        item = menu.Append(
            help='Replay saved commands', id=wx.ID_ANY,
            kind=wx.ITEM_NORMAL,text='Replay log')
        self.Bind(wx.EVT_MENU, log.ReplayLog, item)
        
# End of logging ##############################################################

    def _init_Exports(self,menu):
        '''Find exporter routines and add them into menus
        '''
        # set up the top-level menus
        projectmenu = wx.Menu()
        item = menu.AppendSubMenu(projectmenu,'Entire project as','Export entire project')
        self.ExportNonSeq.append([menu,item.Id])
        
        phasemenu = wx.Menu()
        item = menu.AppendSubMenu(phasemenu,'Phase as','Export phase or sometimes phases')

        powdermenu = wx.Menu()
        item = menu.AppendSubMenu(powdermenu,'Powder data as','Export powder diffraction histogram(s)')
        
        sasdmenu = wx.Menu()
        item = menu.AppendSubMenu(sasdmenu,'Small angle data as','Export small angle histogram(s)')

        refdmenu = wx.Menu()
        item = menu.AppendSubMenu(refdmenu,'Reflectometry data as','Export reflectometry histogram(s)')

        singlemenu = wx.Menu()
        item = menu.AppendSubMenu(singlemenu,'Single crystal data as','Export single crystal histogram(s)')

        imagemenu = wx.Menu()
        item = menu.AppendSubMenu(imagemenu,'Image data as','Export powder image(s) data')

        mapmenu = wx.Menu()
        item = menu.AppendSubMenu(mapmenu,'Maps as','Export density map(s)')

        # sequential exports are handled differently; N.B. enabled in testSeqRefineMode
        seqPhasemenu = wx.Menu()
        item = menu.AppendSubMenu(seqPhasemenu,'Sequential phases','Export phases from sequential fit')
        self.ExportSeq.append([menu,item.Id])
        seqHistmenu = wx.Menu()
        item = menu.AppendSubMenu(seqHistmenu,'Sequential histograms','Export histograms from sequential fit')
        self.ExportSeq.append([menu,item.Id])
        
        # find all the exporter files
        if not self.exporterlist: # this only needs to be done once
            self.exporterlist = G2fil.LoadExportRoutines(self)
                
        # Add submenu item(s) for each Exporter by its self-declared type (can be more than one)
        for obj in self.exporterlist:
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
                elif typ == "sasd":
                    submenu = sasdmenu
                elif typ == "refd":
                    submenu = refdmenu
                # elif typ == "pdf":
                #     submenu = pdfmenu
                else:
                    print("Error, unknown type in "+str(obj))
                    break
                item = submenu.Append(wx.ID_ANY,obj.formatName,obj.longFormatName)
                self.Bind(wx.EVT_MENU, obj.Exporter, id=item.GetId())
                self.ExportLookup[item.GetId()] = typ # lookup table for submenu item
            for lbl,submenu in (('Phase',seqPhasemenu),
                        ('Powder',seqHistmenu),
                        ):
                if lbl.lower() in obj.exporttype:
                    try:
                        obj.Writer
                    except AttributeError:
                        continue
                    # define a unique event handler for this menu item
                    def seqMenuItemEventHandler(event,obj=obj,typ=lbl):
                        'This handler has the needed exporter/type embedded'
                        # lookup sequential table
                        Id = GetGPXtreeItemId(self,self.root,'Sequential results')
                        if not Id:
                            print('Error in Seq seqMenuItemEventHandler for ',typ,'without Seq Res table')
                            return
                        data = self.GPXtree.GetItemPyData(Id)
                        G2IO.ExportSequential(self,data,obj,typ)
                        if '2' in platform.python_version_tuple()[0]:
                            if 'mode' in inspect.getargspec(obj.Writer)[0]:
                                item = submenu.Append(wx.ID_ANY,obj.formatName,obj.longFormatName)
                                self.Bind(wx.EVT_MENU, seqMenuItemEventHandler, item)
                        else:
                            if 'mode' in inspect.getfullargspec(obj.Writer)[0]:
                                item = submenu.Append(wx.ID_ANY,obj.formatName,obj.longFormatName)
                                self.Bind(wx.EVT_MENU, seqMenuItemEventHandler, item)
                        #                    self.SeqExportLookup[item.GetId()] = (obj,lbl) # lookup table for submenu item
                        # Bind is in UpdateSeqResults

        item = imagemenu.Append(wx.ID_ANY,'Multiple image controls and masks',
            'Export image controls and masks for multiple images')
        self.Bind(wx.EVT_MENU, self.OnSaveMultipleImg, id=item.GetId())
        #code to debug an Exporter. hard-code the routine below, to allow a reload before use
        # def DebugExport(event):
        #      print 'start reload'
        #      reload(G2IO)
        #      import G2export_pwdr as dev
        #      reload(dev)
        #      dev.ExportPowderFXYE(self).Exporter(event)
        # item = menu.Append(
        #     wx.ID_ANY,kind=wx.ITEM_NORMAL,
        #     help="debug exporter",text="test Export FXYE")
        # self.Bind(wx.EVT_MENU, DebugExport, id=item.GetId())
        # # #self.ExportLookup[item.GetId()] = 'image'
        # self.ExportLookup[item.GetId()] = 'powder'

###############################################################################
# Exporters
###############################################################################
                            
    def _Add_ExportMenuItems(self,parent):
        # item = parent.Append(
        #     help='Select PWDR item to enable',id=wx.ID_ANY,
        #     kind=wx.ITEM_NORMAL,
        #     text='Export Powder Patterns...')
        # self.ExportPattern.append(item)
        # item.Enable(False)
        # self.Bind(wx.EVT_MENU, self.OnExportPatterns, id=item.GetId())

        item = parent.Append(wx.ID_ANY,'Export All Peak Lists...','')
        self.ExportPeakList.append(item)
        item.Enable(True)
        self.Bind(wx.EVT_MENU, self.OnExportPeakList, id=item.GetId())

        item = parent.Append(wx.ID_ANY,'Export HKLs...','')
        self.ExportHKL.append(item)
        self.Bind(wx.EVT_MENU, self.OnExportHKL, id=item.GetId())

        item = parent.Append(wx.ID_ANY,'Export PDF...','Select PDF item to enable')
        self.ExportPDF.append(item)
        item.Enable(False)
        self.Bind(wx.EVT_MENU, self.OnExportPDF, id=item.GetId())
            
    def FillMainMenu(self,menubar,addhelp=True):
        '''Define contents of the main GSAS-II menu for the (main) data tree window.
        For the mac, this is also called for the data item windows as well so that
        the main menu items are data menu as well.
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
        self._Add_ImportMenu_Image(Import)
        self._Add_ImportMenu_Phase(Import)
        self._Add_ImportMenu_powder(Import)
        self._Add_ImportMenu_Sfact(Import)
        self._Add_ImportMenu_smallangle(Import)
        self._Add_ImportMenu_reflectometry(Import)
        self._Add_ImportMenu_PDF(Import)
        
        item = Import.Append(wx.ID_ANY,'Column metadata test','Test Column (.par) metadata import')
        self.Bind(wx.EVT_MENU, self.OnColMetaTest, id=item.GetId())

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
            # if not rd.ContentsValidator(filename):
            #     print 'not validated'
            #     # make a list of used phase ranId's
            # phaseRIdList = []
            # sub = GetGPXtreeItemId(self,self.root,'Phases')
            # if sub:
            #     item, cookie = self.GPXtree.GetFirstChild(sub)
            #     while item:
            #         phaseName = self.GPXtree.GetItemText(item)
            #         ranId = self.GPXtree.GetItemPyData(item).get('ranId')
            #         if ranId: phaseRIdList.append(ranId)
            #         item, cookie = self.GPXtree.GetNextChild(sub, cookie)
            # if rd.Reader(filename,usedRanIdList=phaseRIdList):
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
        if GSASIIpath.GetConfigValue('Enable_logging'):
            self.MacroMenu = wx.Menu(title='')
            menubar.Append(menu=self.MacroMenu, title='Macro')
            self._init_Macro()
        if addhelp:
            HelpMenu=G2G.MyHelp(self,includeTree=True,
                morehelpitems=[('&Tutorials\tCtrl+T','Tutorials'),])
            menubar.Append(menu=HelpMenu,title='&Help')
            
    def _init_ctrls(self, parent):
        try:
            size = GSASIIpath.GetConfigValue('Main_Size')
            if type(size) is tuple:
                pass
            elif type(size) is str:
                size = eval(size)
            else:
                raise Exception
        except:
            size = wx.Size(700,450)
        wx.Frame.__init__(self, name='GSASII', parent=parent,
            size=size,style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II main window')
        self._init_Imports()
        #initialize Menu item objects (these contain lists of menu items that are enabled or disabled)
        self.MakePDF = []
        self.Refine = []
        self.ExportSeq = []
        self.ExportNonSeq = []
        #self.ExportPattern = []
        self.ExportPeakList = []
        self.ExportHKL = []
        self.ExportPDF = []
        self.ExportPhase = []
        self.ExportCIF = []
        #
        self.MacroStatusList = []  # logging
        self.Status = self.CreateStatusBar()
        self.Status.SetFieldsCount(2)
        # Bob: note different ways to display the SplitterWindow. I like the 3d effect on the Mac
        # as it makes the splitter bar a bit easier to "grab" -- this might need to be platform selected.
        #self.mainPanel = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_BORDER|wx.SP_LIVE_UPDATE)
        #self.mainPanel = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_BORDER|wx.SP_LIVE_UPDATE|wx.SP_3DSASH)
        self.mainPanel = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_LIVE_UPDATE|wx.SP_3D)
        self.mainPanel.SetMinimumPaneSize(100)
        self.treePanel = wx.Panel(self.mainPanel, wx.ID_ANY,
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        
        self.dataWindow = G2DataWindow(self.mainPanel)
        dataSizer = wx.BoxSizer(wx.VERTICAL)
        self.dataWindow.SetSizer(dataSizer)
        self.mainPanel.SplitVertically(self.treePanel, self.dataWindow, 200)
        self.Status.SetStatusWidths([200,-1])   # make these match?
        
        G2G.wxID_GPXTREE = wx.NewId()
        treeSizer = wx.BoxSizer(wx.VERTICAL)
        self.treePanel.SetSizer(treeSizer)
        self.GPXtree = G2G.G2TreeCtrl(id=G2G.wxID_GPXTREE,
            parent=self.treePanel, size=self.treePanel.GetClientSize(),style=wx.TR_DEFAULT_STYLE )
        treeSizer.Add(self.GPXtree,1,wx.EXPAND|wx.ALL,0)
        self.GPXtree.Bind(wx.EVT_TREE_SEL_CHANGED,self.OnDataTreeSelChanged)
        self.GPXtree.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK,self.OnDataTreeSelChanged)
        self.GPXtree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
            self.OnGPXtreeItemCollapsed, id=G2G.wxID_GPXTREE)
        self.GPXtree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
            self.OnGPXtreeItemExpanded, id=G2G.wxID_GPXTREE)
        self.GPXtree.Bind(wx.EVT_TREE_DELETE_ITEM,
            self.OnGPXtreeItemDelete, id=G2G.wxID_GPXTREE)
        self.GPXtree.Bind(wx.EVT_TREE_KEY_DOWN,
            self.OnGPXtreeKeyDown, id=G2G.wxID_GPXTREE)
        self.GPXtree.Bind(wx.EVT_TREE_BEGIN_RDRAG,
            self.OnGPXtreeBeginRDrag, id=G2G.wxID_GPXTREE)        
        self.GPXtree.Bind(wx.EVT_TREE_END_DRAG,
            self.OnGPXtreeEndDrag, id=G2G.wxID_GPXTREE)        
        self.root = self.GPXtree.root        

        try:
            size = GSASIIpath.GetConfigValue('Plot_Size')
            if type(size) is tuple:
                pass
            elif type(size) is str:
                size = eval(size)
            else:
                raise Exception
        except:
            size = wx.Size(700,600)                
        self.plotFrame = wx.Frame(None,-1,'GSASII Plots',size=size,
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(self.plotFrame,G2frame=self)
        self.plotFrame.Show()

        for win,var in ((self,'Main_Pos'),(self.plotFrame,'Plot_Pos')):
            try:
                pos = GSASIIpath.GetConfigValue(var)
                if type(pos) is str: pos = eval(pos)
                win.SetPosition(pos)
                if GetDisplay(pos) is None: win.Center()
            except:
                if GSASIIpath.GetConfigValue(var):
                    print('Value for config {} {} is invalid'.format(var,GSASIIpath.GetConfigValue(var)))
                    win.Center()
################################################################################
#### init_vars
################################################################################
    def init_vars(self):
        # initialize default values for GSAS-II "global" variables (saved in main Frame)
        self.oldFocus = None
        self.undofile = ''
        self.TreeItemDelete = False
        self.Weight = False
        self.IfPlot = False
        self.DDShowAll = False
        self.atmSel = ''
        self.PatternId = 0
        self.PickId = 0
        self.PickIdText = None
        self.PeakTable = []
        self.LimitsTable = []
        self.ifX20 = True   #use M20 /= (1+X20) in powder indexing, etc.
        self.HKL = []
        self.Lines = []     # lines used for data limits & excluded regions
        self.MagLines = []  # lines used for plot magnification
        self.itemPicked = None
        self.Interpolate = 'nearest'
        self.ContourColor = GSASIIpath.GetConfigValue('Contour_color','Paired')
        self.VcovColor = 'RdYlGn'
        self.RamaColor = 'Blues'
        self.Projection = 'equal area'
        self.logPlot = False
        self.plusPlot = True
        self.ErrorBars = False
        self.Contour = False
        self.TforYaxis = False
        self.Legend = False
        self.SinglePlot = True
        self.Waterfall = False
        self.selections= None
        self.PDFselections = None
        self.SubBack = False
        self.seqReverse = False
        self.seqLines = True #draw lines between points
        self.plotView = 0
        self.Image = 0
        self.oldImagefile = '' # the name of the last image file read
        self.oldImageTag = None # the name of the tag for multi-image files
        self.PauseIntegration = False
        self.ImageZ = []
        self.Integrate = 0
        self.imageDefault = {}
        self.IntgOutList = [] # list of integration tree item Ids created in G2IO.SaveIntegration
        self.AutointPWDRnames = [] # list of autoint created PWDR tree item names (to be deleted on a reset)
        self.autoIntFrame = None
        self.IntegratedList = [] # list of already integrated IMG tree items
        self.Sngl = False
        self.ifGetRing = False
        self.MaskKey = ''           #trigger for making image masks
        self.MskDelete = False      #trigger for mask delete
        self.StrainKey = ''         #ditto for new strain d-zeros
        self.EnablePlot = True
        self.hist = ''              # selected histogram in Phase/Data tab
        self.dataDisplayPhaseText = ''
        self.lastTreeSetting = [] # used to track the selected Tree item before a refinement
        self.ExpandingAll = False
        self.SeqTblHideList = None
        self.newGPXfile = ''
        self.lastSelectedPhaseTab = None # track the last tab pressed on a phase window
        self.testRBObjSizers = {}   #rigid body sizer datafile contents
        self.RMCchoice = 'RMCProfile' 
        
    def __init__(self, parent):
        self.ExportLookup = {}
        self.exporterlist = []
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
        if 'phoenix' in wx.version():
            self.SetIcon(wx.Icon(img))
        else:
            self.SetIcon(wx.IconFromBitmap(img))
        self.Bind(wx.EVT_CLOSE, self.ExitMain)
        self.GSASprojectfile = ''
        self.dirname = os.path.abspath(os.path.expanduser('~'))       #start in the users home directory by default; may be meaningless
        self.TutorialImportDir = None  # location to read tutorial files, set when a tutorial is viewed
        self.LastImportDir = None # last-used directory where an import was done
        self.LastGPXdir = None    # directory where a GPX file was last read
        self.LastExportDir = None  # the last directory used for exports, if any.
        self.dataDisplay = None
        self.init_vars()
                        
        arg = sys.argv
        if len(arg) > 1 and arg[1]:
            try:
                self.GSASprojectfile = os.path.splitext(arg[1])[0]+u'.gpx'
            except:
                self.GSASprojectfile = os.path.splitext(arg[1])[0]+'.gpx'
            self.dirname = os.path.abspath(os.path.dirname(arg[1]))
            if self.dirname:
                self.GSASprojectfile = os.path.split(self.GSASprojectfile)[1]
                os.chdir(self.dirname)
                self.LastGPXdir = self.dirname
            try:
                #open the file if possible
                if sys.platform == "darwin": # on Mac delay a bit so GUI can open
                    wx.CallAfter(self.StartProject)
                else:
                    self.StartProject()
                return
            except Exception:
                print ('Error opening or reading file'+arg[1])
                import traceback
                print (traceback.format_exc())
                
        if GSASIIpath.GetConfigValue('Starting_directory'):
            try:
                pth = GSASIIpath.GetConfigValue('Starting_directory')
                pth = os.path.expanduser(pth) 
                os.chdir(pth)
                self.LastGPXdir = pth
            except:
                print('Ignoring Config Starting_directory value: '+
                      GSASIIpath.GetConfigValue('Starting_directory'))
                

    def GetTreeItemsList(self,item):
        return self.GPXtree._getTreeItemsList(item)

    # def OnSize(self,event):
    #     'Called to make GPXtree fill mainPanel'
    #     print 'OnSize'
    #     event.Skip()
        # w,h = self.GetClientSizeTuple()
        # self.dataWindow.SetupScrolling()
        # self.mainPanel.SetSize(wx.Size(w,h))
        # self.GPXtree.SetSize(wx.Size(w,h))
        # self.dataWindow.SetSize(self.dataPanel.GetClientSize())
        
    def SetDataSize(self):
        '''this routine is a placeholder until all G2frame.SetDataSize calls are replaced
        by G2frame.dataWindow.SetDataSize
        '''
        # TOTO: diagnostic patch
        print ('G2frame.SetDataSize called rather than dataWindow.SetDataSize')
        G2obj.HowDidIgetHere(True)
        self.dataWindow.SetDataSize()
                                
    def OnDataTreeSelChanged(self, event):
        '''Called when a data tree item is selected. May be called on item deletion as well.
        '''
        if self.TreeItemDelete:
            self.TreeItemDelete = False
        else:
            if self.ExpandingAll:
                if GSASIIpath.GetConfigValue('debug'): print('Skipping Tree selection due to ExpandAll')
                return
            pltNum = self.G2plotNB.nb.GetSelection()
            if pltNum >= 0:                         #to avoid the startup with no plot!
                self.G2plotNB.nb.GetPage(pltNum)
            item = event.GetItem()
            wx.CallAfter(SelectDataTreeItem,self,item,self.oldFocus)
            #if self.oldFocus: # now done via last parameter on SelectDataTreeItem
            #    wx.CallAfter(self.oldFocus.SetFocus)
        
    def OnGPXtreeItemCollapsed(self, event):
        'Called when a tree item is collapsed - all children will be collapsed'
        self.GPXtree.CollapseAllChildren(event.GetItem())

    def OnGPXtreeItemExpanded(self, event):
        'Called when a tree item is expanded'
        event.Skip()
        
    def OnGPXtreeItemDelete(self, event):
        'Called when a tree item is deleted, inhibit the next tree item selection action'
        self.TreeItemDelete = True

    def OnGPXtreeItemActivated(self, event):
        'Called when a tree item is activated'
        event.Skip()
        
    def OnGPXtreeBeginRDrag(self,event):
        event.Allow()
        self.BeginDragId = event.GetItem()
        self.ParentId = self.GPXtree.GetItemParent(self.BeginDragId)
        DragText = self.GPXtree.GetItemText(self.BeginDragId)
        self.DragData = [[DragText,self.GPXtree.GetItemPyData(self.BeginDragId)],]
        item, cookie = self.GPXtree.GetFirstChild(self.BeginDragId)
        while item:     #G2 data tree has no sub children under a child of a tree item
            name = self.GPXtree.GetItemText(item)
            self.DragData.append([name,self.GPXtree.GetItemPyData(item)])
            item, cookie = self.GPXtree.GetNextChild(self.BeginDragId, cookie)                            
        
    def OnGPXtreeEndDrag(self,event):
        event.Allow()
        self.EndDragId = event.GetItem()
        try:
            NewParent = self.GPXtree.GetItemParent(self.EndDragId)
        except:
            self.EndDragId = self.GPXtree.GetLastChild(self.root)
            NewParent = self.root
        if self.ParentId != NewParent:
            self.ErrorDialog('Drag not allowed','Wrong parent for item dragged')
        else:
            Name,Item = self.DragData[0]
            NewId = self.GPXtree.InsertItem(self.ParentId,self.EndDragId,Name,data=None)
            self.GPXtree.SetItemPyData(NewId,Item)
            for name,item in self.DragData[1:]:     #loop over children
                Id = self.GPXtree.AppendItem(parent=NewId,text=name)
                self.GPXtree.SetItemPyData(Id,item)
            self.GPXtree.Delete(self.BeginDragId)
            SelectDataTreeItem(self,NewId)
        
    def OnGPXtreeKeyDown(self,event): #doesn't exactly work right with Shift key down
        'Allows stepping through the tree with the up/down arrow keys'
        self.oldFocus = wx.Window.FindFocus()
        keyevt = event.GetKeyEvent()
        key = event.GetKeyCode()
        item = self.GPXtree.GetSelection()
        if type(item) is int: return # is this the toplevel in tree?
        name = self.GPXtree.GetItemText(item)
        parent = self.GPXtree.GetItemParent(item)
        if key == wx.WXK_UP:
            if keyevt.GetModifiers() == wx.MOD_SHIFT and parent != self.root:
                if type(parent) is int: return # is this the toplevel in tree?
                prev = self.GPXtree.GetPrevSibling(parent)
                NewId = GetGPXtreeItemId(self,prev,name)
                if NewId:
                    self.GPXtree.Collapse(parent)
                    self.GPXtree.Expand(prev)
                    self.oldFocus = wx.Window.FindFocus()
                    wx.CallAfter(self.GPXtree.SelectItem,NewId)
                else:
                    wx.CallAfter(self.GPXtree.SelectItem,item)
            elif sys.platform == "win32":    
                self.GPXtree.GetPrevSibling(item)
                self.GPXtree.SelectItem(item)
            else:
                item = self.GPXtree.GetPrevSibling(item)
                if item.IsOk(): self.GPXtree.SelectItem(item)
        elif key == wx.WXK_DOWN:
            if keyevt.GetModifiers() == wx.MOD_SHIFT and parent != self.root:
                prev = self.GPXtree.GetNextSibling(parent)
                NewId = GetGPXtreeItemId(self,prev,name)
                if NewId:
                    self.GPXtree.Collapse(parent)
                    self.GPXtree.Expand(prev)
                    self.oldFocus = wx.Window.FindFocus()
                    wx.CallAfter(self.GPXtree.SelectItem,NewId)
                else:
                    wx.CallAfter(self.GPXtree.SelectItem,item)
            elif sys.platform == "win32":    
                self.GPXtree.GetNextSibling(item)
                self.GPXtree.SelectItem(item)
            else:    
                item = self.GPXtree.GetNextSibling(item)
                if item.IsOk(): self.GPXtree.SelectItem(item)
                    
    def OnColMetaTest(self,event):
        'Test the .par/.*lbls pair for contents'
        G2imG.testColumnMetadata(self)
                
    def OnPowderFPA(self,event):
        'Perform FPA simulation/peak fitting'
        G2fpa.GetFPAInput(self)
        
    def OnReadPowderPeaks(self,event):
        'Bound to menu Data/Read Powder Peaks'
        self.CheckNotebook()
        pth = G2G.GetImportPath(self)
        if not pth: pth = '.'
        dlg = wx.FileDialog(self, 'Choose file with peak list', pth, '', 
            'peak files (*.txt)|*.txt|All files (*.*)|*.*',wx.FD_MULTIPLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                for file_ajk in dlg.GetPaths():
                    self.HKL = []
                    self.powderfile = file_ajk
                    comments,peaks,limits,wave = G2IO.GetPowderPeaks(self.powderfile)
                    Id = self.GPXtree.AppendItem(parent=self.root,text='PKS '+os.path.basename(self.powderfile))
                    data = ['PKS',wave,0.0]
                    names = ['Type','Lam','Zero'] 
                    codes = [0,0,0]
                    inst = [G2fil.makeInstDict(names,data,codes),{}]
                    self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Instrument Parameters'),inst)
                    self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Comments'),comments)
                    self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Limits'),[tuple(limits),limits])
                    self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Index Peak List'),[peaks,[]])
                    self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Unit Cells List'),[])             
                    self.GPXtree.Expand(Id)
                    self.GPXtree.SelectItem(Id)
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
        finally:
            dlg.Destroy()
                        
    def CheckNotebook(self):
        '''Make sure the data tree has the minimally expected controls.
        '''
        new = False
        if not GetGPXtreeItemId(self,self.root,'Notebook'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Notebook')
            self.GPXtree.SetItemPyData(sub,[''])
        if not GetGPXtreeItemId(self,self.root,'Controls'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Controls')
            self.GPXtree.SetItemPyData(sub,copy.copy(G2obj.DefaultControls))
        if not GetGPXtreeItemId(self,self.root,'Covariance'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Covariance')
            self.GPXtree.SetItemPyData(sub,{})
        if not GetGPXtreeItemId(self,self.root,'Constraints'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Constraints')
            self.GPXtree.SetItemPyData(sub,{'Hist':[],'HAP':[],'Phase':[]})
        if not GetGPXtreeItemId(self,self.root,'Restraints'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Restraints')
            self.GPXtree.SetItemPyData(sub,{})
        if not GetGPXtreeItemId(self,self.root,'Rigid bodies'):
            new = True
            sub = self.GPXtree.AppendItem(parent=self.root,text='Rigid bodies')
            self.GPXtree.SetItemPyData(sub,{'Vector':{'AtInfo':{}},
                'Residue':{'AtInfo':{}},'RBIds':{'Vector':[],'Residue':[]}})
        if new:
            self.GPXtree.Expand(self.GPXtree.root)
                
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
            dataGridSizer = wx.FlexGridSizer(cols=ncols,hgap=2,vgap=2)
            for Id,item in enumerate(self.data):
                ckbox = wx.CheckBox(panel,Id,item[1])
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
            Id = event.GetId()
            self.data[Id][0] = self.FindWindowById(Id).GetValue()        
            
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
        '''Allows user to supply scale factor(s) when summing data 
        '''
        def __init__(self,parent,title,text,dataType,data,dataList,Limits=None):
            wx.Dialog.__init__(self,parent,-1,title,size=(400,250),
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
            self.plotFrame = wx.Frame(self,-1,'Sum Plots',size=wx.Size(700,600), \
                style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
            self.G2plotNB = G2plt.G2PlotNoteBook(self.plotFrame,G2frame=self)
            self.text = text
            self.data = data
            self.average = False
            self.selectData = copy.copy(data[:-1])
            self.selectVals = len(data)*[0.0,]
            self.dataList = dataList
            self.Limits = Limits
            self.filterlist = range(len(self.dataList)) # list of the choice numbers that have been filtered (list of int indices)
            self.dataType = dataType
            self.filterVal = ''
            self.panel = None
            self.Draw()
            
        def Draw(self):            
            if self.panel:
                self.panel.DestroyChildren()  #safe: wx.Panel
                self.panel.Destroy()
            size = (480,350)
            self.panel = wxscroll.ScrolledPanel(self, wx.ID_ANY,size=size,
                style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
            mainSizer = wx.BoxSizer(wx.VERTICAL)
            mainSizer.Add(wx.StaticText(self.panel,label=self.text),0)
            mainSizer.Add((10,10))
            self.dataGridSizer = wx.FlexGridSizer(cols=2,hgap=2,vgap=2)
            self.dataGridSizer.Add((-1,-1))
            topSizer = wx.BoxSizer(wx.HORIZONTAL)
            topSizer.Add((-1,-1),1,wx.EXPAND,1)
            topSizer.Add(wx.StaticText(self.panel,label='Filter:  '),0,WACV)
            self.timer = wx.Timer()
            self.timer.Bind(wx.EVT_TIMER,self.OnFilter)
            self.filterBox = wx.TextCtrl(self.panel, wx.ID_ANY, self.filterVal,
                                         size=(80,-1),style=wx.TE_PROCESS_ENTER)
            self.filterBox.Bind(wx.EVT_TEXT,self.onChar)
            self.filterBox.Bind(wx.EVT_TEXT_ENTER,self.OnFilter)
            topSizer.Add(self.filterBox,0,WACV)
            self.dataGridSizer.Add(topSizer,1,wx.RIGHT|wx.BOTTOM|wx.EXPAND,1)
            self.dataGridSizer.Add((-1,10))
            self.dataGridSizer.Add((-1,10))
            for Id,item in enumerate(self.selectData):
                name = wx.TextCtrl(self.panel,-1,item,size=wx.Size(300,20))
                name.SetEditable(False)
                scale = G2G.ValidatedTxtCtrl(self.panel,self.selectVals,Id,nDig=(10,3),typeHint=float)
                self.dataGridSizer.Add(scale,0,wx.LEFT,10)
                self.dataGridSizer.Add(name,0,wx.RIGHT,10)
            if self.dataType:
                ScaleAll = wx.Button(self.panel,wx.ID_ANY,'Set all above')
                ScaleAll.Bind(wx.EVT_BUTTON, self.OnAllScale)
                if self.dataType == 'PWDR':
                    self.Avg = wx.CheckBox(self.panel,label=' Make average?')
                    self.Avg.Bind(wx.EVT_CHECKBOX,self.OnAve)
                self.dataGridSizer.Add(ScaleAll,0,wx.LEFT,10)
                if self.dataType == 'PWDR':
                    self.dataGridSizer.Add(self.Avg,0,wx.RIGHT,10)
                self.dataGridSizer.Add(wx.StaticText(self.panel,-1,' Result type: '+self.dataType),1,
                    wx.LEFT|wx.ALIGN_CENTER_VERTICAL,1)
            mainSizer.Add(self.dataGridSizer,0,wx.EXPAND)
            self.name = G2G.ValidatedTxtCtrl(self.panel,self.data,-1,size=wx.Size(300,20))
            mainSizer.Add(self.name,0,wx.RIGHT|wx.TOP,10)
            self.OkBtn = wx.Button(self.panel,-1,"Ok")
            self.OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
            cancelBtn = wx.Button(self.panel,-1,"Cancel")
            cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
            btnSizer = wx.FlexGridSizer(0,3,10,20)
            if self.dataType =='PWDR':
                TestBtn = wx.Button(self.panel,-1,"Test")
                TestBtn.Bind(wx.EVT_BUTTON, self.OnTest)
                btnSizer.Add(TestBtn)
            btnSizer.Add(self.OkBtn)
            btnSizer.Add(cancelBtn)
            btnSizer.Add((5,5))
            
            self.panel.SetSizer(mainSizer)
            self.panel.SetAutoLayout(1)
            self.panel.SetupScrolling()
            mainSizer.Add((10,10),1)
            mainSizer.Add(btnSizer,0,wx.CENTER)
            self.panel.SetSizer(mainSizer)
            self.panel.Fit()
            self.Fit()
            
        def OnAve(self,event):
            self.average = self.Avg.GetValue()

        def OnFilter(self,event):
            '''Read text from filter control and select entries that match.
            '''
            if self.timer.IsRunning():
                self.timer.Stop()
            self.filterVal = txt = self.filterBox.GetValue()
            if txt:
                txt = txt.lower()
                ChoiceList = []
                ChoiceVals = []
                for i,item in enumerate(self.selectData):
                    if item.lower().find(txt) != -1:
                        ChoiceList.append(item)
                        ChoiceVals.append(self.selectVals[i])
                self.selectData = ChoiceList
                self.selectVals = ChoiceVals
            else:
#                self.selectData = copy.copy(self.data[:-1])
                self.selectVals = len(self.data)*[0.0,]
            wx.CallAfter(self.Draw)
            
        def GetData(self):
            if self.dataType == 'PWDR':
                return self.selectData+[self.data[-1],],self.result
            else:
                return self.selectData+[self.data[-1],],self.selectVals
                        
        def onChar(self,event):
            'Respond to keyboard events in the Filter box'
            self.filterVal = self.filterBox.GetValue()
            if self.timer.IsRunning():
                self.timer.Stop()
            self.timer.Start(1000,oneShot=True)
            if event: event.Skip()

        def OnAllScale(self,event):
            dlg = G2G.SingleFloatDialog(self,'New scale',
                                        'Enter new value for all scale factors',1.)
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                val = dlg.GetValue()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
            for Id,item in enumerate(self.selectData):
                self.selectVals[Id] = val
            wx.CallAfter(self.Draw)
            
        def OnTest(self,event):
            lenX = 0
            Xminmax = [0,0]
            XY = []
            Xsum = []
            Ysum = []
            Vsum = []
            for i,item in enumerate(self.selectData):
                name = item
                scale = self.selectVals[i]
                Id = self.data.index(name)
                data = self.dataList[Id]
                if scale:
                    x,y,w,yc,yb,yd = data   #numpy arrays!
                    if self.Limits is not None:
                        xMin = np.searchsorted(x,self.Limits[1][0])
                        xMax = np.searchsorted(x,self.Limits[1][1])
                        x = x[xMin:xMax+1]
                        y = y[xMin:xMax+1]
                        lenX = xMax-xMin+1
                    XY.append([x,scale*y])
                    v = 1./w[xMin:xMax+1]
                    if lenX:
                        if lenX != len(x):
                            self.GetParent().ErrorDialog('Data length error','Data to be summed must have same number of points'+
                                '\nExpected:'+str(lenX)+
                                '\nFound:   '+str(len(x))+'\nfor '+name)
                            return
#                            self.OnCancel(event)
                    else:
                        lenX = len(x)
                    if Xminmax[1]:
                        if Xminmax != [x[0],x[-1]]:
                            self.GetParent().ErrorDialog('Data range error','Data to be summed must span same range'+
                                '\nExpected:'+str(Xminmax[0])+' '+str(Xminmax[1])+
                                '\nFound:   '+str(x[0])+' '+str(x[-1])+'\nfor '+name)
                            return
#                            self.OnCancel(event)
                    else:
                        Xminmax = [x[0],x[-1]]
                        Xsum = x
                    if self.dataType == 'PWDR' and self.average:
                        Ysum.append(scale*y)
                        Vsum.append(abs(scale)*v)
                    else:
                        try:
                            Ysum += scale*y
                            Vsum += abs(scale)*v
                        except ValueError:
                            Ysum = scale*y
                            Vsum = abs(scale)*v
            if self.dataType =='PWDR' and self.average:
                maYsum = ma.masked_equal(Ysum,0)
                Ysum = ma.mean(maYsum,axis=0)
                Wsum = 1./np.array(Ysum)
            else:
                Wsum = 1./Vsum
            YCsum = np.zeros(lenX)
            YBsum = np.zeros(lenX)
            YDsum = np.zeros(lenX)
            XY.append([Xsum,Ysum])
            self.result = [Xsum,Ysum,Wsum,YCsum,YBsum,YDsum]
            # N.B. PlotXY expects the first arg to point to G2frame. In this case, we
            # create a duplicate (temporary) Plot notebook window that is a child of the
            # modal SumDialog dialog (self). This nicely gets deleted when the dialog is destroyed,
            # but the plot window is not fully functional, at least on the Mac.
            if len(XY[0][0]):
                G2plt.PlotXY(self,XY,lines=True,Title='Sum:'+self.data[-1],labelY='Intensity',)
                self.plotFrame.Show()
                return True
                        
        def OnOk(self,event):
            if self.dataType == 'PWDR':
                if not self.OnTest(event): return
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_OK)              
            
        def OnCancel(self,event):
            parent = self.GetParent()
            parent.Raise()
            self.EndModal(wx.ID_CANCEL)              
            
    def OnPwdrSum(self,event):
        'Sum or Average together powder data(?)'
        TextList = []
        DataList = []
        Limits = []
        Names = []
        Inst = None
        Comments = ['Sum/Average equals: \n']
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append(name)
                    DataList.append(self.GPXtree.GetItemPyData(item)[1])    # (x,y,w,yc,yb,yd)
                    if not Inst:
                        Inst = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,item, 'Instrument Parameters'))
                    if not Limits:
                        Limits = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,item, 'Limits'))
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum/average','There must be more than one "PWDR" pattern')
                return
            TextList.append('default_ave_name')                
            dlg = self.SumDialog(self,'Sum/Average data','''
    Enter scale for each pattern to be summed/averaged
    Limits for first pattern used sets range for the sum
    All patterns used must extend over this range
                ''','PWDR',
                TextList,DataList,Limits)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result,sumData = dlg.GetData()
                    Xsum,Ysum,Wsum,YCsum,YBsum,YDsum = sumData
                    Xminmax = [Xsum[0],Xsum[-1]]
                    outname = 'PWDR '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = GetGPXtreeItemId(self,self.root,name)
                                self.GPXtree.Delete(Id)
                        finally:
                            dlg2.Destroy()
                    Id = self.GPXtree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        Sample = G2obj.SetDefaultSample()
                        Ymin = np.min(Ysum)
                        Ymax = np.max(Ysum)
                        valuesdict = {
                            'wtFactor':1.0,
                            'Dummy':False,
                            'ranId':ran.randint(0,sys.maxsize),
                            'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-.1*Ymax,'refDelt':0.1*Ymax,
                            'Yminmax':[Ymin,Ymax]
                            }
                        self.GPXtree.SetItemPyData(Id,[valuesdict,[np.array(Xsum),np.array(Ysum),np.array(Wsum),
                            np.array(YCsum),np.array(YBsum),np.array(YDsum)]])
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Comments'),Comments)                    
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Limits'),[tuple(Xminmax),Xminmax])
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Background'),[['chebyschev-1',True,3,1.0,0.0,0.0],
                            {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}])
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Instrument Parameters'),Inst)
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Sample Parameters'),Sample)
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Peak List'),{'peaks':[],'sigDict':{}})
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Index Peak List'),[[],[]])
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Unit Cells List'),[])             
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Reflection Lists'),{})             
                        self.GPXtree.SelectItem(Id)
                        self.GPXtree.Expand(Id)
            finally:
                dlg.Destroy()

    def OnImageSum(self,event):
        'Sum together image data'
        TextList = []
        DataList = []
        IdList = []
        Names = []
        Comments = ['Sum equals: \n']
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                Names.append(name)
                if 'IMG' in name:
                    TextList.append(name)
                    DataList.append(self.GPXtree.GetImageLoc(item))        #Size,Image,Tag
                    IdList.append(item)
                    Data = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,item,'Image Controls'))
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
            if len(TextList) < 2:
                self.ErrorDialog('Not enough data to sum','There must be more than one "IMG" pattern')
                return
            TextList.append('default_sum_name')                
            dlg = self.SumDialog(self,'Sum data',' Enter scale for each image to be summed','IMG',
                TextList,DataList)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    imSize = 0
                    result,scales = dlg.GetData()
                    First = True
                    Found = False
                    for name,scale in zip(result,scales):
                        if scale:
                            Found = True                                
                            Comments.append("%10.3f %s" % (scale,' * '+name))
                            i = TextList.index(name)
                            Npix,imagefile,imagetag = DataList[i]
                            imagefile = G2IO.GetCheckImageFile(self,IdList[i])[1]
                            image = G2IO.GetImageData(self,imagefile,imageOnly=True,ImageTag=imagetag)
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
                        
                        
                    newImage = np.array(newImage,dtype=np.int32)                        
                    outname = 'IMG '+result[-1]
                    Id = 0
                    if outname in Names:
                        dlg2 = wx.MessageDialog(self,'Overwrite data?','Duplicate data name',wx.OK|wx.CANCEL)
                        try:
                            if dlg2.ShowModal() == wx.ID_OK:
                                Id = GetGPXtreeItemId(self,self.root,name)
                        finally:
                            dlg2.Destroy()
                    else:
                        Id = self.GPXtree.AppendItem(parent=self.root,text=outname)
                    if Id:
                        pth = os.path.split(os.path.abspath(imagefile))[0]
#                        pth = G2G.GetExportPath(self)
                        dlg = wx.FileDialog(self, 'Choose sum image filename', pth,outname.split('IMG ')[1], 
                            'G2img files (*.G2img)|*.G2img', 
                            wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                        if dlg.ShowModal() == wx.ID_OK:
                            newimagefile = dlg.GetPath()
                            newimagefile = G2IO.FileDlgFixExt(dlg,newimagefile)
                            G2IO.PutG2Image(newimagefile,Comments,Data,Npix,newImage)
                            Imax = np.amax(newImage)
                            Imin = np.amin(newImage)
                            newImage = []
                            self.GPXtree.SetItemPyData(Id,[imSize,newimagefile])
                            self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Comments'),Comments)
                        del(newImage)
                        if self.imageDefault:
                            Data = copy.copy(self.imageDefault)
                        Data['formatName'] = 'GSAS-II image'
                        Data['showLines'] = True
                        Data['ring'] = []
                        Data['rings'] = []
                        Data['cutoff'] = 10
                        Data['pixLimit'] = 20
                        Data['ellipses'] = []
                        Data['calibrant'] = ''
                        Data['range'] = [(Imin,Imax),[Imin,Imax]]
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Image Controls'),Data)                                            
                        Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],
                            'Frames':[],'Thresholds':[(Imin,Imax),[Imin,Imax]],
                                     'SpotMask':{'esdMul':2.,'spotMask':None}}
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Masks'),Masks)
                        self.GPXtree.SetItemPyData(self.GPXtree.AppendItem(Id,text='Stress/Strain'),
                            {'Type':'True','d-zero':[],'Sample phi':0.0,'Sample z':0.0,'Sample load':0.0})
                        self.GPXtree.SelectItem(Id)
                        self.GPXtree.Expand(Id)
                        self.PickId = GetGPXtreeItemId(self,self.root,outname)
                        self.Image = self.PickId
            finally:
                dlg.Destroy()
                      
    def OnAddPhase(self,event):
        'Add a new, empty phase to the tree. Called by Data/Add Phase menu'
        self.CheckNotebook()
        if not GetGPXtreeItemId(self,self.root,'Phases'):
            sub = self.GPXtree.AppendItem(parent=self.root,text='Phases')
        else:
            sub = GetGPXtreeItemId(self,self.root,'Phases')
        PhaseName = ''
        dlg = wx.TextEntryDialog(None,'Enter a name for this phase','Phase Name Entry','New phase',
            style=wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            PhaseName = dlg.GetValue()
        dlg.Destroy()
        if not GetGPXtreeItemId(self,self.root,'Restraints'):
            subr = self.GPXtree.AppendItem(parent=self.root,text='Restraints')
            self.GPXtree.SetItemPyData(subr,{PhaseName:{}})
        else:
            subr = GetGPXtreeItemId(self,self.root,'Restraints')
            self.GPXtree.GetItemPyData(subr).update({PhaseName:{}})
        self.GPXtree.AppendItem(parent=subr,text=PhaseName)
        newphase = self.GPXtree.AppendItem(parent=sub,text=PhaseName)
        E,SGData = G2spc.SpcGroup('P 1')
        self.GPXtree.SetItemPyData(newphase,G2obj.SetNewPhase(Name=PhaseName,SGData=SGData))
        self.GPXtree.Expand(sub)
        SelectDataTreeItem(self,newphase) #bring up new phase General tab
        
    def OnDeletePhase(self,event):
        '''Delete one or more phases from the tree. Called by Data/Delete Phase menu.
        Also delete this phase from Reflection Lists for each PWDR histogram; 
        removes the phase from restraints and deletes any constraints 
        with variables from the phase.
        If any deleted phase is marked as Used in a histogram, a more rigorous
        "deep clean" is done and histogram refinement results are cleared, as well as 
        the covariance information and all plots are deleted
        '''
        selItem = self.GPXtree.GetSelection()
        if self.dataWindow:
            self.dataWindow.ClearData() 
        TextList = []
        DelList = []
        DelItemList = []
        consDeleted = 0
        usedPhase = False
        if GetGPXtreeItemId(self,self.root,'Phases'):
            sub = GetGPXtreeItemId(self,self.root,'Phases')
        else:
            return
        if GetGPXtreeItemId(self,self.root,'Restraints'):
            subr = GetGPXtreeItemId(self,self.root,'Restraints')
        else:
            subr = 0
        if GetGPXtreeItemId(self,self.root,'Constraints'):
            id = GetGPXtreeItemId(self,self.root,'Constraints')
            constr = self.GPXtree.GetItemPyData(id)
        else:
            constr = {}
            
        item, cookie = self.GPXtree.GetFirstChild(sub)
        while item:
            TextList.append(self.GPXtree.GetItemText(item))
            item, cookie = self.GPXtree.GetNextChild(sub, cookie)                
        dlg = wx.MultiChoiceDialog(self, 'Which phase to delete?', 'Delete phase', TextList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: DelList.append([i,TextList[i]])
                item, cookie = self.GPXtree.GetFirstChild(sub)
                i = 0
                while item:
                    if [i,self.GPXtree.GetItemText(item)] in DelList: DelItemList.append(item)
                    item, cookie = self.GPXtree.GetNextChild(sub, cookie)
                    i += 1
                for item in DelItemList:
                    phase = self.GPXtree.GetItemPyData(item)
                    for h in phase['Histograms']:
                        if 'Use' not in phase['Histograms'][h]: continue
                        if phase['Histograms'][h]['Use']:
                            usedPhase = True
                            break
                    if 'pId' in phase:
                        p = phase['pId']
                    else:
                        p = '?'
                    self.GPXtree.Delete(item)
                    if item == selItem: selItem = self.root
                    # look for constraints to remove
                    for key in constr:
                        delThis = []
                        if key.startswith('_'): continue
                        for i,cons in enumerate(constr[key]):
                            for var in cons[0:-3]:
                                if str(var[1]).startswith(str(p)):
                                    delThis.append(i)
                                    break
                        for i in reversed(delThis):
                            consDeleted += 1
                            del constr[key][i]
                # delete refinement results from histograms
                item, cookie = self.GPXtree.GetFirstChild(self.root)
                while item:
                    name = self.GPXtree.GetItemText(item)
                    if 'PWDR' in name:
                        data = self.GPXtree.GetItemPyData(item)
                        if usedPhase: # remove r-factors
                            dellist = [value for value in data[0] if ':' in value]
                            for v in dellist+['Durbin-Watson', 'R', 'wR', 'Rb',
                                                  'wRb', 'wRmin','Nobs']:
                                if v in data[0]: del data[0][v]
                            # could wipe out computed & difference patterns, but does not work
                            #data[1][3] = np.zeros_like(data[1][3])
                            #data[1][5] = np.zeros_like(data[1][5])
                        # always get rid of reflection lists 
                        Id = GetGPXtreeItemId(self,item, 'Reflection Lists')
                        refList = self.GPXtree.GetItemPyData(Id)
                        if len(refList):
                            for i,item in DelList:
                                if item in refList:
                                    del(refList[item])
                    elif 'HKLF' in name and usedPhase: # probably not needed if phase is not used
                        data = self.GPXtree.GetItemPyData(item)
                        data[0] = {}

                    item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        finally:
            dlg.Destroy()
        if usedPhase: # clear info from last refinement for "deep clean" if a used phase is deleted
            id = GetGPXtreeItemId(self,self.root,'Covariance')
            if DelItemList and id:
                self.GPXtree.SetItemPyData(id,{})
            id = GetGPXtreeItemId(self,self.root,'Sequential results')
            if DelItemList and id:
                self.GPXtree.Delete(id)
                if id == selItem: selItem = self.root
            # delete all plots
            for lbl in self.G2plotNB.plotList:
                self.G2plotNB.Delete(lbl)
        if subr and DelList:     #remove restraints for deleted phase
            DelList = [itm[1] for itm in DelList]
            item, cookie = self.GPXtree.GetFirstChild(subr)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name in DelList:
                    self.GPXtree.Delete(item)
                    if item == selItem: selItem = self.root
                item, cookie = self.GPXtree.GetNextChild(subr, cookie)
        # force redisplay of current tree item if it was not deleted
        self.PickId = 0
        self.PatternId = 0
        self.PickIdText = None
        SelectDataTreeItem(self,selItem)
        wx.CallAfter(self.GPXtree.SelectItem,selItem)
        if consDeleted:
            print('\n',consDeleted,'constraints were deleted')
        
    def OnRenameData(self,event):
        '''Renames an existing histogram. Called by Data/Rename Phase menu.
        Must be used before a histogram is used in a phase.
        '''
        name = self.GPXtree.GetItemText(self.PickId)      
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        if name in Histograms:
            G2G.G2MessageBox(self,
                'Histogram is used. You must remove it from all phases before it can be renamed',
                'Rename not allowed')
            return
        if 'PWDR' in name or 'HKLF' in name or 'IMG' in name:
            if 'Bank' in name:
                names = name.split('Bank')
                names[1] = ' Bank'+names[1]
            elif 'Azm' in name:
                names = name.split('Azm')
                names[1] = ' Azm'+names[1]
            else:
                names = [name,'']
            dataType = names[0][:names[0].index(' ')+1]                 #includes the ' '
            dlg = G2G.SingleStringDialog(self,'Change tree name',
                    'Data name: '+name,names[0][names[0].index(' ')+1:])
            #if dlg.ShowModal() == wx.ID_OK:
            if dlg.Show():
                name = dataType+dlg.GetValue().strip()+names[1]
                self.GPXtree.SetItemText(self.PickId,name)
                if 'PWDR' in name:
                    self.GPXtree.GetItemPyData(self.PickId)[2] = name
            dlg.Destroy()
        
    def GetFileList(self,fileType,skip=None):        #potentially useful?
        'Appears unused. Note routine of same name in GSASIIpwdGUI'
        fileList = []
        Source = ''
        Id, cookie = self.GPXtree.GetFirstChild(self.root)
        while Id:
            name = self.GPXtree.GetItemText(Id)
            if fileType in name:
                if Id == skip:
                    Source = name
                else:
                    fileList.append([False,name,Id])
            Id, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        if skip:
            return fileList,Source
        else:
            return fileList
            
    def OnDataDelete(self, event):
        '''Delete one or more histograms from data tree. Called by the
        Data/DeleteData menu
        '''
        TextList = []
        DelList = []
        DelItemList = []
        nItems = {'PWDR':0,'SASD':0,'REFD':0,'IMG':0,'HKLF':0,'PDF':0}
        PDFnames = []
        selItem = self.GPXtree.GetSelection()
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        if not self.GPXtree.GetCount():
            G2G.G2MessageBox(self,'No tree items to be deleted',
                                 'Nothing to delete')
            return            
        item, cookie = self.GPXtree.GetFirstChild(self.root)
        used = False
        while item:
            name = self.GPXtree.GetItemText(item)
            if name not in ['Notebook','Controls','Covariance','Constraints',
                'Restraints','Phases','Rigid bodies'] and 'Sequential' not in name:
                if 'PWDR' in name[:4]:
                    nItems['PWDR'] += 1
                    if name in Histograms:
                        used = True
                        item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
                        continue
                if 'SASD' in name[:4]: nItems['SASD'] += 1
                if 'REFD' in name[:4]: nItems['REFD'] += 1
                if 'IMG' in name[:3]:  nItems['IMG'] += 1
                if 'HKLF' in name[:4]:
                    nItems['HKLF'] += 1
                    if name in Histograms:
                        used = True
                        item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
                        continue
                if 'PDF' in name[:3]:
                    PDFnames.append(name)
                    nItems['PDF'] += 1
                TextList.append(name)
            item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        for pdfName in PDFnames:
            try:
                TextList.remove('PWDR'+pdfName[4:])
            except ValueError:
                print (u'PWDR'+pdfName[4:]+u' for '+pdfName+u' not found')
        if len(TextList) == 0 and used:
            G2G.G2MessageBox(self,'All histograms are used. You must remove them from phases before they can be deleted',
                                 'Nothing to delete')
            return
        elif len(TextList) == 0:
            G2G.G2MessageBox(self,'None of the tree items are allowed to be deleted',
                                 'Nothing to delete')
            return
        
        dlg = G2G.G2MultiChoiceDialog(self, 'Which data to delete?', 'Delete data', TextList, wx.CHOICEDLG_STYLE)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                result = dlg.GetSelections()
                for i in result: DelList.append(TextList[i])
                item, cookie = self.GPXtree.GetFirstChild(self.root)
                while item:
                    itemName = self.GPXtree.GetItemText(item)
                    if itemName in DelList:
                        if 'PWDR' in itemName[:4]: nItems['PWDR'] -= 1
                        elif 'SASD' in itemName[:4]: nItems['SASD'] -= 1
                        elif 'REFD' in itemName[:4]: nItems['REFD'] -= 1
                        elif 'IMG' in itemName[:3]: nItems['IMG'] -= 1
                        elif 'HKLF' in itemName[:4]: nItems['HKLF'] -= 1
                        elif 'PDF' in itemName[:3]: nItems['PDF'] -= 1
                        DelItemList.append(item)
                    item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
                for item in DelItemList:
                    self.GPXtree.Delete(item)
                    if item == selItem: selItem = self.root
                if DelList:
                    self.PickId = 0
                    self.PickIdText = None
                    self.PatternId = 0
                    if nItems['PWDR']:
                        wx.CallAfter(G2plt.PlotPatterns,self,True)
                    else:
                        self.G2plotNB.Delete('Powder Patterns')
                        self.lastPlotType = None
                    if not nItems['IMG']:
                        self.G2plotNB.Delete('2D Powder Image')
                    if not nItems['HKLF']:
                        self.G2plotNB.Delete('Structure Factors')
                        if '3D Structure Factors' in self.G2plotNB.plotList:
                            self.G2plotNB.Delete('3D Structure Factors')
        finally:
            dlg.Destroy()
        if DelList:
            SelectDataTreeItem(self,selItem)
            wx.CallAfter(self.GPXtree.SelectItem,selItem)
                
    def OnPlotDelete(self,event):
        '''Delete one or more plots from plot window. Called by the
        Data/DeletePlots menu
        '''
        plotNames = self.G2plotNB.plotList
        if len(plotNames):
            dlg = G2G.G2MultiChoiceDialog(self, 'Which plots to delete?', 'Delete plots', plotNames, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    result = dlg.GetSelections()
                    result.sort(reverse=True)
                    for i in result: 
                        self.G2plotNB.Delete(plotNames[i])
            finally:
                dlg.Destroy()

    def OnFileReopen(self, event):
        files = GSASIIpath.GetConfigValue('previous_GPX_files')
        if not files:
            print('no previous projects found')
            return
        sellist = []
        for f in files:
            dirname,filroot = os.path.split(f)
            if os.path.exists(f) and '.gpx' in f:
                sellist.append("{} from {}".format(filroot,dirname))
#            else:
#                sellist.append("not found: {}".format(f))
        
        dlg = G2G.G2SingleChoiceDialog(self,
                                           'Select previous project to open',
                                           'Select project',sellist)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelection()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        filroot,dirname = sellist[sel].split(' from ')
        f = os.path.join(dirname,filroot)
        if os.path.exists(f):
            self.OnFileOpen(event, filename=f)
            self.LastGPXdir = dirname
        else:
            print('file not found',f)        
        
    def OnFileOpen(self, event, filename=None):
        '''Gets a GSAS-II .gpx project file in response to the
        File/Open Project menu button
        '''
        def SaveOld():
            '''See if we should save current project and continue 
            to read another.
            returns True if the project load should continue
            '''
            if self.dataWindow:
                self.dataWindow.ClearData()
            dlg = wx.MessageDialog(self,
                    'Do you want to save and replace the current project?\n(Use No to read without saving or Cancel to continue with current project)',
                    'Save & Overwrite?',
                    wx.YES|wx.NO|wx.CANCEL)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result == wx.ID_NO:
                result = True
            elif result == wx.ID_CANCEL:
                return False
            else:
                if not self.OnFileSave(None): return False
            self.GPXtree.DeleteChildren(self.root)
            self.GSASprojectfile = ''
            self.HKL = []
            if self.G2plotNB.plotList:
                self.G2plotNB.clear()
            return True
        
        def GetGPX():
            if self.LastGPXdir:
                pth = self.LastGPXdir
            else:
                pth = '.'
            #if GSASIIpath.GetConfigValue('debug'): print('debug: open from '+pth)
            dlg = wx.FileDialog(self, 'Choose GSAS-II project file', pth, 
                wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN)
            try:
                if dlg.ShowModal() != wx.ID_OK: return
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.LastGPXdir = dlg.GetDirectory()
            finally:
                dlg.Destroy()
            
        self.EnablePlot = False
        if self.GPXtree.GetChildrenCount(self.root,False):
            if not SaveOld(): return

        if not filename:
            GetGPX()
            filename = self.GSASprojectfile
        else:
            try:
                self.GSASprojectfile = os.path.splitext(filename)[0]+u'.gpx'
            except:
                self.GSASprojectfile = os.path.splitext(filename)[0]+'.gpx'
            self.dirname = os.path.split(filename)[0]

        self.init_vars()
        try:
            self.StartProject()         #open the file if possible
        except:
            print ('\nError opening file '+filename)
            import traceback
            print (traceback.format_exc())

    def StartProject(self):
        '''Opens a GSAS-II project file & selects the 1st available data set to 
        display (PWDR, HKLF, REFD or SASD)
        '''

        Id = 0
        phaseId = None
        G2IO.ProjFileOpen(self)
        self.GPXtree.SetItemText(self.root,'Project: '+self.GSASprojectfile)
        self.GPXtree.Expand(self.root)
        self.HKL = []
        item, cookie = self.GPXtree.GetFirstChild(self.root)
        while item:
            name = self.GPXtree.GetItemText(item)
            if name[:4] in ['PWDR','HKLF','IMG ','PDF ','SASD','REFD']:
                if not Id:
                    if name[:4] == 'IMG ':
                        Id = GetGPXtreeItemId(self,item,'Image Controls')
                    else:
                        Id = item
            elif name == "Phases":
                phaseId = item
            elif name == 'Controls':
                data = self.GPXtree.GetItemPyData(item)
                if data:
                    for item in self.Refine: item.Enable(True)
            item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        if phaseId: # show all phases
            self.GPXtree.Expand(phaseId)
        if Id:
            self.EnablePlot = True
            self.GPXtree.Expand(Id)
            SelectDataTreeItem(self,Id)
            self.GPXtree.SelectItem(Id)  # needed on OSX or item is not selected in tree; perhaps not needed elsewhere
        elif phaseId:
            Id = phaseId
            # open 1st phase
            Id, unused = self.GPXtree.GetFirstChild(phaseId)
            SelectDataTreeItem(self,Id)
            self.GPXtree.SelectItem(Id) # as before for OSX
        self.CheckNotebook()
        if self.dirname: os.chdir(self.dirname)           # to get Mac/Linux to change directory!
        pth = os.path.split(os.path.abspath(self.GSASprojectfile))[0]
        if GSASIIpath.GetConfigValue('Save_paths'):
            G2G.SaveGPXdirectory(pth,write=False)
        config = G2G.GetConfigValsDocs()
        GSASIIpath.addPrevGPX(self.GSASprojectfile,config)
        G2G.SaveConfigVars(config)
        self.LastGPXdir = pth

    def OnFileClose(self, event):
        '''Clears the data tree in response to the
        File/New Project menu button. User is given option to save
        the project.
        '''
        dlg = wx.MessageDialog(self,
                    'Do you want to save the current project and start with an empty one?\n(Use No to clear without saving or Cancel to continue with current project)',
                    'Save & Clear?',
                    wx.YES | wx.NO | wx.CANCEL)
        try:
            result = dlg.ShowModal()
            if result == wx.ID_OK:
                self.OnFileSaveMenu(event)
            if result != wx.ID_CANCEL:
                self.GSASprojectfile = ''
                self.GPXtree.SetItemText(self.root,'Project: ')
                self.GPXtree.DeleteChildren(self.root)
                self.dataWindow.ClearData()
                if len(self.HKL): self.HKL = []
                if self.G2plotNB.plotList:
                    self.G2plotNB.clear()
                self.SetTitleByGPX()
                self.EnableRefineCommand()
                self.init_vars()
        finally:
            dlg.Destroy()

    def OnFileSave(self, event):
        '''Save the current project in response to the
        File/Save Project menu button
        '''
        
        if self.GSASprojectfile: 
            self.GPXtree.SetItemText(self.root,'Project: '+self.GSASprojectfile)
            self.CheckNotebook()
            G2IO.ProjFileSave(self)
            return True
        else:
            return self.OnFileSaveas(event)
            
    def OnNewGSASII(self, event):
        '''Gets a GSAS-II .gpx project file in response to the
        File/Open new window menu button. Runs only on Mac.
        '''
        if self.LastGPXdir:
            pth = self.LastGPXdir
        else:
            pth = '.'
        GSASprojectfile = ''
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file', pth, 
                wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK: 
                GSASprojectfile = dlg.GetPath()
                GSASprojectfile = G2IO.FileDlgFixExt(dlg,GSASprojectfile)
                self.LastGPXdir = dlg.GetDirectory()
        finally:
            dlg.Destroy()
        G2script = os.path.join(os.path.split(__file__)[0],'GSASII.py')
        GSASIIpath.MacStartGSASII(G2script,GSASprojectfile)

    def SetTitleByGPX(self):
        '''Set the title for the two window frames
        '''
        projName = os.path.split(self.GSASprojectfile)[1]
        if not projName: projName = "<unnamed project>"
        if self.testSeqRefineMode():
            s = u' (sequential refinement)'
        else:
            s = u''
        self.SetTitle("GSAS-II project: "+projName + s)
        self.plotFrame.SetTitle("GSAS-II plots: "+projName)

    def OnFileSaveas(self, event):
        '''Save the current project in response to the
        File/Save as menu button
        '''
        if GSASIIpath.GetConfigValue('Starting_directory'):
            pth = GSASIIpath.GetConfigValue('Starting_directory')
            pth = os.path.expanduser(pth) 
        elif self.LastGPXdir:
            pth = self.LastGPXdir
        else:
            pth = '.'
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file name', pth, self.newGPXfile, 
            'GSAS-II project file (*.gpx)|*.gpx',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:     #TODO: what about Cancel?
                self.GSASprojectfile = dlg.GetPath()
                self.GSASprojectfile = G2IO.FileDlgFixExt(dlg,self.GSASprojectfile)
                self.GPXtree.SetItemText(self.root,'Project: '+self.GSASprojectfile)
                self.CheckNotebook()
                G2IO.ProjFileSave(self)
                self.SetTitleByGPX()
                os.chdir(dlg.GetDirectory())           # to get Mac/Linux to change directory!
                config = G2G.GetConfigValsDocs()
                GSASIIpath.addPrevGPX(self.GSASprojectfile,config)
                return True
            else:
                return False
        finally:
            dlg.Destroy()
            
    def ExpandAll(self,event):
        '''Expand all tree items or those of a single type
        '''
        txt = self.GetMenuBar().GetLabel(event.Id)
        if txt == 'all':
            self.ExpandingAll = True
            try:
                self.GPXtree.ExpandAll()
            finally:
                self.ExpandingAll = False
        else:
            self.ExpandingAll = True
            try:
                item, cookie = self.GPXtree.GetFirstChild(self.root)
                while item:
                    name = self.GPXtree.GetItemText(item)
                    if name.startswith(txt+' '): self.GPXtree.Expand(item)
                    item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
            finally:
                self.ExpandingAll = False

    def MoveTreeItems(self,event):
        '''Move tree items of a single type to the end of the tree
        '''
        txt = self.GetMenuBar().GetLabel(event.Id)
        # make a list of items to copy
        copyList = []
        item, cookie = self.GPXtree.GetFirstChild(self.root)
        while item:
            if self.GPXtree.GetItemText(item).startswith(txt+' '):
                copyList.append(item)
            item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        
        self.ExpandingAll = True
        try:
            for item in copyList:
                name = self.GPXtree.GetItemText(item)
                newId = self.GPXtree.AppendItem(self.root,name)
                self.GPXtree.SetItemPyData(newId,self.GPXtree.GetItemPyData(item))
                chld, chldcookie = self.GPXtree.GetFirstChild(item)
                while chld:
                    chname = self.GPXtree.GetItemText(chld)
                    newCh = self.GPXtree.AppendItem(newId,chname)
                    self.GPXtree.SetItemPyData(newCh,self.GPXtree.GetItemPyData(chld))
                    chld, chldcookie = self.GPXtree.GetNextChild(item, chldcookie)
                self.GPXtree.Delete(item)
        finally:
            self.ExpandingAll = False
        SelectDataTreeItem(self,self.root)
            
    def ExitMain(self, event):
        '''Called if exit selected or the main window is closed
        rescord last position of data & plot windows; saved to config.py file
        NB: not called if console window closed
        '''
        if self.GPXtree.GetCount() > 1:
            dlg = wx.MessageDialog(self,
                    'Do you want to save and exit?\n(Use No to exit without save or Cancel to prevent exiting)',
                    'Confirm exit/save?',
                    wx.YES|wx.NO|wx.CANCEL)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
        else:
            result = wx.ID_NO
        if result == wx.ID_NO:
            pass
        elif result == wx.ID_CANCEL:
            return
        else:
            if not self.OnFileSave(event): return
        FrameInfo = {'Main_Pos':tuple(self.GetPosition()),
                     'Main_Size':tuple(self.GetSize()),
                     'Plot_Pos':tuple(self.plotFrame.GetPosition()),
                     'Plot_Size':tuple(self.plotFrame.GetSize())}
        GSASIIpath.SetConfigValue(FrameInfo)
#        FramePos = {'Main_Pos':tuple(self.GetPosition()),'Plot_Pos':tuple(self.plotFrame.GetPosition())}
#        GSASIIpath.SetConfigValue(FramePos)
        config = G2G.GetConfigValsDocs()
        G2G.SaveConfigVars(config)
        if self.G2plotNB:
            self.G2plotNB.Destroy()
        if self.undofile:
            os.remove(self.undofile)
        sys.exit()
        
    def OnExportPeakList(self,event):
        pth = G2G.GetExportPath(self)
        dlg = wx.FileDialog(self, 'Choose output peak list file name', pth, '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.GPXtree.GetFirstChild(self.root)
                while item:
                    name = self.GPXtree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.GPXtree.GetFirstChild(item)
                        wave = 0.0
                        while item2:
                            name2 = self.GPXtree.GetItemText(item2)
                            if name2 == 'Instrument Parameters':
                                Inst = self.GPXtree.GetItemPyData(item2)[0]
                                Type = Inst['Type'][0]
                                if 'T' not in Type:
                                    wave = G2mth.getWave(Inst)
                            elif name2 == 'Peak List':
                                pkdata = self.GPXtree.GetItemPyData(item2)
                                peaks = pkdata['peaks']
                                sigDict = pkdata['sigDict']
                            item2, cookie2 = self.GPXtree.GetNextChild(item, cookie2)                            
                        file.write("#%s \n" % (name+' Peak List'))
                        if wave:
                            file.write('#wavelength = %10.6f\n'%(wave))
                        if 'T' in Type:
                            file.write('#%9s %10s %10s %12s %10s %10s %10s %10s %10s %10s\n'%('pos','dsp','esd','int','esd','alp','bet','sig','gam','FWHM'))                                    
                        else:
                            file.write('#%9s %10s %10s %12s %10s %10s %10s %10s\n'%('pos','dsp','esd','int','esd','sig','gam','FWHM'))
                        for ip,peak in enumerate(peaks):
                            dsp = G2lat.Pos2dsp(Inst,peak[0])
                            if 'T' in Type:  #TOF - more cols
                                esds = {'pos':0.,'int':0.,'alp':0.,'bet':0.,'sig':0.,'gam':0.}
                                for name in list(esds.keys()):
                                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                                sig = np.sqrt(peak[8])
                                gam = peak[10]
                                esddsp = abs(G2lat.Pos2dsp(Inst,peak[0]-esds['pos'])-G2lat.Pos2dsp(Inst,peak[0]+esds['pos']))/2.
                                FWHM = G2pwd.getgamFW(gam,sig) +(peak[4]+peak[6])*np.log(2.)/(peak[4]*peak[6])     #to get delta-TOF from Gam(peak)
                                file.write("%10.2f %10.5f %10.5f %12.2f%10.2f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % \
                                    (peak[0],dsp,esddsp,peak[2],esds['int'],peak[4],peak[6],peak[8],peak[10],FWHM))
                            else:               #CW
                                #get esds from sigDict for each peak & put in output - esds for sig & gam from UVWXY?
                                esds = {'pos':0.,'int':0.,'sig':0.,'gam':0.}
                                for name in list(esds.keys()):
                                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                                sig = np.sqrt(peak[4]) #var -> sig
                                gam = peak[6]
                                esddsp = abs(G2lat.Pos2dsp(Inst,peak[0]-esds['pos'])-G2lat.Pos2dsp(Inst,peak[0]+esds['pos']))/2.
                                FWHM = G2pwd.getgamFW(gam,sig)      #to get delta-2-theta in deg. from Gam(peak)
                                file.write("%10.4f %10.5f %10.5f %12.2f %10.2f %10.5f %10.5f %10.5f \n" % \
                                    (peak[0],dsp,esddsp,peak[2],esds['int'],np.sqrt(max(0.0001,peak[4]))/100.,peak[6]/100.,FWHM/100.)) #convert to deg
                    item, cookie = self.GPXtree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportHKL(self,event):
        pth = G2G.GetExportPath(self)
        dlg = wx.FileDialog(self, 'Choose output reflection list file name', pth, '', 
            '(*.*)|*.*',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.peaklistfile = dlg.GetPath()
                self.peaklistfile = G2IO.FileDlgFixExt(dlg,self.peaklistfile)
                file = open(self.peaklistfile,'w')                
                item, cookie = self.GPXtree.GetFirstChild(self.root)
                while item:
                    name = self.GPXtree.GetItemText(item)
                    if 'PWDR' in name:
                        item2, cookie2 = self.GPXtree.GetFirstChild(item)
                        while item2:
                            name2 = self.GPXtree.GetItemText(item2)
                            if name2 == 'Reflection Lists':
                                data = self.GPXtree.GetItemPyData(item2)
                                phases = data.keys()
                                for phase in phases:
                                    peaks = data[phase]
                                    I100 = peaks['RefList'].T[8]*np.array([refl[11] for refl in peaks['RefList']])
                                    Imax = np.max(I100)
                                    if Imax:
                                        I100 *= 100.0/Imax
                                    file.write("%s %s %s \n" % (name,phase,' Reflection List'))
                                    if 'T' in peaks.get('Type','PXC'):
                                        file.write('%s \n'%('   h   k   l   m   d-space       TOF       wid     Fo**2     Fc**2     Icorr      Prfo     Trans      ExtP      I100'))
                                    else:                
                                        file.write('%s \n'%('   h   k   l   m   d-space   2-theta       wid     Fo**2     Fc**2     Icorr      Prfo     Trans      ExtP      I100'))
                                    for ipk,peak in enumerate(peaks['RefList']):
                                        if 'T' in peaks.get('Type','PXC'):
                                            sig = np.sqrt(peak[6])
                                            gam = peak[7]
                                            FWHM = G2pwd.getgamFW(gam,sig)
                                            file.write(" %3d %3d %3d %3d%10.5f%10.2f%10.5f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" % \
                                                (int(peak[0]),int(peak[1]),int(peak[2]),int(peak[3]),peak[4],peak[5],FWHM,peak[8],
                                                peak[9],peak[11],peak[12],peak[13],peak[14],I100[ipk]))
                                        else:
                                            sig = np.sqrt(peak[6])
                                            gam = peak[7]
                                            FWHM = G2pwd.getgamFW(gam,sig)
                                            file.write(" %3d %3d %3d %3d%10.5f%10.5f%10.5f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" % \
                                                (int(peak[0]),int(peak[1]),int(peak[2]),int(peak[3]),peak[4],peak[5],FWHM/100.,
                                                peak[8],peak[9],peak[11],peak[12],peak[13],peak[14],I100[ipk]))
                            item2, cookie2 = self.GPXtree.GetNextChild(item, cookie2)                            
                    item, cookie = self.GPXtree.GetNextChild(self.root, cookie)                            
                file.close()
        finally:
            dlg.Destroy()
        
    def OnExportPDF(self,event):
        'Save S(Q), G(R),... as selected by user'
        def PDFSave(G2frame,exports,PDFsaves):
            'Save a PDF I(Q), S(Q), F(Q), G(r) and g(r)  in column formats'
            if len(exports) > 1:
                dirname = G2G.askSaveDirectory(G2frame)
                if not dirname: return
            else:
                defnam = exports[0].replace(' ','_')[5:]
                filename = G2G.askSaveFile(G2frame,defnam,'.gr','G(r) file, etc.')
                if not filename: return
                dirname,filename = os.path.split(filename)
                filename = os.path.splitext(filename)[0]
            Inst = None
            Limits = None
            for export in exports:
                if len(exports) > 1:
                    filename = export.replace(' ','_')[5:]
                PickId = GetGPXtreeItemId(G2frame, G2frame.root, export)
                PDFControls = G2frame.GPXtree.GetItemPyData(
                    GetGPXtreeItemId(G2frame, PickId,'PDF Controls'))
                if PDFsaves[4]: #pdfGUI file for G(R)
                    pId = GetGPXtreeItemId(G2frame, G2frame.root, 'PWDR'+export[4:])
                    Inst = G2frame.GPXtree.GetItemPyData(
                        GetGPXtreeItemId(G2frame, pId,'Instrument Parameters'))[0]
                    Limits = G2frame.GPXtree.GetItemPyData(
                        GetGPXtreeItemId(G2frame, pId,'Limits'))
                G2fil.PDFWrite(export,os.path.join(dirname,filename),
                            PDFsaves,PDFControls,Inst,Limits)
        
        names = G2pdG.GetFileList(self,'PDF')
        exports = []
        if names:
            od = {'label_1':'Export I(Q)','value_1':False,'label_2':'Export S(Q)','value_2':False,
                  'label_3':'Export F(Q)','value_3':False,'label_4':'Export G(R)','value_4':True,
                  'label_5':'Make G(R) for pdfGUI','value_5':False,
                  'label_6':'Make F(Q) & g(r) for RMCProfile','value_6':False}
            dlg = G2G.G2MultiChoiceDialog(self,'Select','PDF patterns to export',names,extraOpts=od)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for x in sel:
                    exports.append(names[x])
            dlg.Destroy()
        if exports:
            PDFsaves = [od['value_1'],od['value_2'],od['value_3'],od['value_4'],od['value_5'],od['value_6']]
            PDFSave(self,exports,PDFsaves)
        
    def OnMakePDFs(self,event):
        '''Sets up PDF data structure filled with defaults; if found chemical formula is inserted
        so a default PDF can be made.
        '''
        sind = lambda x: math.sin(x*math.pi/180.)
        tth2q = lambda t,w:4.0*math.pi*sind(t/2.0)/w
        tof2q = lambda t,C:2.0*math.pi*C/t
        TextList = []
        ElLists = []
        Qlimits = []
        Names = []
        if self.GPXtree.GetCount():
            Id, cookie = self.GPXtree.GetFirstChild(self.root)
            while Id:
                name = self.GPXtree.GetItemText(Id)
                Names.append(name)
                if 'PWDR' in name:
                    TextList.append(name)
                    Data = self.GPXtree.GetItemPyData(Id)[1]
                    pwdrMin = np.min(Data[1])
                    Comments = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,Id,'Comments'))
                    Parms = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,Id,'Instrument Parameters'))[0]
                    fullLimits = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,Id,'Limits'))[0]
                    if 'C' in Parms['Type'][0]:
                        wave = G2mth.getWave(Parms)
                        qMax = tth2q(fullLimits[1],wave)
                    else:   #'T'of
                        qMax = tof2q(fullLimits[0],Parms['difC'][1])
                    Qlimits.append([0.9*qMax,qMax])
                    ElList = {}
                    sumnum = 1.
                    for item in Comments:           #grab chemical formula from Comments
                        if 'formula' in item[:15].lower():
                            formula = item.split('=')[1].strip('"\n').split()
                            try:
                                if len(formula) < 2:
                                    raise ValueError
                                elems = formula[::2]
                                nums = formula[1::2]
                                Formula = zip(elems,nums)
                                sumnum = 0.
                                for [elem,num] in Formula:
                                    ElData = G2elem.GetElInfo(elem,Parms)
                                    ElData['FormulaNo'] = float(num)
                                    sumnum += float(num)
                                    ElList[elem] = ElData
                                
                            except ValueError:
                                G2G.G2MessageBox(self,'Carbon-based (and wrong) PDF will be generated','Missing chemical formula')
                                ElData = G2elem.GetElInfo('C',Parms)
                                sumnum = 1.0
                                ElData['FormulaNo'] = 1.0
                                ElList['C'] = ElData
                    ElLists.append(ElList)
                Id, cookie = self.GPXtree.GetNextChild(self.root, cookie)
            if len(TextList) < 1:
                self.ErrorDialog('Nothing to make PDFs for','There must be at least one "PWDR" pattern')
                return
            dlg = G2G.G2MultiChoiceDialog(self,'Make PDF controls','Make PDF controls for:',TextList, wx.CHOICEDLG_STYLE)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    for i in dlg.GetSelections():
                        PDFnames = GetGPXtreeDataNames(self,['PDF ',])
                        G2obj.CreatePDFitems(self,TextList[i],ElLists[i],Qlimits[i],sumnum,pwdrMin,PDFnames)
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
            PWDRdata.update(self.GPXtree.GetItemPyData(PWDRname)[0])            #wtFactor + ?
        except ValueError:
            PWDRdata['wtFactor'] = 1.0
        PWDRdata['Data'] = self.GPXtree.GetItemPyData(PWDRname)[1]          #powder data arrays
        PWDRdata['Limits'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,PWDRname,'Limits'))
        PWDRdata['Background'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,PWDRname,'Background'))
        PWDRdata['Instrument Parameters'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,PWDRname,'Instrument Parameters'))
        PWDRdata['Sample Parameters'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,PWDRname,'Sample Parameters'))
        PWDRdata['Reflection Lists'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,PWDRname,'Reflection Lists'))
        if 'ranId' not in PWDRdata:  # patch, add a random Id
            PWDRdata['ranId'] = ran.randint(0,sys.maxsize)
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
        HKLFdata.update(self.GPXtree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        try:
#            HKLFdata.update(self.GPXtree.GetItemPyData(HKLFname)[0])            #wtFactor + ?
#        except ValueError:
#            HKLFdata['wtFactor'] = 1.0
        HKLFdata['Data'] = self.GPXtree.GetItemPyData(HKLFname)[1]
        HKLFdata['Instrument Parameters'] = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,HKLFname,'Instrument Parameters'))
        return HKLFdata
        
    def GetPhaseData(self):
        '''Returns a dict with defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseData` also exists to
        get same info from GPX file.
        '''
        phaseData = {}
        if GetGPXtreeItemId(self,self.root,'Phases'):
            sub = GetGPXtreeItemId(self,self.root,'Phases')
        else:
#            print 'no phases found in GetPhaseData'
            sub = None
        if sub:
            item, cookie = self.GPXtree.GetFirstChild(sub)
            while item:
                phaseName = self.GPXtree.GetItemText(item)
                phaseData[phaseName] =  self.GPXtree.GetItemPyData(item)
                if 'ranId' not in phaseData[phaseName]:
                    phaseData[phaseName]['ranId'] = ran.randint(0,sys.maxsize)          
                item, cookie = self.GPXtree.GetNextChild(sub, cookie)
        return phaseData

    def GetPhaseInfofromTree(self):
        '''Get the phase names and their rId values,
        also the histograms used in each phase. 

        :returns: (phaseRIdList, usedHistograms) where 

          * phaseRIdList is a list of random Id values for each phase
          * usedHistograms is a dict where the keys are the phase names
            and the values for each key are a list of the histogram names
            used in each phase. 
        '''
        phaseRIdList = []
        usedHistograms = {}
        sub = GetGPXtreeItemId(self,self.root,'Phases')
        if sub:
            item, cookie = self.GPXtree.GetFirstChild(sub)
            while item:
                phaseName = self.GPXtree.GetItemText(item)
                ranId = self.GPXtree.GetItemPyData(item).get('ranId')
                if ranId: phaseRIdList.append(ranId)
                data = self.GPXtree.GetItemPyData(item)
                UseList = data['Histograms']
                usedHistograms[phaseName] = list(UseList.keys())
                item, cookie = self.GPXtree.GetNextChild(sub, cookie)
        return phaseRIdList,usedHistograms

    def GetPhaseNames(self):
        '''Returns a list of defined phases. 
        Note routine :func:`GSASIIstrIO.GetPhaseNames` also exists to
        get same info from GPX file.
        '''
        phaseNames = []
        if GetGPXtreeItemId(self,self.root,'Phases'):
            sub = GetGPXtreeItemId(self,self.root,'Phases')
        else:
#            print 'no phases found in GetPhaseNames'
            sub = None
        if sub:
            item, cookie = self.GPXtree.GetFirstChild(sub)
            while item:
                phase = self.GPXtree.GetItemText(item)
                phaseNames.append(phase)
                item, cookie = self.GPXtree.GetNextChild(sub, cookie)
        return phaseNames
    
    def GetHistogramNames(self,hType):
        """ Returns a list of histogram names found in the GSASII data tree
        Note routine :func:`GSASIIstrIO.GetHistogramNames` also exists to
        get same info from GPX file.
        
        :param str hType: list of histogram types
        :return: list of histogram names
        
        """
        HistogramNames = []
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name[:4] in hType:
                    HistogramNames.append(name)        
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        return HistogramNames
    
    def GetHistogramNamesID(self,hType):
        """ Returns a list of histogram names found in the GSASII data tree
        and a lookup table of their Id values. Should replace GetHistogramNames
        since that will not be much faster (and there may be real speed gains from
        caching the Ids rather than keep searching for them). 

        N.B routine :func:`GSASIIstrIO.GetHistogramNames` also exists to
        get same info, but from GPX file.
        
        :param str hType: list of histogram types
        :return: list of histogram names and a dict of histogram Ids 
           keyed by histogram name.
        """
        HistogramNames = []
        HistogramIds = {}
        if self.GPXtree.GetCount():
            item, cookie = self.GPXtree.GetFirstChild(self.root)
            while item:
                name = self.GPXtree.GetItemText(item)
                if name[:4] in hType:
                    HistogramNames.append(name)        
                    HistogramIds[name] = item
                item, cookie = self.GPXtree.GetNextChild(self.root, cookie)
        return HistogramNames,HistogramIds
                    
    def GetUsedHistogramsAndPhasesfromTree(self):
        ''' Returns all histograms that are found in any phase
        and any phase that uses a histogram.
        This also assigns numbers to used phases and histograms by the
        order they appear in the file. 
        Note routine :func:`GSASIIstrIO.GetUsedHistogramsAndPhases` also exists to
        get same info from GPX file.

        :returns: (Histograms,Phases)

            * Histograms = dictionary of histograms as {name:data,...}
            * Phases = dictionary of phases that use histograms
        '''
        Histograms = {}
        Phases = {}
        phaseNames = self.GetPhaseNames()
        phaseData = self.GetPhaseData()
        histoList,histIdList = self.GetHistogramNamesID(['PWDR','HKLF'])

        badnum = 0
        for phase in phaseData:
            Phase = phaseData[phase]
            pId = phaseNames.index(phase)
            Phase['pId'] = pId
            if Phase['Histograms']:
                for hist in Phase['Histograms']:
                    if 'Use' not in Phase['Histograms'][hist]:      #patch: add Use flag as True
                        Phase['Histograms'][hist]['Use'] = True
                    if Phase['Histograms'][hist]['Use'] and phase not in Phases:
                        Phases[phase] = Phase
                    if (hist not in Histograms) and Phase['Histograms'][hist]['Use']:
                        if hist not in histIdList:
                            if badnum == 0:
                                print('Error: hist {} not found in histIdList. Deleted?'.format(hist))
                            badnum += 1
                            continue
                        item = histIdList[hist]
                        if item:
                            if 'PWDR' in hist[:4]: 
                                Histograms[hist] = self.GetPWDRdatafromTree(item)
                            elif 'HKLF' in hist[:4]:
                                Histograms[hist] = self.GetHKLFdatafromTree(item)
                            hId = histoList.index(hist)
                            Histograms[hist]['hId'] = hId
                        else: # would happen if a referenced histogram were renamed or deleted
                            print(u'For phase "'+phase+
                                  u'" unresolved reference to histogram "'+hist+u'"')
        if badnum > 1: print('  ...hist not in histIdList error occured {} times'.format(badnum))
        G2obj.IndexAllIds(Histograms=Histograms,Phases=phaseData)
        return Histograms,Phases
        
    def MakeLSParmDict(self,seqHist=None):
        '''Load all parameters used for computation from the tree into a
        dict of paired values [value, refine flag]. Note that this is
        different than the parmDict used in the refinement, which only has
        values.

        Note that similar things are done in
        :meth:`GSASIIIO.ExportBaseclass.loadParmDict` (from the tree) and 
        :func:`GSASIIstrMain.Refine` and :func:`GSASIIstrMain.SeqRefine` (from
        a GPX file).

        :param dict seqHist: defines a specific histogram to be loaded for a sequential
          refinement, if None (default) all are loaded.
          Note: at present this parameter is not used anywhere.

        :returns: (parmDict,varyList) where:

         * parmDict is a dict with values and refinement flags
           for each parameter and
         * varyList is a list of variables (refined parameters).
        '''
        parmDict = {}
        Histograms,Phases = self.GetUsedHistogramsAndPhasesfromTree()
        if seqHist:
            histDict = {seqHist:Histograms[seqHist]}
        else:
            histDict = Histograms
        for phase in Phases:
            if 'pId' not in Phases[phase]:
                self.ErrorDialog('View parameter error','You must run least squares at least once')
                raise Exception('No pId for phase '+phase)
        rigidbodyDict = self.GPXtree.GetItemPyData(   
            GetGPXtreeItemId(self,self.root,'Rigid bodies'))
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtable,BLtable,MFtable,maxSSwave = G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)        
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,histDict,Print=False,resetRefList=False)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(histDict,Print=False)
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
        # for i in parmDict: print i,'\t',parmDict[i]
        # fl = open('parmDict.dat','wb')
        # cPickle.dump(parmDict,fl,1)
        # fl.close()
        return parmDict,varyList

    def OnShowLSParms(self,event):
        '''Displays a window showing all parameters in the refinement.
        Called from the Calculate/View LS Parms menu.

        This could potentially be sped up by loading only the histogram that is needed
        for a sequential fit. 
        '''
        G2cnstG.CheckAllScalePhaseFractions(self)
        try:
            parmDict,varyList = self.MakeLSParmDict()
        except:
            print('Error retrieving parameters')
            return
        parmValDict = {}
        for i in parmDict:
            parmValDict[i] = parmDict[i][0]
            
        reqVaryList = tuple(varyList) # save requested variables
        wx.BeginBusyCursor()
        try:
            # process constraints
            sub = GetGPXtreeItemId(self,self.root,'Constraints') 
            Constraints = self.GPXtree.GetItemPyData(sub)
            constList = []
            for item in Constraints:
                if item.startswith('_'): continue
                constList += Constraints[item]
            G2mv.InitVars()
            constrDict,fixedList,ignored = G2stIO.ProcessConstraints(constList)
            msg = G2mv.EvaluateMultipliers(constrDict,parmValDict)
            if msg:
                print('Unable to interpret multiplier(s):',msg)
            else:
                G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmValDict)
                G2mv.Map2Dict(parmValDict,varyList)
        except G2mv.ConstraintException:
            pass
        Controls = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,self.root, 'Controls'))
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        wx.EndBusyCursor()
        # debug stuff
        #if GSASIIpath.GetConfigValue('debug'):
        #    print('reloading',G2G)
        #    import imp
        #    imp.reload(G2G)
        # end debug stuff            
        dlg = G2G.ShowLSParms(self,'Least Squares Parameters',parmValDict,
                    varyList,reqVaryList,Controls)
        dlg.ShowModal()
        dlg.Destroy()

    def OnRefine(self,event):
        '''Perform a refinement or a sequential refinement (depending on controls setting)
        Called from the Calculate/Refine menu.
        '''
        if self.testSeqRefineMode():
            self.OnSeqRefine(event)
            return
        G2cnstG.CheckAllScalePhaseFractions(self) # can be slow for sequential fits, skip
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.GSASprojectfile)
        if errmsg:
            self.ErrorDialog('Refinement error',errmsg+'\nCheck console output for more information')
            return
        if warnmsg:
            print('Conflict between refinement flag settings and constraints:\n'+
                warnmsg+'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                'Conflict between refinement flag settings and constraints:\n'+
                warnmsg+'\nRefinement not possible')
            return
        dlg = wx.ProgressDialog('Residual','All data Rw =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.STAY_ON_TOP,parent=self)
        Size = dlg.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            dlg.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        dlg.CenterOnParent()
        dlg.Raise()
        Rw = 100.00
        self.SaveTreeSetting() # save the current tree selection
        self.GPXtree.SaveExposedItems()             # save the exposed/hidden tree items
        if self.PatternId and self.GPXtree.GetItemText(self.PatternId).startswith('PWDR '):
            refPlotUpdate = G2plt.PlotPatterns(self,refineMode=True) # prepare for plot updating
        else:
            refPlotUpdate = None
        try:
            OK,Msg = G2stMn.Refine(self.GSASprojectfile,dlg,refPlotUpdate=refPlotUpdate)    #Msg is Rvals dict if Ok=True
        finally:
            dlg.Update(101.) # forces the Auto_Hide; needed after move w/Win & wx3.0
            dlg.Destroy()
#            wx.Yield()
        if OK:
            Rw = Msg['Rwp']
            rtext = 'LS Refinement: Rw = %.3f%%, GOF = %.2f, Nobs = %d, Max delt/sig = %.3f'%(Msg['Rwp'],Msg['GOF'],Msg['Nobs'],Msg['Max shft/sig'])
            lamMax = Msg.get('lamMax',0.001)
            lst = os.path.splitext(os.path.abspath(self.GSASprojectfile))[0]
            text = 'Detailed results are in ' + lst + '.lst\n'
            if 'GOF0' in Msg and 'GOF' in Msg:
                text += '\nReduced Chi^2 before refinement={:.3f} and after={:.3f}\n'.format(
                    Msg['GOF0']**2,Msg['GOF']**2)
            if 'Max shft/sig' in Msg:
                text += '\nMax shift/sigma={:.3f}\n'.format(Msg['Max shft/sig'])
            if 'msg' in Msg: text += '\n' + Msg['msg'] + '\n'
            text += '\nLoad new result?'
            if lamMax >= 10.:
                text += '\nWARNING: Steepest descents dominates;'+   \
                ' minimum may not have been reached\nor result may be false minimum.'+  \
                ' You should reconsider your parameter suite'
            dlg2 = wx.MessageDialog(self,text,'Refinement results, Rw =%.3f'%(Rw),wx.OK|wx.CANCEL)
            try:
                if dlg2.ShowModal() == wx.ID_OK:
                    if refPlotUpdate: refPlotUpdate({},restore=True)
                    wx.CallAfter(self.reloadFromGPX,rtext)
                else:
                    if refPlotUpdate: refPlotUpdate({},restore=True)
            finally:
                dlg2.Destroy()
        else:
            self.ErrorDialog('Refinement error',Msg)
            
    def reloadFromGPX(self,rtext=None):
        '''Deletes current data tree & reloads it from GPX file (after a 
        refinemnt.) Done after events are completed to avoid crashes.
        :param rtext str: string info from cller to be put in Notebook after reload
        '''
        self.GPXtree.DeleteChildren(self.root)
        self.HKL = []
        G2IO.ProjFileOpen(self,False)
        self.TreeItemDelete = False  # tree has been repopulated; ignore previous deletions
        self.GPXtree.RestoreExposedItems() # reset exposed/hidden tree items
        if rtext is not None:
            self.AddToNotebook(rtext)
        self.ResetPlots()        
        
    def SaveTreeSetting(self):
        'Save the current selected tree item by name (since the id will change)'
        oldId =  self.GPXtree.GetSelection()        #retain current selection
        oldPath = self.GetTreeItemsList(oldId)
        self.lastTreeSetting = oldPath
        
    def ResetPlots(self):
        '''This reloads the current tree item, often drawing a plot. It 
        also refreshes any plots that have registered a refresh routine 
        (see G2plotNB.RegisterRedrawRoutine) and deletes all plots that 
        have not been refreshed and require one (see G2plotNB.SetNoDelete).
        '''
        for lbl,win in zip(self.G2plotNB.plotList,self.G2plotNB.panelList):
            win.plotInvalid = True  # mark all current plots as invalid so we can tell what has been updated
            
        oldPath = self.lastTreeSetting  # reload last selected tree item, triggers window and possibly plot redraw
        Id = self.root
        for txt in oldPath:
            Id = GetGPXtreeItemId(self, Id, txt)
        pltText = None
        self.PickIdText = None   # forces reload of page when selected
        if Id:
            self.PickId = Id
            self.GPXtree.SelectItem(Id)
            # selection above appears to trigger a tree event (all platforms?),
            # but if not call SelectDataTreeItem(self,Id)
#            wx.Yield()
            time.sleep(0.1)
            pltNumber = self.G2plotNB.nb.GetSelection()
            if pltNumber >= 0:
                pltText = self.G2plotNB.nb.GetPageText(pltNumber)
            else:
                pltText = None
        # update plots where a routine is supplied
        for lbl,win in zip(self.G2plotNB.plotList,self.G2plotNB.panelList):
            if win.plotInvalid and win.replotFunction:
                if GSASIIpath.GetConfigValue('debug'):
                    print('updating',lbl,'by calling',str(win.replotFunction))
                try:
                    win.replotFunction(*win.replotArgs,**win.replotKWargs)
                except Exception as msg:
                    if GSASIIpath.GetConfigValue('debug'):
                        print('Error calling',win.replotFunction,'with args',
                                  win.replotArgs,win.replotKWargs)
                    if GSASIIpath.GetConfigValue('debug'):
                        print(msg)
                        GSASIIpath.IPyBreak()
        # delete any remaining plots unless they have been
        # updated (win.plotInvalid is False) or are tagged as not
        # needing a refresh (win.plotRequiresRedraw is False)
        for lbl,win in zip(self.G2plotNB.plotList,self.G2plotNB.panelList):
            if win.plotInvalid and win.plotRequiresRedraw:
                if GSASIIpath.GetConfigValue('debug'):
                    print('Closing out-of-date plot',lbl)
                if lbl == 'Powder Patterns':
                    self.lastPlotType = None
                self.G2plotNB.Delete(lbl)
        # put the previously last-raised plot tab on top, if present.
        # Search by label text, since tab number may have changed
        for i in range(self.G2plotNB.nb.GetPageCount()):
            if self.G2plotNB.nb.GetPageText(i) == pltText:
                self.G2plotNB.nb.SetSelection(i)
                break
        
    def OnSeqRefine(self,event):
        '''Perform a sequential refinement.
        Called from self.OnRefine (Which is called from the Calculate/Refine menu)
        '''
        seqList = self.testSeqRefineMode()
        if not seqList:
            self.OnRefine(event)
            return
        #plotHist = self.GPXtree.GetItemText(self.PatternId)
        Id = GetGPXtreeItemId(self,self.root,'Sequential results')
        if not Id:
            Id = self.GPXtree.AppendItem(self.root,text='Sequential results')
            self.GPXtree.SetItemPyData(Id,{})            
        self.G2plotNB.Delete('Sequential refinement')    #clear away probably invalid plot
        Controls = self.GPXtree.GetItemPyData(GetGPXtreeItemId(self,self.root, 'Controls'))
        Controls['ShowCell'] = True
        self.OnFileSave(event)
        # check that constraints are OK here
        errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.GSASprojectfile)
        #errmsg, warnmsg = G2stIO.ReadCheckConstraints(self.GSASprojectfile,seqList[0]) # this would be faster, but at present it might not catch all errors
        if errmsg:
            self.ErrorDialog('Refinement error',errmsg)
            return
        if warnmsg:
            print(u'Conflict between refinment flag settings and constraints:\n'+
                  warnmsg+u'\nRefinement not possible')
            self.ErrorDialog('Refinement Flag Error',
                             u'Conflict between refinment flag settings and constraints:\n'+
                             warnmsg+u'\nRefinement not possible')
            return
        self.GPXtree.SaveExposedItems()        
        dlg = wx.ProgressDialog('Residual for histogram 0','Powder profile Rwp =',101.0, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT,
            parent=self)            
        Size = dlg.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            dlg.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        dlg.CenterOnParent()
        dlg.Show()
        # find 1st histogram to be refined
        if 'Seq Data' in Controls: 
            histNames = Controls['Seq Data']
        else: # patch from before Controls['Seq Data'] was implemented
            histNames = G2stIO.GetHistogramNames(self.GSASprojectfile,['PWDR',])
        if Controls.get('Reverse Seq'):
            histNames.reverse()
        # select it
        self.PatternId = GetGPXtreeItemId(self,self.root,histNames[0])
        if self.GPXtree.GetItemText(self.PatternId).startswith('PWDR '):
            refPlotUpdate = G2plt.PlotPatterns(self,refineMode=True) # prepare for plot updating
        else:
            refPlotUpdate = None
        try:
            OK,Msg = G2stMn.SeqRefine(self.GSASprojectfile,dlg,refPlotUpdate) #Msg is Rvals dict if Ok=True
        finally:
            dlg.Update(101.) # forces the Auto_Hide; needed after move w/Win & wx3.0
            dlg.Destroy()
            wx.Yield()
        if OK:
            lst = os.path.splitext(os.path.abspath(self.GSASprojectfile))[0]
            text = 'Detailed results are in ' + lst + '.lst\n'
            if Msg.get('Frozen'):
                text += '\n' +  Msg['Frozen']
            if Msg.get('steepestNum',0) > 0:
                text += '\nNote that {} histograms had extreme correlations where steepest descents dominates\n'.format(Msg['steepestNum'])
            if len(Msg.get('maxshift/sigma',[])) > 0:
                avg = np.average(Msg['maxshift/sigma'])
                mx = np.max(Msg['maxshift/sigma'])
                text += '\nBiggest Max shft/sig was {:.3f} (average across histograms {:.3f})\n'.format(mx,avg)
            text += '\nLoad new result?'                
            dlg = wx.MessageDialog(self,text,'Refinement results',wx.OK|wx.CANCEL)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    if refPlotUpdate: refPlotUpdate({},restore=True)
                    self.PickIdText = None  #force reload of PickId contents
                    self.GPXtree.DeleteChildren(self.root)
                    if len(self.HKL): self.HKL = []
                    if self.G2plotNB.plotList:
                        self.G2plotNB.clear()
                    G2IO.ProjFileOpen(self,False)
                    self.GPXtree.RestoreExposedItems()
                    self.ResetPlots()
                    #self.PatternId = GetGPXtreeItemId(self,self.root,plotHist)
                    #SelectDataTreeItem(self,self.PatternId)
                    sId = GetGPXtreeItemId(self,self.root,'Sequential results')
                    SelectDataTreeItem(self,sId)
                    self.GPXtree.SelectItem(sId)
                else:
                    if refPlotUpdate: refPlotUpdate({},restore=True)
            finally:
                dlg.Destroy()
            
#            self.SeqTblHideList = []
        else:
            self.ErrorDialog('Sequential refinement error',Msg)
            
    def OnRunFprime(self,event):
        import fprime
        self.fprime = fprime.Fprime(self)
        self.fprime.Show()
        
    def OnRunAbsorb(self,event):
        import Absorb
        self.absorb = Absorb.Absorb(self)
        self.absorb.Show()
        
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
    
    def OnSaveMultipleImg(self,event):
        '''Select and save multiple image parameter and mask files
        '''
        G2IO.SaveMultipleImg(self)

################################################################################
# Data window side of main GUI
################################################################################
class G2DataWindow(wx.ScrolledWindow):      #wxscroll.ScrolledPanel):
    '''Create the data item window as well as the menu. Note that 
    the same core menu items are used in all menus, but different items may be
    added depending on what data tree item (and for phases, the phase tab).

    Note that while the menus are created here, 
    the binding for the menus is done later in various GSASII*GUI modules,
    where the functions to be called are defined.

    Use of the dataWindow scrolled panel:
    
    dataWindow has a “master” vertical BoxSizer: find it with
    G2frame.dataWindow.GetSizer() and always use it. A call to
    dataWindow.SetSizer() should not be needed. 

    When placing a widget in the sizer that has its own scrolling
    (e.g. G2G.GSNoteBook, anything else?) that one widget should be placed
    in the sizer as

         G2frame.dataWindow.GetSizer().Add(G2frame.<obj>,1,wx.ALL|wx.EXPAND)

    [is wx.ALL superfluous here?] so that it consumes the full size of the
    panel and so that the NoteBook widget does the scrolling.

    For other uses, one will likely place a bunch of widgets and (other
    [sub-]sizers) into the master sizer. In this case, DO NOT use wx.EXPAND,
    as this will result in the widget resizing/repositioning as the window
    resizes. Possible exceptions might be for widgets going into a fixed-size
    panel that is inside the dataWindow (probably not being done). A call to
    Sizer.Fit(dataWindow) will do bad things, though a call to
    SubSizer.Fit(dataWindow.subpanel) could make sense.

    Initial GUI draws to dataWindow will go through
    GSASIIdataGUI.SelectDataTreeItem(), which is called after any changes to
    data tree selection. SelectDataTreeItem places items in dataWindow or
    calls that do that. Before it calls those routines, it calls 

        G2frame.dataWindow.ClearData()

    which deletes the contents of the master sizer. After the contents are
    posted a call is made to 

        G2frame.dataWindow.SetDataSize()

    which repaints the window. For routines [such as GSASIIpwdGUI.UpdatePeakGrid()]
    that are called repeatedly to update the entire contents of dataWindow
    themselves, it is important to add calls to 

        G2frame.dataWindow.ClearData()

    and

    	 G2frame.dataWindow.SetDataSize()

    at the beginning and end respectively to clear and refresh. This is not
    needed for GSNoteBook repaints, which seem to be working mostly
    automatically. If there is a problem, a call like 

         wx.CallAfter(G2frame.phaseDisplay.SendSizeEvent)

    might be needed. There are some calls to G2frame.dataWindow.SendSizeEvent()
    that may be doing the same thing. 
    '''

    def __init__(self,parent):
#        wxscroll.ScrolledPanel.__init__(self,parent,wx.ID_ANY,size=parent.GetSize())
        wx.ScrolledWindow.__init__(self,parent,wx.ID_ANY,size=parent.GetSize())
        self.parent = parent
        self._initMenus()
        self.currentGrids = []
        self.helpKey = ''  # defines help entry for current item selected in data tree
        
    def ClearData(self):
        '''Initializes the contents of the dataWindow panel
        '''
        self.Unbind(wx.EVT_SIZE)
        #self.SetBackgroundColour(wx.WHITE)
        self.SetBackgroundColour(VERY_LIGHT_GREY)  # BHT: I prefer a gray background. Makes TextCtrls stand out, but
        # a bit lighter than the splitter bar
        Sizer = self.GetSizer()
        if Sizer:
            try:
                Sizer.Clear(True)
            except:
                pass
        else:
            print ('No sizer in dataWindow')
            if GSASIIpath.GetConfigValue('debug'): raise Exception

    def OnResize(self,event):
        'Used for grids to match ScrolledWindow size'
        event.Skip()
        Sizer = self.GetSizer()
        if not Sizer: return
        #if Sizer.GetItemCount() == 1: # not wx 2.8
        if len(Sizer.GetChildren()) == 1: # if there is a single grid, resize it
            if isinstance(Sizer.GetItem(0).GetWindow(), G2G.GSGrid): 
                Sizer.GetItem(0).GetWindow().SetSize(self.GetSize())
                    
    def SetDataSize(self):
        '''Sizes the contents of the dataWindow panel
        '''
        Sizer = self.GetSizer()
        if not Sizer:
            print ('No sizer in dataWindow')
            if GSASIIpath.GetConfigValue('debug'): raise Exception
            return
        #if Sizer.GetItemCount() == 1: # not wx 2.8
        if len(Sizer.GetChildren()) == 1: # handle cases with a single grid in DataWindow differently
            # note that Grid's scroll bars must be turned on with .SetScrollRate(10,10)
            # just after the call to .GSGrid()
            if isinstance(Sizer.GetItem(0).GetWindow(), G2G.GSGrid):
                self.Bind(wx.EVT_SIZE,self.OnResize)
                Sizer.GetItem(0).GetWindow().SetSize(self.GetSize())
                self.SetAutoLayout(False)
                self.SetScrollRate(0,0) 
                self.SendSizeEvent()
                return
        self.SetAutoLayout(True)
        self.SetScrollRate(10,10)
        self.SendSizeEvent()

    def PrefillDataMenu(self,menu,empty=False):
        '''Create the "standard" part of data frame menus & add the dataWindow menu headings
        This menu duplicates the tree menu, but adds an extra help command for the current
        data item and a separator. 
        '''
        self.datamenu = menu
        self.GetTopLevelParent().FillMainMenu(menu,addhelp=False) # add the data tree menu items to the main menu
        if not empty:
            menu.Append(wx.Menu(title=''),title='|') # add a separator
        
    def PostfillDataMenu(self,empty=False):
        '''Add the help menu to the menus associated with data tree items.
        '''
        menu = self.datamenu
        G2frame = self.GetTopLevelParent()
        if not empty:
            menu.Append(wx.Menu(title=''),title='|') # add another separator
        HelpMenu=G2G.MyHelp(G2frame,includeTree=True,
            morehelpitems=[('&Tutorials','Tutorials'),])
        menu.Append(menu=HelpMenu,title='&Help')

    def _initMenus(self):
        '''define all GSAS-II data window menus. 
        NB: argument order conforms to both classic & phoenix variants for wx.
        Do not use argument= for these as the argument names are different for classic & phoenix
        '''
        
#        G2G.Define_wxId('wxID_MCRON', 'wxID_MCRLIST', 'wxID_MCRSAVE', 'wxID_MCRPLAY',)
##### GSAS-II Menu items
        # Main menu
        G2frame = self.GetTopLevelParent()
        G2frame.GSASIIMenu = wx.MenuBar()
        G2frame.dataMenuBars = [G2frame.GSASIIMenu] # list of all menus (added to in PrefillDataMenu)
        G2frame.FillMainMenu(G2frame.GSASIIMenu)
        G2frame.SetMenuBar(G2frame.GSASIIMenu)

        # Controls
        self.ControlsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ControlsMenu,empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Notebook
        self.DataNotebookMenu = wx.MenuBar() 
        self.PrefillDataMenu(self.DataNotebookMenu,empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Comments
        self.DataCommentsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataCommentsMenu,empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Constraints
        G2G.Define_wxId('wxID_CONSTRAINTADD', 'wxID_EQUIVADD', 'wxID_HOLDADD', 'wxID_FUNCTADD',
                        'wxID_ADDRIDING', 'wxID_CONSPHASE', 'wxID_CONSHIST', 'wxID_CONSHAP',
                        'wxID_CONSGLOBAL', 'wxID_CONSSYM', 'wxID_EQUIVALANCEATOMS',)
        self.ConstraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ConstraintMenu)
        self.ConstraintTab = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintTab, title='Select tab')
        for Id,txt in (
                (G2G.wxID_CONSPHASE,'Phase'),
                (G2G.wxID_CONSHAP,'Histogram/Phase'),
                (G2G.wxID_CONSHIST,'Histogram'),
                (G2G.wxID_CONSGLOBAL,'Global'),
                (G2G.wxID_CONSSYM,'Sym-Generated'),
                ):
            self.ConstraintTab.Append(Id,txt,'Select '+txt+' constraint editing tab')
        self.ConstraintEdit = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintEdit, title='Edit Constr.') # renamed from Edit due to Mac adding extra items to menu
        self.ConstraintEdit.Append(G2G.wxID_HOLDADD,'Add hold','Prevent refinement of parameter values')
        self.ConstraintEdit.Append(G2G.wxID_EQUIVADD,'Add equivalence','Force parameter values to be equivalent')
        self.ConstraintEdit.Append(G2G.wxID_CONSTRAINTADD,'Add constraint equation',
            'Add a constraint equation to apply to parameter values')
        self.ConstraintEdit.Append(G2G.wxID_FUNCTADD,'Add New Var',
            'Create a variable composed of existing parameters')
        self.ConstraintEdit.Append(G2G.wxID_EQUIVALANCEATOMS,'Make atoms equivalent',
            'Force atom parameter values to be equivalent')
        self.ConstraintEdit.Enable(G2G.wxID_EQUIVALANCEATOMS,False)
#        self.ConstraintEdit.Append(id=G2G.wxID_ADDRIDING, kind=wx.ITEM_NORMAL,text='Add H riding constraints',
#            help='Add H atom riding constraints between atom parameter values')
#        self.ConstraintEdit.Enable(G2G.wxID_ADDRIDING,False)
        G2G.Define_wxId('wxID_SHOWISO')
        self.ConstraintEdit.Append(G2G.wxID_SHOWISO,'Show ISODISTORT modes',
                'Show ISODISTORT mode values for all phases')
        self.ConstraintEdit.Enable(G2G.wxID_SHOWISO,False)
        # for DEBUG only
        #def OnShowISODISTORT(event):
        #    import imp
        #    imp.reload(G2cnstG)  # for testing changes to GSASIIconstrGUI
        #    G2cnstG.ShowIsoDistortCalc(G2frame)
        #G2frame.Bind(wx.EVT_MENU, OnShowISODISTORT, id=G2G.wxID_SHOWISO)
        # end DEBUG
        
        self.PostfillDataMenu()

        # Rigid bodies
        G2G.Define_wxId('wxID_RIGIDBODYADD', 'wxID_DRAWDEFINERB', 'wxID_RIGIDBODYIMPORT', 'wxID_RESIDUETORSSEQ',
            'wxID_VECTORBODYADD', 'wxID_RIGIDBODYSAVE',)        
        self.RigidBodyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodyMenu)
        self.ResidueRBMenu = wx.Menu(title='')
        self.ResidueRBMenu.Append(G2G.wxID_RIGIDBODYIMPORT,'Import XYZ','Import rigid body XYZ from file')
        G2G.Define_wxId('wxID_RIGIDBODYIMP')
        self.ResidueRBMenu.Append(G2G.wxID_RIGIDBODYIMP,'Extract from file',
                                'Extract rigid body from phase file')
        self.ResidueRBMenu.Append(G2G.wxID_RIGIDBODYSAVE,'Save as PDB','Save rigid body to PDB file')        
        self.ResidueRBMenu.Append(G2G.wxID_RESIDUETORSSEQ,'Define torsion','Define torsion sequence')
        self.ResidueRBMenu.Append(G2G.wxID_RIGIDBODYADD,'Import residues','Import residue rigid bodies from macro file')
        G2G.Define_wxId('wxID_RESBODYSAV')
        self.ResidueRBMenu.Append(G2G.wxID_RESBODYSAV,'Save rigid body','Write a rigid body to a file')
        G2G.Define_wxId('wxID_RESBODYRD')
        self.ResidueRBMenu.Append(G2G.wxID_RESBODYRD,'Read rigid body','Read a rigid body from a file')
        self.RigidBodyMenu.Append(menu=self.ResidueRBMenu, title='Edit Residue Body')
        self.PostfillDataMenu()

        self.VectorBodyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.VectorBodyMenu)
        self.VectorRBEdit = wx.Menu(title='')
        self.VectorRBEdit.Append(G2G.wxID_VECTORBODYADD,'Add rigid body','Add vector rigid body')
        self.VectorBodyMenu.Append(menu=self.VectorRBEdit, title='Edit Vector Body')
        G2G.Define_wxId('wxID_VECTORBODYIMP')
        self.VectorRBEdit.Append(G2G.wxID_VECTORBODYIMP,'Extract from file',
                                'Extract rigid body from phase file')
        G2G.Define_wxId('wxID_VECTORBODYSAV')
        self.VectorRBEdit.Append(G2G.wxID_VECTORBODYSAV,'Save rigid body','Write a rigid body to a file')
        G2G.Define_wxId('wxID_VECTORBODYRD')
        self.VectorRBEdit.Append(G2G.wxID_VECTORBODYRD,'Read rigid body','Read a rigid body from a file')
        G2G.Define_wxId('wxID_VECTORBODYEXTD')
        self.VectorRBEdit.Append(G2G.wxID_VECTORBODYEXTD,'Add translation',
                                'Add translation to existing rigid body')
        self.PostfillDataMenu()

        # Restraints
        G2G.Define_wxId('wxID_RESTRAINTADD', 'wxID_RESTDELETE', 'wxID_RESRCHANGEVAL',
            'wxID_RESTCHANGEESD', 'wxID_AARESTRAINTADD', 'wxID_AARESTRAINTPLOT','wxID_USEMOGUL')
        self.RestraintTab = wx.Menu(title='')
        self.RestraintEdit = wx.Menu(title='')
        self.RestraintEdit.Append(G2G.wxID_RESTRAINTADD,'Add restraints','Add restraints')
        self.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)    #gets disabled if macromolecule phase
        self.RestraintEdit.Append(G2G.wxID_AARESTRAINTADD,'Add residue restraints',
            'Add residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(G2G.wxID_AARESTRAINTADD,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(G2G.wxID_USEMOGUL,'Add MOGUL restraints',
            'Add restraints from MOGUL csv file')
        self.RestraintEdit.Enable(G2G.wxID_USEMOGUL,False)    #gets enabled if bonds or angles
        self.RestraintEdit.Append(G2G.wxID_AARESTRAINTPLOT,'Plot residue restraints',
            'Plot selected residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(G2G.wxID_AARESTRAINTPLOT,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(G2G.wxID_RESRCHANGEVAL,'Change value','Change observed value')
        self.RestraintEdit.Append(G2G.wxID_RESTCHANGEESD,'Change esd','Change esd in observed value')
        self.RestraintEdit.Append(G2G.wxID_RESTDELETE,'Delete restraints','Delete selected restraints')

        self.RestraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RestraintMenu)
        self.RestraintMenu.Append(menu=self.RestraintTab, title='Select tab')
        self.RestraintMenu.Append(menu=self.RestraintEdit, title='Edit Restr.')
        self.PostfillDataMenu()
            
        # Sequential results
        G2G.Define_wxId('wxID_RENAMESEQSEL', 'wxID_SAVESEQSEL', 'wxID_SAVESEQCSV', 'wxID_SAVESEQSELCSV', 'wxID_PLOTSEQSEL',
          'wxID_ORGSEQSEL', 'wxID_ADDSEQVAR', 'wxID_DELSEQVAR', 'wxID_EDITSEQVAR', 'wxID_COPYPARFIT', 'wxID_AVESEQSEL','wxID_SELECTUSE',
          'wxID_ADDPARFIT', 'wxID_DELPARFIT', 'wxID_EDITPARFIT', 'wxID_DOPARFIT', 'wxID_ADDSEQDIST', 'wxID_ADDSEQANGLE', 'wxID_ORGSEQINC',)
        self.SequentialMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SequentialMenu)
        self.SequentialFile = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialFile, title='Columns/Rows')
        self.SequentialFile.Append(G2G.wxID_SELECTUSE,'Set used...',
            'Dialog to select rows for plots/equation fitting')
        G2G.Define_wxId('wxID_UPDATESEQSEL')
        self.SequentialFile.Append(G2G.wxID_UPDATESEQSEL,'Update phase from selected',
            'Update phase information from selected row')
        self.SequentialFile.AppendSeparator()
        self.SequentialFile.Append(G2G.wxID_PLOTSEQSEL,'Plot selected cols',
            'Plot selected sequential refinement columns')
        self.SequentialFile.Append(G2G.wxID_RENAMESEQSEL,'Rename selected cols',
            'Rename selected sequential refinement columns')
        self.SequentialFile.Append(G2G.wxID_SAVESEQSEL,'Save selected as text',
            'Save selected sequential refinement results as a text file')
        self.SequentialFile.Append(G2G.wxID_SAVESEQSELCSV,'Save selected as CSV',
            'Save selected sequential refinement columns as a CSV spreadsheet file')
        self.SequentialFile.Append(G2G.wxID_AVESEQSEL,'Compute average',
            'Compute average for selected parameter')            
        self.SequentialFile.Append(G2G.wxID_ORGSEQINC,'Hide columns...',
            'Select columns to remove from displayed table')
        self.SequentialFile.AppendSeparator()       
        self.SequentialFile.Append(G2G.wxID_SAVESEQCSV,'Save all as CSV',
            'Save all sequential refinement results as a CSV spreadsheet file')
#        self.SequentialFile.Append(id=G2G.wxID_ORGSEQSEL, kind=wx.ITEM_NORMAL,text='Reorganize',
#            help='Reorganize variables where variables change')
        self.SequentialPvars = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialPvars, title='Pseudo Vars')
        self.SequentialPvars.Append(G2G.wxID_ADDSEQVAR,'Add Formula','Add a new custom pseudo-variable')
        self.SequentialPvars.Append(G2G.wxID_ADDSEQDIST,'Add Distance','Add a new bond distance pseudo-variable')
        self.SequentialPvars.Append(G2G.wxID_ADDSEQANGLE,'Add Angle','Add a new bond angle pseudo-variable')
        self.SequentialPvars.Append(G2G.wxID_DELSEQVAR,'Delete','Delete an existing pseudo-variable')
        self.SequentialPvars.Append(G2G.wxID_EDITSEQVAR,'Edit','Edit an existing pseudo-variable')

        self.SequentialPfit = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialPfit, title='Parametric Fit')
        self.SequentialPfit.Append(G2G.wxID_ADDPARFIT,'Add equation','Add a new equation to minimize')
        self.SequentialPfit.Append(G2G.wxID_COPYPARFIT,'Copy equation','Copy an equation to minimize - edit it next')
        self.SequentialPfit.Append(G2G.wxID_DELPARFIT,'Delete equation','Delete an equation for parametric minimization')
        self.SequentialPfit.Append(G2G.wxID_EDITPARFIT,'Edit equation','Edit an existing parametric minimization equation')
        self.SequentialPfit.Append(G2G.wxID_DOPARFIT,'Fit to equation(s)','Perform a parametric minimization')
        # fill sequential Export menu
        # for an exporter to be used for sequential exports, it must have a Writer method and
        # that Writer method must offer a mode argument.
        #============================================================
        # N.B. this largely duplicates menu items now in Export
        #============================================================
        self.SeqExportLookup = {}
        self.SequentialEx = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialEx, title='Seq Export')
        for lbl,txt in (('Phase','Export selected phase(s)'),
                        ('Project','Export entire sequential fit'),
                        ('Powder','Export selected powder histogram(s)'),
                        ('sasd','Export selected small angle histogram(s)')
                        ):
            objlist = []
            for obj in self.parent.GetTopLevelParent().exporterlist:
                if lbl.lower() in obj.exporttype:
                    try:
                        obj.Writer
                    except AttributeError:
                        continue
                    if '2' in platform.python_version_tuple()[0]:
                        if 'mode' in inspect.getargspec(obj.Writer)[0]:
                            objlist.append(obj)
                    else:
                        if 'mode' in inspect.getfullargspec(obj.Writer)[0]:
                            objlist.append(obj)
            if objlist:
                submenu = wx.Menu()
                item = self.SequentialEx.AppendSubMenu(submenu,lbl+' as',txt)
                for obj in objlist:
                    item = submenu.Append(wx.ID_ANY,obj.formatName,obj.longFormatName)
                    self.SeqExportLookup[item.GetId()] = (obj,lbl) # lookup table for submenu item
                    # Bind is in UpdateSeqResults
        
        self.PostfillDataMenu()
            
        # PWDR & SASD
        G2G.Define_wxId('wxID_PWDANALYSIS','wxID_PWDCOPY','wxID_PLOTCTRLCOPY','wxID_MERGEHKL',
            'wxID_PWDHKLPLOT', 'wxID_PWD3DHKLPLOT','wxID_3DALLHKLPLOT','wxID_1DHKLSTICKPLOT',)            
        self.PWDRMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PWDRMenu)
        self.ErrorAnal = wx.Menu(title='')
        self.PWDRMenu.Append(menu=self.ErrorAnal,title='Commands')
        self.ErrorAnal.Append(G2G.wxID_PWDANALYSIS,'Error Analysis','Error analysis on powder pattern')
        self.ErrorAnal.Append(G2G.wxID_PWDCOPY,'Copy params','Copy of PWDR parameters')
        self.ErrorAnal.Append(G2G.wxID_PLOTCTRLCOPY,'Copy plot controls','Copy of PWDR plot controls')
        self.moveDiffCurve = self.ErrorAnal.Append(wx.ID_ANY,'Move diff. curve',
            'Click on position where difference curve is placed')
        self.moveTickLoc = self.ErrorAnal.Append(wx.ID_ANY,'Move ticks','Move mouse to where tick marks should be positioned')
        self.moveTickSpc = self.ErrorAnal.Append(wx.ID_ANY,'Set tick space','Click to set spacing between phase tick marks')
        self.setPlotLim = self.ErrorAnal.Append(wx.ID_ANY,'Set plot limits...','Allows entry of plot min & max values')
        self.setPlotFmt = self.ErrorAnal.Append(wx.ID_ANY,'Set plot formatting...','Allows changes to text size and line widths, etc.')
        self.PostfillDataMenu()
            
        # HKLF - wxIDs defined in PWDR & SASD above
        self.HKLFMenu = wx.MenuBar()
        self.PrefillDataMenu(self.HKLFMenu)
        self.ErrorAnal = wx.Menu(title='')
        self.HKLFMenu.Append(menu=self.ErrorAnal,title='Commands')
        self.ErrorAnal.Append(G2G.wxID_PWDANALYSIS,'Error Analysis','Error analysis on single crystal data')
        self.ErrorAnal.Append(G2G.wxID_MERGEHKL,'Merge HKLs','Transform & merge HKLF data to new histogram')
        self.ErrorAnal.Append(G2G.wxID_1DHKLSTICKPLOT,'Plot 1D HKLs','Plot of HKLs from single crystal data in 1D')
        self.ErrorAnal.Append(G2G.wxID_PWD3DHKLPLOT,'Plot 3D HKLs','Plot HKLs from single crystal data in 3D')
        self.ErrorAnal.Append(G2G.wxID_3DALLHKLPLOT,'Plot all 3D HKLs','Plot HKLs from all single crystal data in 3D')
        self.ErrorAnal.Append(G2G.wxID_PWDCOPY,'Copy params','Copy of HKLF parameters')
        self.PostfillDataMenu()
            
        # PWDR / Limits
        G2G.Define_wxId('wxID_LIMITCOPY', 'wxID_ADDEXCLREGION',)
        self.LimitMenu = wx.MenuBar()
        self.PrefillDataMenu(self.LimitMenu)
        self.LimitEdit = wx.Menu(title='')
        self.LimitMenu.Append(menu=self.LimitEdit, title='Edit Limits')
        self.LimitEdit.Append(G2G.wxID_LIMITCOPY,'Copy','Copy limits to other histograms')
        self.LimitEdit.Append(G2G.wxID_ADDEXCLREGION,'Add exclude',
            'Add excluded region - select a point on plot; drag to adjust')            
        self.PostfillDataMenu()
            
        # PDR / Background
        G2G.Define_wxId('wxID_BACKCOPY', 'wxID_BACKFLAGCOPY','wxID_MAKEBACKRDF', 
            'wxID_RESCALEALL','wxID_BACKPEAKSMOVE','wxID_BACKSAVE','wxID_BACKLOAD')
        self.BackMenu = wx.MenuBar()
        self.PrefillDataMenu(self.BackMenu)
        self.BackEdit = wx.Menu(title='')
        self.BackMenu.Append(menu=self.BackEdit, title='File')
        self.BackEdit.Append(G2G.wxID_BACKCOPY,'Copy','Copy background parameters to other histograms')
        self.BackEdit.Append(G2G.wxID_BACKFLAGCOPY,'Copy flags',
            'Copy background refinement flags to other histograms')
        self.BackEdit.Append(G2G.wxID_BACKSAVE,'Save ...','Save background parameters to file')
        self.BackEdit.Append(G2G.wxID_BACKLOAD,'Load ...','Load background parameters from file')
        self.BackEdit.Append(G2G.wxID_BACKPEAKSMOVE,'Move peaks','Move background peaks to Peak List')
        self.BackEdit.Append(G2G.wxID_MAKEBACKRDF,'Plot RDF','Plot radial distribution from differences')
        self.BackFixed = wx.Menu(title='') # fixed background point menu
        self.BackMenu.Append(menu=self.BackFixed, title='Fixed Points')
        self.wxID_BackPts = {}
        self.wxID_BackPts['Add'] = wx.NewId() # N.B. not using wxID_ global as for other menu items
        self.BackFixed.AppendRadioItem(self.wxID_BackPts['Add'],'Add','Add fixed background points with mouse clicks')
        self.wxID_BackPts['Move'] = wx.NewId() 
        item = self.BackFixed.AppendRadioItem(self.wxID_BackPts['Move'],'Move','Move selected fixed background points with mouse drags')
        item.Check(True)
        self.wxID_BackPts['Del'] = wx.NewId()
        self.BackFixed.AppendRadioItem(self.wxID_BackPts['Del'],'Delete','Delete fixed background points with mouse clicks')
        self.wxID_BackPts['Clear'] = wx.NewId() 
        self.BackFixed.Append(self.wxID_BackPts['Clear'],'Clear','Clear fixed background points')
        self.wxID_BackPts['Fit'] = wx.NewId() 
        self.BackFixed.Append(self.wxID_BackPts['Fit'],'Fit background',
            'Fit background function to fixed background points')
        self.PostfillDataMenu()
            
        # PDR / Instrument Parameters
        G2G.Define_wxId('wxID_INSTPRMRESET','wxID_INSTCOPY','wxID_INSTFLAGCOPY','wxID_INSTLOAD',
            'wxID_INSTSAVE', 'wxID_INST1VAL', 'wxID_INSTCALIB', 'wxID_INSTSAVEALL',)
        self.InstMenu = wx.MenuBar()
        self.PrefillDataMenu(self.InstMenu)
        self.InstEdit = wx.Menu(title='')
        self.InstMenu.Append(menu=self.InstEdit, title='Operations')
        self.InstEdit.Append(G2G.wxID_INSTCALIB,'Calibrate','Calibrate from indexed peaks')
        self.InstEdit.Append(G2G.wxID_INSTPRMRESET,'Reset profile','Reset instrument profile parameters to default')
        self.InstEdit.Append(G2G.wxID_INSTLOAD,'Load profile...','Load instrument profile parameters from file')
        self.InstEdit.Append(G2G.wxID_INSTSAVE,'Save profile...','Save instrument profile parameters to file')
        self.InstEdit.Append(G2G.wxID_INSTSAVEALL,'Save all profile...','Save all instrument profile parameters to one file')
        self.InstEdit.Append(G2G.wxID_INSTCOPY,'Copy','Copy instrument profile parameters to other histograms')
        self.InstEdit.Append(G2G.wxID_INSTFLAGCOPY,'Copy flags','Copy instrument parameter refinement flags to other histograms')
        self.InstEdit.Append(G2G.wxID_INST1VAL,'Set one value','Set one instrument parameter value across multiple histograms')
        self.PostfillDataMenu()
        
        # PDR / Sample Parameters
        G2G.Define_wxId('wxID_SAMPLECOPY', 'wxID_SAMPLECOPYSOME', 'wxID_SAMPLEFLAGCOPY','wxID_SAMPLESAVE',
             'wxID_SAMPLELOAD', 'wxID_SETSCALE', 'wxID_SAMPLE1VAL', 'wxID_ALLSAMPLELOAD',)            
        self.SampleMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SampleMenu)
        self.SampleEdit = wx.Menu(title='')
        self.SampleMenu.Append(menu=self.SampleEdit, title='Command')
        self.SetScale = self.SampleEdit.Append(G2G.wxID_SETSCALE,'Set scale','Set scale by matching to another histogram')
        self.SampleEdit.Append(G2G.wxID_SAMPLELOAD,'Load','Load sample parameters from file')
        self.SampleEdit.Append(G2G.wxID_SAMPLESAVE,'Save','Save sample parameters to file')
        self.SampleEdit.Append(G2G.wxID_SAMPLECOPY,'Copy','Copy refinable and most other sample parameters to other histograms')
        self.SampleEdit.Append(G2G.wxID_SAMPLECOPYSOME,'Copy selected...','Copy selected sample parameters to other histograms')
        self.SampleEdit.Append(G2G.wxID_SAMPLEFLAGCOPY,'Copy flags','Copy sample parameter refinement flags to other histograms')
        self.SampleEdit.Append(G2G.wxID_SAMPLE1VAL,'Set one value','Set one sample parameter value across multiple histograms')
        self.SampleEdit.Append(G2G.wxID_ALLSAMPLELOAD,'Load all','Load sample parameters over multiple histograms')
        self.SampleEdit.Append(G2G.wxID_RESCALEALL,'Rescale all','Rescale all data with selected range')
        self.PostfillDataMenu()
        self.SetScale.Enable(False)

        # PDR / Peak List
        G2G.Define_wxId('wxID_UNDO', 'wxID_LSQPEAKFIT', 'wxID_LSQONECYCLE', 'wxID_RESETSIGGAM', 
            'wxID_CLEARPEAKS', 'wxID_AUTOSEARCH','wxID_PEAKSCOPY', 'wxID_SEQPEAKFIT','wxID_PEAKLOAD','wxID_PEAKSAVE')
        self.PeakMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PeakMenu)
        self.PeakEdit = wx.Menu(title='')
        self.PeakMenu.Append(menu=self.PeakEdit, title='Peak Fitting')
        self.peaksSel = self.PeakEdit.Append(wx.ID_ANY,'Set sel. ref flags...','Set refinement flags for selected peaks')
        self.peaksAll = self.PeakEdit.Append(wx.ID_ANY,'Set all ref flags...','Set refinement flags for all peaks')
        self.AutoSearch = self.PeakEdit.Append(G2G.wxID_AUTOSEARCH,'Auto search','Automatic peak search')
        self.UnDo = self.PeakEdit.Append(G2G.wxID_UNDO,'UnDo','Undo last least squares refinement')
        self.PeakFit = self.PeakEdit.Append(G2G.wxID_LSQPEAKFIT,'Peakfit','Peak fitting' )
        self.PFOneCycle = self.PeakEdit.Append(G2G.wxID_LSQONECYCLE,'Peakfit one cycle','One cycle of Peak fitting' )
        self.PeakEdit.Append(G2G.wxID_RESETSIGGAM,'Reset sig and gam','Reset sigma and gamma to global fit' )
        self.PeakCopy = self.PeakEdit.Append(G2G.wxID_PEAKSCOPY,'Peak copy','Copy peaks to other histograms')
        self.PeakEdit.Append(G2G.wxID_PEAKLOAD,'Load peaks...','Load peak list from file')
        self.PeakEdit.Append(G2G.wxID_PEAKSAVE,'Save peaks...','Save peak list to file')
        self.SeqPeakFit = self.PeakEdit.Append(G2G.wxID_SEQPEAKFIT,'Seq PeakFit', 
            'Sequential Peak fitting for all histograms' )
        self.PeakEdit.Append(G2G.wxID_CLEARPEAKS,'Clear peaks','Clear the peak list' )
        self.movePeak = self.PeakEdit.Append(wx.ID_ANY,'Move selected peak',
            'Select a peak in the table, then use this to move it with the mouse.')
        self.PostfillDataMenu()
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.PFOneCycle.Enable(False)
        self.AutoSearch.Enable(True)
        
        # PDR / Index Peak List
        G2G.Define_wxId('wxID_INDXRELOAD','wxID_INDEXSAVE',)
        self.IndPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndPeaksMenu)
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndPeaksMenu.Append(menu=self.IndPeaksEdit,title='Operations')
        self.IndPeaksEdit.Append(G2G.wxID_INDXRELOAD,'Load/Reload','Load/Reload index peaks from peak list')
        self.IndPeaksEdit.Append(G2G.wxID_INDEXSAVE,'Save','Save index peaks to CSV file')
        self.PostfillDataMenu()
        
        # PDR / Unit Cells List
        G2G.Define_wxId('wxID_INDEXPEAKS', 'wxID_REFINECELL', 'wxID_COPYCELL', 'wxID_MAKENEWPHASE',
            'wxID_EXPORTCELLS','wxID_LOADCELL','wxID_IMPORTCELL','wxID_TRANSFORMCELL','wxID_RUNSUB','wxID_RUNSUBMAG','wxID_LATSYM')
        self.IndexMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndexMenu)
        self.IndexEdit = wx.Menu(title='')
        self.IndexMenu.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        self.IndexPeaks = self.IndexEdit.Append(G2G.wxID_INDEXPEAKS,'Index Cell',
            'Find cells that index fitted peaks')
        self.IndexEdit.Append(G2G.wxID_LATSYM,'Cell Symmetry Search',
            'Run Bilboa "Lattice Symmetry" to find higher symmetry cells')
        self.RunSubGroups = self.IndexEdit.Append(G2G.wxID_RUNSUB,'Run SUBGROUPS',
            'If disabled, do Load Phase first')
        self.RunSubGroupsMag = self.IndexEdit.Append(G2G.wxID_RUNSUBMAG,'Run k-SUBGROUPMAG',
            'If disabled, do Load Phase first')
        self.CopyCell = self.IndexEdit.Append(G2G.wxID_COPYCELL,'Copy Cell', 
            'Copy selected unit cell from indexing to cell refinement fields')
        self.LoadCell = self.IndexEdit.Append(G2G.wxID_LOADCELL,'Load Phase', 
            'Load unit cell from a phase tree entry')
        self.ImportCell = self.IndexEdit.Append(G2G.wxID_IMPORTCELL,'Import Cell', 
            'Import unit cell from file')
        self.TransposeCell = self.IndexEdit.Append(G2G.wxID_TRANSFORMCELL,'Transform Cell', 
            'Transform unit cell')
        self.RefineCell = self.IndexEdit.Append(G2G.wxID_REFINECELL,'Refine Cell',
            'Refine unit cell parameters from indexed peaks')
        self.MakeNewPhase = self.IndexEdit.Append(G2G.wxID_MAKENEWPHASE,'Make new phase',
            'Make new phase from selected unit cell')
        self.ExportCells = self.IndexEdit.Append(G2G.wxID_EXPORTCELLS,'Export cell list','Export cell list to csv file')
        self.PostfillDataMenu()
        self.LoadCell.Enable(False)
        self.IndexPeaks.Enable(False)
        self.RunSubGroups.Enable(False)
        self.RunSubGroupsMag.Enable(False)
        self.CopyCell.Enable(False)
        self.RefineCell.Enable(False)
        self.MakeNewPhase.Enable(False)
        
        # PDR / Reflection Lists
        G2G.Define_wxId('wxID_SELECTPHASE', ) #some wxIDs defined above in PWDR & SASD
        self.ReflMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ReflMenu)
        self.ReflEdit = wx.Menu(title='')
        self.ReflMenu.Append(menu=self.ReflEdit, title='Reflection List')
        self.SelectPhase = self.ReflEdit.Append(G2G.wxID_SELECTPHASE,'Select phase','Select phase for reflection list')
        self.ReflEdit.Append(G2G.wxID_1DHKLSTICKPLOT,'Plot 1D HKLs','Plot of HKLs in 1D')
        self.ReflEdit.Append(G2G.wxID_PWDHKLPLOT,'Plot HKLs','Plot HKLs in 2D')
        self.ReflEdit.Append(G2G.wxID_PWD3DHKLPLOT,'Plot 3D HKLs','Plot HKLs in 3D')
        self.PostfillDataMenu()
        
        # SASD / Instrument Parameters
        G2G.Define_wxId('wxID_SASDINSTCOPY',)
        self.SASDInstMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SASDInstMenu)
        self.SASDInstEdit = wx.Menu(title='')
        self.SASDInstMenu.Append(menu=self.SASDInstEdit, title='Operations')
        self.SASDInstEdit.Append(G2G.wxID_SASDINSTCOPY,'Copy','Copy instrument profile parameters to other histograms')
        self.PostfillDataMenu()
        
        #SASD & REFL/ Substance editor
        G2G.Define_wxId('wxID_LOADSUBSTANCE','wxID_RELOADSUBSTANCES','wxID_ADDSUBSTANCE','wxID_COPYSUBSTANCE',
            'wxID_DELETESUBSTANCE','wxID_ELEMENTADD', 'wxID_ELEMENTDELETE',)    
        self.SubstanceMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SubstanceMenu)
        self.SubstanceEdit = wx.Menu(title='')
        self.SubstanceMenu.Append(menu=self.SubstanceEdit, title='Edit substance')
        self.SubstanceEdit.Append(G2G.wxID_LOADSUBSTANCE,'Load substance','Load substance from file')
        self.SubstanceEdit.Append(G2G.wxID_RELOADSUBSTANCES,'Reload substances','Reload all substances from file')
        self.SubstanceEdit.Append(G2G.wxID_ADDSUBSTANCE,'Add substance','Add new substance to list')
        self.SubstanceEdit.Append(G2G.wxID_COPYSUBSTANCE,'Copy substances','Copy substances')
        self.SubstanceEdit.Append(G2G.wxID_DELETESUBSTANCE,'Delete substance','Delete substance from list')            
        self.SubstanceEdit.Append(G2G.wxID_ELEMENTADD,'Add elements','Add elements to substance')
        self.SubstanceEdit.Append(G2G.wxID_ELEMENTDELETE,'Delete elements','Delete elements from substance')
        self.PostfillDataMenu()
        
        # SASD/ Models
        G2G.Define_wxId('wxID_MODELCOPY', 'wxID_MODELFIT', 'wxID_MODELADD','wxID_MODELUNDO', 
            'wxID_MODELFITALL', 'wxID_MODELCOPYFLAGS','wxID_MODELPLOT',)    
        self.ModelMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ModelMenu)
        self.ModelEdit = wx.Menu(title='')
        self.ModelMenu.Append(menu=self.ModelEdit, title='Models')
        self.ModelEdit.Append(G2G.wxID_MODELADD,'Add','Add new term to model')
        self.ModelEdit.Append(G2G.wxID_MODELFIT,'Fit','Fit model parameters to data')
        self.SasdUndo = self.ModelEdit.Append(G2G.wxID_MODELUNDO,'Undo','Undo model fit')
        self.SasdUndo.Enable(False)            
        self.SasSeqFit = self.ModelEdit.Append(G2G.wxID_MODELFITALL,'Sequential fit','Sequential fit of model parameters to all SASD data')
        self.SasSeqFit.Enable(False)
        self.ModelEdit.Append(G2G.wxID_MODELCOPY,'Copy','Copy model parameters to other histograms')
        self.ModelEdit.Append(G2G.wxID_MODELCOPYFLAGS,'Copy flags','Copy model refinement flags to other histograms')
        self.PostfillDataMenu()
        
        # REFD/ Models - wxIDs as for SASD/Models
        self.REFDModelMenu = wx.MenuBar()
        self.PrefillDataMenu(self.REFDModelMenu)
        self.REFDModelEdit = wx.Menu(title='')
        self.REFDModelMenu.Append(menu=self.REFDModelEdit, title='Models')
        self.REFDModelEdit.Append(G2G.wxID_MODELFIT,'Fit','Fit model parameters to data')
        self.REFDUndo = self.REFDModelEdit.Append(G2G.wxID_MODELUNDO,'Undo','Undo model fit')
        self.REFDUndo.Enable(False)            
        self.REFDModelEdit.Append(G2G.wxID_MODELFITALL,'Sequential fit','Sequential fit of model parameters to all REFD data')
        self.REFDModelEdit.Append(G2G.wxID_MODELCOPY,'Copy','Copy model parameters to other histograms')
        self.REFDModelEdit.Append(G2G.wxID_MODELPLOT,'Plot','Plot model SDL for selected histograms')
        self.PostfillDataMenu()

        # IMG / Image Controls
        G2G.Define_wxId('wxID_IMCALIBRATE', 'wxID_IMRECALIBRATE', 'wxID_IMINTEGRATE', 'wxID_IMCLEARCALIB', 'wxID_IMRECALIBALL', 
            'wxID_IMCOPYCONTROLS', 'wxID_INTEGRATEALL', 'wxID_IMSAVECONTROLS', 'wxID_IMLOADCONTROLS', 'wxID_IMAUTOINTEG',
            'wxID_IMCOPYSELECTED', 'wxID_SAVESELECTEDCONTROLS', 'wxID_IMXFERCONTROLS', 'wxID_IMRESETDIST',
            'wxID_LOADELECTEDCONTROLS')
        self.ImageMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ImageMenu)
        self.ImageEdit = wx.Menu(title='')
        self.ImageMenu.Append(menu=self.ImageEdit, title='Calibration')
        self.ImageEdit.Append(G2G.wxID_IMCALIBRATE,'Calibrate','Calibrate detector by fitting to calibrant lines')
        self.ImageEdit.Append(G2G.wxID_IMRECALIBRATE,'Recalibrate','Recalibrate detector by fitting to calibrant lines')
        self.ImageEdit.Append(G2G.wxID_IMRECALIBALL,'Recalibrate all','Recalibrate all images by fitting to calibrant lines')
        G2G.Define_wxId('wxID_IMDISTRECALIB')
        self.ImageEdit.Append(G2G.wxID_IMDISTRECALIB,'Multi-distance Recalibrate','Recalibrate all images varying delta-distance and fitting wavelength')
        self.ImageEdit.Append(G2G.wxID_IMCLEARCALIB,'Clear calibration','Clear calibration data points and rings')
        
        ImageIntegrate = wx.Menu(title='')
        self.ImageMenu.Append(menu=ImageIntegrate, title='Integration')
        ImageIntegrate.Append(G2G.wxID_IMINTEGRATE,'Integrate','Integrate selected image')
        ImageIntegrate.Append(G2G.wxID_INTEGRATEALL,'Integrate all','Integrate all images selected from list')
        ImageIntegrate.Append(G2G.wxID_IMAUTOINTEG,'Auto Integrate','Open Auto-integration window to integrate a series of images')
        G2G.Define_wxId('wxID_IMINTEGPDFTOOL')
        ImageIntegrate.Append(G2G.wxID_IMINTEGPDFTOOL,'Integrate/PDF app (in dev)','Start Integration/PDF task (in development)')

        ImageParams = wx.Menu(title='')
        self.ImageMenu.Append(menu=ImageParams, title='Parms')
        ImageParams.Append(G2G.wxID_IMCOPYCONTROLS,'Copy Controls','Copy image controls to other images')
        ImageParams.Append(G2G.wxID_IMCOPYSELECTED,'Copy Selected','Copy selected image controls to other images')
        ImageParams.Append(G2G.wxID_IMSAVECONTROLS,'Save Controls','Save image controls to file')
        ImageParams.Append(G2G.wxID_SAVESELECTEDCONTROLS,'Save Multiple Controls','Save controls from selected images to file')
        ImageParams.Append(G2G.wxID_IMLOADCONTROLS,'Load Controls','Load image controls from file')
        ImageParams.Append(G2G.wxID_LOADELECTEDCONTROLS,'Load Multiple Controls','Load multiple image controls from multiple files')
        ImageParams.Append(G2G.wxID_IMXFERCONTROLS,'Xfer controls','Transfer integration controls to other detector distances')
        ImageParams.Append(G2G.wxID_IMRESETDIST,'Reset dist','Reset all detector dist to set dist')
        self.PostfillDataMenu()
            
        # IMG / Masks
        G2G.Define_wxId('wxID_MASKCOPY', 'wxID_MASKSAVE', 'wxID_MASKLOAD', 'wxID_NEWMASKSPOT', 'wxID_NEWMASKARC', 'wxID_NEWMASKRING',
            'wxID_NEWMASKFRAME', 'wxID_NEWMASKPOLY','wxID_NEWMASKXLINE','wxID_NEWMASKYLINE','wxID_MASKLOADNOT',
            'wxID_FINDSPOTS', 'wxID_AUTOFINDSPOTS', 'wxID_DELETESPOTS',)
        self.MaskMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MaskMenu)
        self.MaskEdit = wx.Menu(title='')
        self.MaskMenu.Append(menu=self.MaskEdit, title='Operations')
        submenu = wx.Menu()
        self.MaskEdit.AppendSubMenu( submenu,'Create new','')
        self.MaskEdit.Append(G2G.wxID_MASKCOPY,'Copy mask','Copy mask to other images')
        self.MaskEdit.Append(G2G.wxID_MASKSAVE,'Save mask','Save mask to file')
        self.MaskEdit.Append(G2G.wxID_MASKLOADNOT,'Load mask','Load mask from file; ignoring threshold')
        self.MaskEdit.Append(G2G.wxID_MASKLOAD,'Load mask w/threshold','Load mask from file keeping the threshold value')
        self.MaskEdit.Append(G2G.wxID_FINDSPOTS,'Spot mask search','Search for spot mask; NB: slow')
        self.MaskEdit.Append(G2G.wxID_AUTOFINDSPOTS,'Auto spot mask search','Auto spot mask ssearch; NB: slow')
        self.MaskEdit.Append(G2G.wxID_DELETESPOTS,'Delete spot masks','Delete all spot masks')
        submenu.Append(G2G.wxID_NEWMASKARC,'Arc mask','Create an arc mask with mouse input')
        submenu.Append(G2G.wxID_NEWMASKFRAME,'Frame mask','Create a frame mask with mouse input')
        submenu.Append(G2G.wxID_NEWMASKPOLY,'Polygon mask','Create a polygon mask with mouse input')
        submenu.Append(G2G.wxID_NEWMASKRING,'Ring mask','Create a ring mask with mouse input')
        submenu.Append(G2G.wxID_NEWMASKSPOT,'Spot mask','Create spot masks with mouse input')
        submenu.Append(G2G.wxID_NEWMASKXLINE,'X line mask','Create line masks with mouse input')
        submenu.Append(G2G.wxID_NEWMASKYLINE,'Y line mask','Create line masks with mouse input')
        self.PostfillDataMenu()
            
        # IMG / Stress/Strain
        G2G.Define_wxId('wxID_STRSTACOPY', 'wxID_STRSTAFIT', 'wxID_STRSTASAVE', 'wxID_STRSTALOAD', 'wxID_STRSTSAMPLE',
            'wxID_APPENDDZERO', 'wxID_STRSTAALLFIT', 'wxID_UPDATEDZERO', 'wxID_STRSTAPLOT', 'wxID_STRRINGSAVE',)        
        self.StrStaMenu = wx.MenuBar()
        self.PrefillDataMenu(self.StrStaMenu)
        self.StrStaEdit = wx.Menu(title='')
        self.StrStaMenu.Append(menu=self.StrStaEdit, title='Operations')
        self.StrStaEdit.Append(G2G.wxID_APPENDDZERO,'Append d-zero','Append d-zero for one ring')
        self.StrStaEdit.Append(G2G.wxID_STRSTAFIT,'Fit stress/strain','Fit stress/strain data')
        self.StrStaEdit.Append(G2G.wxID_STRSTAPLOT,'Plot intensity distribution','Plot intensity distribution')
        self.StrStaEdit.Append(G2G.wxID_STRRINGSAVE,'Save intensity distribution','Save intensity distribution')
        self.StrStaEdit.Append(G2G.wxID_UPDATEDZERO,'Update d-zero','Update d-zero from ave d-zero')        
        self.StrStaEdit.Append(G2G.wxID_STRSTAALLFIT,'All image fit','Fit stress/strain data for all images')
        self.StrStaEdit.Append(G2G.wxID_STRSTACOPY,'Copy stress/strain','Copy stress/strain data to other images')
        self.StrStaEdit.Append(G2G.wxID_STRSTASAVE,'Save stress/strain','Save stress/strain data to file')
        self.StrStaEdit.Append(G2G.wxID_STRSTALOAD,'Load stress/strain','Load stress/strain data from file')
        self.StrStaEdit.Append(G2G.wxID_STRSTSAMPLE,'Load sample data','Load sample data from file')
        self.PostfillDataMenu()
            
        # PDF / PDF Controls
        G2G.Define_wxId('wxID_PDFCOPYCONTROLS', 'wxID_PDFSAVECONTROLS', 'wxID_PDFLOADCONTROLS', 'wxID_PDFCOMPUTE',
            'wxID_PDFCOMPUTEALL', 'wxID_PDFADDELEMENT', 'wxID_PDFDELELEMENT',)
        self.PDFMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PDFMenu)
        self.PDFEdit = wx.Menu(title='')
        self.PDFMenu.Append(menu=self.PDFEdit, title='PDF Controls')
        self.PDFEdit.Append(G2G.wxID_PDFADDELEMENT,'Add elements','Add one or more elements to sample composition')
        self.PDFEdit.Append(G2G.wxID_PDFDELELEMENT,'Delete element','Delete element from sample composition')
        self.PDFEdit.Append(G2G.wxID_PDFCOPYCONTROLS,'Copy controls','Copy PDF controls')
        self.PDFEdit.Append(G2G.wxID_PDFLOADCONTROLS,'Load Controls','Load PDF controls from file')
        self.PDFEdit.Append(G2G.wxID_PDFSAVECONTROLS,'Save controls','Save PDF controls to file')
        self.PDFEdit.Append(G2G.wxID_PDFCOMPUTE,'Compute PDF','Compute PDF')
        self.PDFEdit.Append(G2G.wxID_PDFCOMPUTEALL,'Compute all PDFs','Compute all PDFs with or w/o optimization')
        self.PostfillDataMenu()
        
        # PDF / PDF Peaks
        G2G.Define_wxId('wxID_PDFPKSFIT','wxID_PDFPKSFITALL', 'wxID_PDFCOPYPEAKS', 'wxID_CLEARPDFPEAKS',)
        self.PDFPksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PDFPksMenu)
        self.PDFPksEdit = wx.Menu(title='')
        self.PDFPksMenu.Append(menu=self.PDFPksEdit, title='PDF Peaks')
        self.PDFPksEdit.Append(G2G.wxID_PDFPKSFIT,'PDF peak fit','Fit PDF peaks')
        self.PDFPksEdit.Append(G2G.wxID_PDFPKSFITALL,'Seq PDF peak fit','Sequential Peak fitting for all PDFs')
        self.PDFPksEdit.Append(G2G.wxID_PDFCOPYPEAKS,'Copy peaks','Copy PDF peaks')
        self.PDFPksEdit.Append(G2G.wxID_CLEARPDFPEAKS,'Clear peaks','Clear PDF peaks')        
        self.PostfillDataMenu()

        # Phase / General tab
        G2G.Define_wxId('wxID_FOURCALC', 'wxID_FOURSEARCH', 'wxID_FOURCLEAR','wxID_CHARGEFLIP','wxID_VALIDPROTEIN', 
            'wxID_MULTIMCSA','wxID_SINGLEMCSA', 'wxID_4DCHARGEFLIP', 'wxID_TRANSFORMSTRUCTURE','wxID_USEBILBAOMAG',
            'wxID_COMPARESTRUCTURE',)
        self.DataGeneral = wx.MenuBar()
        self.PrefillDataMenu(self.DataGeneral)
        self.DataGeneral.Append(menu=wx.Menu(title=''),title='Select tab')
        self.GeneralCalc = wx.Menu(title='')
        self.DataGeneral.Append(menu=self.GeneralCalc,title='Compute')
        self.GeneralCalc.Append(G2G.wxID_FOURCALC,'Fourier map','Compute Fourier map',)
        self.GeneralCalc.Append(G2G.wxID_FOURSEARCH,'Search map','Search Fourier map')
        self.GeneralCalc.Append(G2G.wxID_CHARGEFLIP,'Charge flipping','Run charge flipping')
        self.GeneralCalc.Append(G2G.wxID_4DCHARGEFLIP,'4D Charge flipping','Run 4D charge flipping')
        self.GeneralCalc.Enable(G2G.wxID_4DCHARGEFLIP,False)   
        self.GeneralCalc.Append(G2G.wxID_FOURCLEAR,'Clear map','Clear map')
        self.GeneralCalc.Append(G2G.wxID_SINGLEMCSA,'MC/SA','Run Monte Carlo - Simulated Annealing')
        self.GeneralCalc.Append(G2G.wxID_MULTIMCSA,'Multi MC/SA','Run Monte Carlo - Simulated Annealing on multiprocessors')
        self.GeneralCalc.Append(G2G.wxID_TRANSFORMSTRUCTURE,'Transform','Transform crystal structure')
        self.GeneralCalc.Append(G2G.wxID_COMPARESTRUCTURE,'Compare','Compare polyhedra to ideal octahedra/tetrahedra')
        self.GeneralCalc.Enable(G2G.wxID_COMPARESTRUCTURE,False)   
        self.GeneralCalc.Append(G2G.wxID_USEBILBAOMAG,'Select magnetic/subgroup phase','If disabled, make in PWDR/Unit Cells')        
        self.GeneralCalc.Append(G2G.wxID_VALIDPROTEIN,'Protein quality','Protein quality analysis')
        self.PostfillDataMenu()
        
        # Phase / Data tab
        G2G.Define_wxId('wxID_DATACOPY', 'wxID_DATACOPYFLAGS', 'wxID_DATASELCOPY', 'wxID_DATAUSE',
            'wxID_PWDRADD', 'wxID_HKLFADD','wxID_DATADELETE', )
        self.DataMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataMenu)
        self.DataMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DataEdit = wx.Menu(title='')
        self.DataMenu.Append(menu=self.DataEdit, title='Edit Phase')
        self.DataEdit.Append(G2G.wxID_DATACOPY,'Copy data','Copy phase data to other histograms')
        self.DataEdit.Append(G2G.wxID_DATACOPYFLAGS,'Copy flags','Copy phase data flags to other histograms')
        self.DataEdit.Append(G2G.wxID_DATASELCOPY,'Copy selected data','Copy selected phase data to other histograms')
        self.DataEdit.Append(G2G.wxID_DATAUSE,'Select used data','Select all histograms to use')
        self.DataEdit.Append(G2G.wxID_PWDRADD,'Add powder histograms','Select new powder histograms to be used for this phase')
        self.DataEdit.Append(G2G.wxID_HKLFADD,'Add single crystal histograms','Select new single crystal histograms to be used for this phase')
        self.DataEdit.Append(G2G.wxID_DATADELETE,'Remove histograms','Remove histograms from use for this phase')
        G2G.Define_wxId('wxID_DATADIJ')
        self.DataEdit.Append(G2G.wxID_DATADIJ,'Apply Strain to Lattice Constants',
                             'Shift cell by Dij of selected histogram')
        self.PostfillDataMenu()
            
        # Phase / Atoms tab
        G2G.Define_wxId('wxID_ATOMSEDITADD', 'wxID_ATOMSEDITINSERT', 'wxID_ATOMSEDITDELETE',
            'wxID_ATOMSMODIFY', 'wxID_ATOMSTRANSFORM', 'wxID_ATOMSVIEWADD', 'wxID_ATOMVIEWINSERT',
            'wxID_RELOADDRAWATOMS', 'wxID_ATOMSDISAGL', 'wxID_ATOMMOVE', 'wxID_MAKEMOLECULE',
            'wxID_ATOMSPDISAGL', 'wxID_ISODISP', 'wxID_ADDHATOM', 'wxID_UPDATEHATOM',
            'wxID_ATOMSROTATE', 'wxID_ATOMSDENSITY','wxID_ATOMSBNDANGLHIST', 
            'wxID_ATOMSSETALL', 'wxID_ATOMSSETSEL',)
        self.AtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.AtomsMenu)
        self.AtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.AtomEdit = wx.Menu(title='')
        self.AtomCompute = wx.Menu(title='')
        self.AtomsMenu.Append(menu=self.AtomEdit, title='Edit Atoms')
        self.AtomsMenu.Append(menu=self.AtomCompute, title='Compute')
        submenu = wx.Menu()
        self.AtomEdit.AppendSubMenu(submenu,'On selected atoms...','Set/Act on selected atoms')
        submenu.Append(G2G.wxID_ATOMSSETSEL,'Refine selected','Set refinement flags for selected atoms')
        submenu.Append(G2G.wxID_ATOMSMODIFY,'Modify parameters','Modify parameters values for all selected atoms')
        submenu.Append(G2G.wxID_ATOMSEDITINSERT,'Insert atom','Inserts an H atom before all selected atoms')
        submenu.Append(G2G.wxID_ADDHATOM,'Calc H atoms','Insert H atoms in expected bonding positions for selected atoms')
        submenu.Append(G2G.wxID_ATOMSEDITDELETE,'Delete atom','Delete selected atoms')
        submenu.Append(G2G.wxID_ATOMSTRANSFORM,'Transform atoms','Symmetry transform selected atoms')
        submenu.Append(G2G.wxID_ATOMSSETALL,'Select All','Set refinement flags for all atoms')
        
        self.AtomEdit.Append(G2G.wxID_ATOMSEDITADD,'Append atom','Appended as an H atom')
        self.AtomEdit.Append(G2G.wxID_ATOMSVIEWADD,'Append view point','Appended as an H atom')
        self.AtomEdit.Append(G2G.wxID_ATOMVIEWINSERT,'Insert view point','Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(G2G.wxID_UPDATEHATOM,'Update H atoms','Update H atoms in standard positions')
        self.AtomEdit.Append(G2G.wxID_ATOMMOVE,'Move selected atom to view point','Select a single atom to be moved to view point in plot')
        self.AtomEdit.Append(G2G.wxID_MAKEMOLECULE,'Assemble molecule','Select a single atom to assemble as a molecule from scattered atom positions')
        self.AtomEdit.Append(G2G.wxID_RELOADDRAWATOMS,'Reload draw atoms','Reload atom drawing list')
        submenu = wx.Menu()
        self.AtomEdit.AppendSubMenu(submenu,'Reimport atoms','Reimport atoms from file; sequence must match')
        # setup a cascade menu for the formats that have been defined
        self.ReImportMenuId = {}  # points to readers for each menu entry
        for reader in self.parent.GetTopLevelParent().ImportPhaseReaderlist:
            item = submenu.Append(wx.ID_ANY,'reimport coordinates from '+reader.formatName+' file',reader.longFormatName)
            self.ReImportMenuId[item.GetId()] = reader
        item = submenu.Append(wx.ID_ANY,'guess format from file','Reimport coordinates, try to determine format from file')
        self.ReImportMenuId[item.GetId()] = None # try all readers
        self.AtomCompute.Append(G2G.wxID_ATOMSDISAGL,'Show Distances && Angles','Compute distances & angles for selected atoms')
        self.AtomCompute.Append(G2G.wxID_ATOMSPDISAGL,'Save Distances && Angles','Compute distances & angles for selected atoms')
        self.AtomCompute.Append(G2G.wxID_ATOMSBNDANGLHIST,'Histogram Bonds && Angles','Histogram bonds & angles for selected atoms')
        self.AtomCompute.Append(G2G.wxID_ATOMSDENSITY,'Density','Compute density for current phase')
        self.AtomCompute.ISOcalc = self.AtomCompute.Append(G2G.wxID_ISODISP,'ISODISTORT mode values',
            'Compute values of ISODISTORT modes from atom parameters')
        self.PostfillDataMenu()
        
        # Phase / Imcommensurate "waves" tab 
        G2G.Define_wxId('wxID_WAVEVARY',)
        self.WavesData = wx.MenuBar()
        self.PrefillDataMenu(self.WavesData)
        self.WavesData.Append(menu=wx.Menu(title=''),title='Select tab')
        self.WavesDataEdit = wx.Menu(title='')
        self.WavesData.Append(menu=self.WavesDataEdit, title='Edit Wave')
        self.WavesDataEdit.Append(G2G.wxID_WAVEVARY,'Global wave vary','Global setting of wave vary flags')
        self.PostfillDataMenu()
        
        #Phase / Dysnomia (Maximum Entropy Method) tab
        G2G.Define_wxId('wxID_LOADDYSNOMIA', 'wxID_SAVEDYSNOMIA', 'wxID_RUNDYSNOMIA', )       
        self.MEMMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MEMMenu)
        self.MEMMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MEMDataEdit = wx.Menu(title='')
        self.MEMMenu.Append(menu=self.MEMDataEdit, title='Operations')
        self.MEMDataEdit.Append(G2G.wxID_LOADDYSNOMIA,'Load from Dysnomia file','Load MEM info from Dysnomia file')
        self.MEMDataEdit.Append(G2G.wxID_SAVEDYSNOMIA,'Save Dysnomia file','Save MEM info in Dysnomia file')
        self.MEMDataEdit.Append(G2G.wxID_RUNDYSNOMIA,'Run Dysonmia','Run Dysnomia to make new Fobs map')
        self.PostfillDataMenu()
        
        #Phase / fullrmc & RMCprofile (Reverse Monte Carlo method) tab
        G2G.Define_wxId('wxID_SETUPRMC','wxID_RUNRMC','wxID_VIEWRMC',
                            'wxID_STOPRMC')       
        self.FRMCMenu = wx.MenuBar()
        self.PrefillDataMenu(self.FRMCMenu)
        self.FRMCMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.FRMCDataEdit = wx.Menu(title='')
        self.FRMCMenu.Append(menu=self.FRMCDataEdit, title='Operations')
        self.FRMCDataEdit.Append(G2G.wxID_SETUPRMC,'Setup RMC','Setup new fullrmc or RMCprofile file')
        self.FRMCDataEdit.Append(G2G.wxID_RUNRMC,'Execute','Run fullrmc or RMCprofile file')
        self.FRMCDataEdit.Append(G2G.wxID_STOPRMC,'Stop run','Stop fullrmc run')
        self.FRMCDataEdit.Append(G2G.wxID_VIEWRMC,'Plot','View fullrmc or RMCprofile results')
        self.PostfillDataMenu()
        
        # Phase / Layer tab 
        G2G.Define_wxId('wxID_LOADDIFFAX', 'wxID_LAYERSIMULATE', 'wxID_SEQUENCESIMULATE', 'wxID_LAYERSFIT', 'wxID_COPYPHASE',)       
        self.LayerData = wx.MenuBar()
        self.PrefillDataMenu(self.LayerData)
        self.LayerData.Append(menu=wx.Menu(title=''),title='Select tab')
        self.LayerDataEdit = wx.Menu(title='')
        self.LayerData.Append(menu=self.LayerDataEdit, title='Operations')
        self.LayerDataEdit.Append(G2G.wxID_LOADDIFFAX,'Load from DIFFaX file','Load layer info from DIFFaX file')
        self.LayerDataEdit.Append(G2G.wxID_COPYPHASE,'Copy phase cell','Copy phase cell from another project')
        self.LayerDataEdit.Append(G2G.wxID_LAYERSIMULATE,'Simulate pattern','Simulate diffraction pattern from layer stacking')
        self.LayerDataEdit.Append(G2G.wxID_LAYERSFIT,'Fit pattern','Fit diffraction pattern with layer stacking model')
        self.LayerDataEdit.Append(G2G.wxID_SEQUENCESIMULATE,'Sequence simulations','Sequence simulation changing one parameter')
        self.PostfillDataMenu()
                 
        # Phase / Draw Options tab
        self.DataDrawOptions = wx.MenuBar()
        self.PrefillDataMenu(self.DataDrawOptions)
        self.DataDrawOptions.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PostfillDataMenu()
        
        # Phase / Draw Atoms tab 
        G2G.Define_wxId('wxID_DRAWATOMSTYLE', 'wxID_DRAWATOMLABEL', 'wxID_DRAWATOMCOLOR', 'wxID_DRAWATOMRESETCOLOR',
            'wxID_DRAWVIEWPOINT', 'wxID_DRAWTRANSFORM', 'wxID_DRAWDELETE', 'wxID_DRAWFILLCELL',
            'wxID_DRAWADDEQUIV', 'wxID_DRAWFILLCOORD', 'wxID_DRAWDISAGLTOR', 'wxID_DRAWPLANE',
            'wxID_DRAWDISTVP', 'wxID_DRAWADDSPHERE', 'wxID_DRWAEDITRADII',)
        G2G.Define_wxId('wxID_DRAWRESTRBOND', 'wxID_DRAWRESTRANGLE', 'wxID_DRAWRESTRPLANE', 'wxID_DRAWRESTRCHIRAL',)        
        self.DrawAtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DrawAtomsMenu)
        self.DrawAtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DrawAtomEdit = wx.Menu(title='')
        self.DrawAtomCompute = wx.Menu(title='')
        self.DrawAtomRestraint = wx.Menu(title='')
        self.DrawAtomRigidBody = wx.Menu(title='')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomEdit, title='Edit Figure')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomCompute,title='Compute')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRestraint, title='Restraints')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRigidBody, title='Rigid body')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWATOMSTYLE,'Atom style','Select atoms first')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWATOMLABEL,'Atom label','Select atoms first')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWATOMCOLOR,'Atom color','Select atoms first')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWATOMRESETCOLOR,'Reset atom colors','Resets all atom colors to defaults')
        self.DrawAtomEdit.Append(G2G.wxID_DRWAEDITRADII,'Edit atom radii','Edit drawing atom radii') # TODO: removed until it can be made to do something
        self.DrawAtomEdit.Append(G2G.wxID_DRAWVIEWPOINT,'View point','View point is 1st atom selected')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWADDEQUIV,'Add atoms','Add symmetry & cell equivalents to drawing set from selected atoms')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWADDSPHERE,'Add sphere of atoms','Add atoms within sphere of enclosure')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWTRANSFORM,'Transform draw atoms','Transform selected atoms by symmetry & cell translations')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWFILLCOORD,'Fill CN-sphere','Fill coordination sphere for selected atoms')            
        self.DrawAtomEdit.Append(G2G.wxID_DRAWFILLCELL,'Fill unit cell','Fill unit cell with selected atoms')
        G2G.Define_wxId('wxID_DRAWADDMOLECULE')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWADDMOLECULE,'Complete molecule','Cyclicly add atoms bonded to selected atoms')
        G2G.Define_wxId('wxID_DRAWVOIDMAP')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWVOIDMAP,'Create void map','Create a map of locations outside of any VDW radius')
        self.DrawAtomEdit.Append(G2G.wxID_DRAWDELETE,'Delete atoms','Delete atoms from drawing set')
        G2G.Define_wxId('wxID_RELOADATOMS')
        self.DrawAtomEdit.Append(G2G.wxID_RELOADATOMS,'Reload draw atoms','Reload atom drawing list')
        
        self.DrawAtomCompute.Append(G2G.wxID_DRAWDISTVP,'View pt. dist.','Compute distance of selected atoms from view point')   
        self.DrawAtomCompute.Append(G2G.wxID_DRAWDISAGLTOR,'Dist. Ang. Tors.',
            'Compute distance, angle or torsion for 2-4 selected atoms')   
        self.DrawAtomCompute.Append(G2G.wxID_DRAWPLANE,'Best plane','Compute best plane for 4+ selected atoms')   
        self.DrawAtomRestraint.Append(G2G.wxID_DRAWRESTRBOND,'Add bond restraint','Add bond restraint for selected atoms (2)')
        self.DrawAtomRestraint.Append(G2G.wxID_DRAWRESTRANGLE,'Add angle restraint',
            'Add angle restraint for selected atoms (3: one end 1st)')
        self.DrawAtomRestraint.Append(G2G.wxID_DRAWRESTRPLANE,'Add plane restraint',
            'Add plane restraint for selected atoms (4+)')
        self.DrawAtomRestraint.Append(G2G.wxID_DRAWRESTRCHIRAL,'Add chiral restraint',
            'Add chiral restraint for selected atoms (4: center atom 1st)')
        self.DrawAtomRigidBody.Append(G2G.wxID_DRAWDEFINERB,'Define rigid body','Define rigid body with selected atoms')
        self.PostfillDataMenu()

        # Phase / MCSA tab
        G2G.Define_wxId('wxID_ADDMCSAATOM', 'wxID_ADDMCSARB', 'wxID_CLEARMCSARB', 'wxID_MOVEMCSA', 'wxID_MCSACLEARRESULTS',)        
        self.MCSAMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MCSAMenu)
        self.MCSAMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MCSAEdit = wx.Menu(title='')
        self.MCSAMenu.Append(menu=self.MCSAEdit, title='MC/SA')
        self.MCSAEdit.Append(G2G.wxID_ADDMCSAATOM,'Add atom','Add single atom to MC/SA model')
        self.MCSAEdit.Append(G2G.wxID_ADDMCSARB,'Add rigid body','Add rigid body to MC/SA model' )
        self.MCSAEdit.Append(G2G.wxID_CLEARMCSARB,'Clear rigid bodies','Clear all atoms & rigid bodies from MC/SA model' )
        self.MCSAEdit.Append(G2G.wxID_MOVEMCSA,'Move MC/SA solution','Move MC/SA solution to atom list' )
        self.MCSAEdit.Append(G2G.wxID_MCSACLEARRESULTS,'Clear results','Clear table of MC/SA results' )
        self.PostfillDataMenu()
            
        # Phase / Texture tab
        G2G.Define_wxId('wxID_CLEARTEXTURE', 'wxID_REFINETEXTURE',)        
        self.TextureMenu = wx.MenuBar()
        self.PrefillDataMenu(self.TextureMenu)
        self.TextureMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.TextureEdit = wx.Menu(title='')
        self.TextureMenu.Append(menu=self.TextureEdit, title='Texture')
        self.TextureEdit.Append(G2G.wxID_REFINETEXTURE,'Refine texture','Refine the texture coefficients from sequential results')
        self.PostfillDataMenu()
            
        # Phase / Pawley tab
        G2G.Define_wxId('wxID_PAWLEYLOAD', 'wxID_PAWLEYESTIMATE', 'wxID_PAWLEYUPDATE', 'wxID_PAWLEYSELALL', 
            'wxID_PAWLEYSELNONE','wxID_PAWLEYSELTOGGLE', 'wxID_PAWLEYSET',)
        self.PawleyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PawleyMenu)
        self.PawleyMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PawleyEdit = wx.Menu(title='')
        self.PawleyMenu.Append(menu=self.PawleyEdit,title='Operations')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYSET,'Pawley settings','Change Pawley refinement settings')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYLOAD,'Pawley create','Initialize Pawley reflection list')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYESTIMATE,'Pawley estimate','Estimate initial Pawley intensities')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYUPDATE,'Pawley update','Update negative Pawley intensities with -0.5*Fobs and turn off refinement')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYSELALL,'Refine all','Refine Fsq of all reflections')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYSELNONE,'Refine none','No reflection Fsq refinement')
        self.PawleyEdit.Append(G2G.wxID_PAWLEYSELTOGGLE,'Toggle Selection','Toggle Selection flag for all reflections to opposite setting')
        self.PostfillDataMenu()
            
        # Phase / Map peaks tab
        G2G.Define_wxId('wxID_PEAKSMOVE', 'wxID_PEAKSCLEAR','wxID_PEAKSUNIQUE', 'wxID_PEAKSDELETE','wxID_PEAKSSAVE','wxID_PEAKSDA',
            'wxID_PEAKSDISTVP', 'wxID_PEAKSVIEWPT', 'wxID_FINDEQVPEAKS', 'wxID_SHOWBONDS','wxID_INVERTPEAKS','wxID_ROLLMAP')
        self.MapPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MapPeaksMenu)
        self.MapPeaksMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MapPeaksEdit = wx.Menu(title='')
        self.MapPeaksMenu.Append(menu=self.MapPeaksEdit, title='Map peaks')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSMOVE,'Move peaks','Move selected peaks to atom list')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSVIEWPT,'View point','View point is 1st peak selected')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSDISTVP,'View pt. dist.','Compute distance of selected peaks from view point')   
        self.MapPeaksEdit.Append(G2G.wxID_SHOWBONDS,'Hide bonds','Hide or show bonds between peak positions')   
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSDA,'Calc dist/ang','Calculate distance or angle for selection')
        self.MapPeaksEdit.Append(G2G.wxID_FINDEQVPEAKS,'Equivalent peaks','Find equivalent peaks')
        self.MapPeaksEdit.Append(G2G.wxID_INVERTPEAKS,'Invert peak positions','Inverts map & peak positions')
        self.MapPeaksEdit.Append(G2G.wxID_ROLLMAP,'Roll charge flip map','Roll charge flip map by specified steps')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSUNIQUE,'Unique peaks','Select unique set')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSSAVE,'Save peaks','Save peaks to csv file')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSDELETE,'Delete peaks','Delete selected peaks')
        self.MapPeaksEdit.Append(G2G.wxID_PEAKSCLEAR,'Clear peaks','Clear the map peak list')
        self.PostfillDataMenu()

        # Phase / Rigid bodies tab
        G2G.Define_wxId(
            'wxID_ASSIGNATMS2RB','wxID_GLOBALRESREFINE','wxID_RBREMOVEALL','wxID_AUTOFINDRESRB',
            'wxID_COPYRBPARMS','wxID_GLOBALTHERM',)
        self.RigidBodiesMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodiesMenu)
        self.RigidBodiesMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.RigidBodiesEdit = wx.Menu(title='')
        self.RigidBodiesMenu.Append(menu=self.RigidBodiesEdit, title='Edit Body')
        self.RigidBodiesEdit.Append(G2G.wxID_ASSIGNATMS2RB,'Locate && Insert Rigid Body',
            'Locate rigid body in structure mapping to existing atoms')
        self.RigidBodiesEdit.Append(G2G.wxID_AUTOFINDRESRB,'Auto find residues',
            'Auto find of residue RBs in macromolecule')
        self.RigidBodiesEdit.Append(G2G.wxID_COPYRBPARMS,'Copy rigid body parms',
            'Copy rigid body location & TLS parameters')
        self.RigidBodiesEdit.Append(G2G.wxID_GLOBALTHERM,'Global thermal motion',
            'Global setting of residue thermal motion models')
        self.RigidBodiesEdit.Append(G2G.wxID_GLOBALRESREFINE,'Global residue refine',
            'Global setting of residue RB refinement flags')
        self.RigidBodiesEdit.Append(G2G.wxID_RBREMOVEALL,'Remove all rigid bodies',
            'Remove all rigid body assignment for atoms')
        self.PostfillDataMenu()
    # end of GSAS-II menu definitions
    
################################################################################
#####  Notebook Tree Item editor
################################################################################                  
def UpdateNotebook(G2frame,data):
    '''Called when the data tree notebook entry is selected. Allows for
    editing of the text in that tree entry
    '''
    def OnNoteBook(event):
        event.Skip()
        data = text.GetValue().split('\n')
        G2frame.GPXtree.SetItemPyData(GetGPXtreeItemId(G2frame,G2frame.root,'Notebook'),data)
        if 'nt' not in os.name:
            text.AppendText('\n')
                    
    text = wx.TextCtrl(G2frame.dataWindow,wx.ID_ANY,
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    text.Bind(wx.EVT_TEXT_ENTER,OnNoteBook)
    text.Bind(wx.EVT_KILL_FOCUS,OnNoteBook)
    for line in data:
        text.AppendText(line+"\n")
    text.AppendText('Notebook entry @ '+time.ctime()+"\n")
    G2frame.dataWindow.GetSizer().Add(wx.StaticText(G2frame.dataWindow,-1,' Add notes on project here: '),0)
    G2frame.dataWindow.GetSizer().Add(text,1,wx.ALL|wx.EXPAND)

################################################################################
#####  Comments
################################################################################       
def UpdateComments(G2frame,data):
    '''Place comments into the data window
    '''
    lines = ""
    for line in data:
        if 'phoenix' in wx.version() and hasattr(line,'decode'):
            lines += line.decode('latin-1').rstrip()+'\n'
        else:
            lines += line.rstrip()+'\n'
    try:
        text = wx.StaticText(G2frame.dataWindow,wx.ID_ANY,lines)
    except:
        text = wx.StaticText(G2frame.dataWindow,wx.ID_ANY,
                                 G2obj.StripUnicode(lines))
    G2frame.dataWindow.GetSizer().Add(text,1,wx.ALL|wx.EXPAND)

            
################################################################################
#####  Controls Tree Item editor
################################################################################           
def UpdateControls(G2frame,data):
    '''Edit overall GSAS-II controls in main Controls data tree entry
    '''
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytic Hessian'
        data['min dM/M'] = 0.001
        data['shift factor'] = 1.
        data['max cyc'] = 3        
        data['F**2'] = False
    if 'SVDtol' not in data:
        data['SVDtol'] = 1.e-6
    if 'shift factor' not in data:
        data['shift factor'] = 1.
    if 'max cyc' not in data:
        data['max cyc'] = 3
    if 'F**2' not in data:
        data['F**2'] = False
    if 'Author' not in data:
        data['Author'] = 'no name'
    if 'FreePrm1' not in data:
        data['FreePrm1'] = 'Sample humidity (%)'
    if 'FreePrm2' not in data:
        data['FreePrm2'] = 'Sample voltage (V)'
    if 'FreePrm3' not in data:
        data['FreePrm3'] = 'Applied load (MN)'
    if 'Copy2Next' not in data:
        data['Copy2Next'] = False
    if 'Reverse Seq' not in data:
        data['Reverse Seq'] = False
    if 'UsrReject' not in data:
        data['UsrReject'] = {'minF/sig':0.,'MinExt':0.01,'MaxDF/F':20.,'MaxD':500.,'MinD':0.05}
    if 'HatomFix' not in data:
        data['HatomFix'] = False
    if 'Marquardt' not in data:
        data['Marquardt'] = -3
    
    #end patch

    def SeqSizer():
        
        def OnSelectData(event):
            choices = GetGPXtreeDataNames(G2frame,['PWDR','HKLF',])
            sel = []
            try:
                if 'Seq Data' in data:
#                    for item in data['Seq Data']:
#                        sel.append(choices.index(item))
                    sel = [choices.index(item) for item in data['Seq Data']]
            except ValueError:  #data changed somehow - start fresh
                sel = []
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Sequential refinement',
                'Select dataset to include',choices)
            dlg.SetSelections(sel)
            names = []
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    names.append(choices[sel])
                data['Seq Data'] = names                
            dlg.Destroy()
            G2frame.SetTitleByGPX()
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnReverse(event):
            data['Reverse Seq'] = reverseSel.GetValue()
            
        def OnCopySel(event):
            data['Copy2Next'] = copySel.GetValue()
            
        def OnClrSeq(event):
            sId = GetGPXtreeItemId(G2frame,G2frame.root,'Sequential results')
            if sId:
                dlg = wx.MessageDialog(G2frame,'Are you sure?','Delete sequential results table',wx.OK|wx.CANCEL)
                try:
                    if dlg.ShowModal() == wx.ID_OK:
                        G2frame.GPXtree.Delete(sId)
                finally:
                    dlg.Destroy()
                    
        seqSizer = wx.BoxSizer(wx.VERTICAL)
        dataSizer = wx.BoxSizer(wx.HORIZONTAL)
        dataSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Sequential Refinement: '),0,WACV)
        selSeqData = wx.Button(G2frame.dataWindow,-1,label=' Select data')
        selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
        dataSizer.Add(selSeqData,0,WACV)
        SeqData = data.get('Seq Data',[])
        if not SeqData:
            lbl = ' (no data selected)'
        else:
            lbl = ' ('+str(len(SeqData))+' dataset(s) selected)'

        dataSizer.Add(wx.StaticText(G2frame.dataWindow,label=lbl),0,WACV)
        seqSizer.Add(dataSizer,0)
        if SeqData:
            selSizer = wx.BoxSizer(wx.HORIZONTAL)
            reverseSel = wx.CheckBox(G2frame.dataWindow,-1,label=' Reverse order?')
            reverseSel.Bind(wx.EVT_CHECKBOX,OnReverse)
            reverseSel.SetValue(data['Reverse Seq'])
            selSizer.Add(reverseSel,0,WACV)
            copySel =  wx.CheckBox(G2frame.dataWindow,-1,label=' Copy results to next histogram?')
            copySel.Bind(wx.EVT_CHECKBOX,OnCopySel)
            copySel.SetValue(data['Copy2Next'])
            selSizer.Add(copySel,0,WACV)
            clrSeq = wx.Button(G2frame.dataWindow,label='Clear previous seq.results')
            clrSeq.Bind(wx.EVT_BUTTON,OnClrSeq)
            selSizer.Add(clrSeq,0,WACV)
            seqSizer.Add(selSizer,0)
        return seqSizer
        
    def LSSizer():        
        
        def OnDerivType(event):
            data['deriv type'] = derivSel.GetValue()
            derivSel.SetValue(data['deriv type'])
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnMaxCycles(event):
            data['max cyc'] = int(maxCyc.GetValue())
            maxCyc.SetValue(str(data['max cyc']))
            
        def OnMarqLam(event):
            data['Marquardt'] = int(marqLam.GetValue())
            marqLam.SetValue(str(data['Marquardt']))
                        
        def OnFactor(event):
            event.Skip()
            try:
                value = min(max(float(Factr.GetValue()),0.00001),100.)
            except ValueError:
                value = 1.0
            data['shift factor'] = value
            Factr.SetValue('%.5f'%(value))
            
        def OnFsqRef(event):
            data['F**2'] = fsqRef.GetValue()
            
        LSSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
        LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Refinement derivatives: '),0,WACV)
        Choice=['analytic Jacobian','numeric','analytic Hessian','Hessian SVD']   #TODO +'SVD refine' - what flags will it need?
        derivSel = wx.ComboBox(parent=G2frame.dataWindow,value=data['deriv type'],choices=Choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        derivSel.SetValue(data['deriv type'])
        derivSel.Bind(wx.EVT_COMBOBOX, OnDerivType)
            
        LSSizer.Add(derivSel,0,WACV)
        LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Min delta-M/M: '),0,WACV)
        LSSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'min dM/M',nDig=(10,2,'g'),xmin=1.e-9,xmax=1.),0,WACV)
        if 'Hessian' in data['deriv type']:
            LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Max cycles: '),0,WACV)
            Choice = ['0','1','2','3','5','10','15','20']
            maxCyc = wx.ComboBox(parent=G2frame.dataWindow,value=str(data['max cyc']),choices=Choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            maxCyc.Bind(wx.EVT_COMBOBOX, OnMaxCycles)
            LSSizer.Add(maxCyc,0,WACV)
            if 'SVD' not in data['deriv type']:
                LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Initial lambda = 10**'),0,WACV)
                MarqChoice = ['-3','-2','-1','0','1','2','3','4']
                marqLam = wx.ComboBox(parent=G2frame.dataWindow,value=str(data['Marquardt']),choices=MarqChoice,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                marqLam.Bind(wx.EVT_COMBOBOX,OnMarqLam)
                LSSizer.Add(marqLam,0,WACV)
            LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' SVD zero tolerance:'),0,WACV)
            LSSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'SVDtol',nDig=(10,1,'g'),xmin=1.e-9,xmax=.01),0,WACV)
        else:       #TODO what for SVD refine?
            LSSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Initial shift factor: '),0,WACV)
            Factr = G2G.ValidatedTxtCtrl(G2frame.dataWindow,data,'shift factor',nDig=(10,5),xmin=1.e-5,xmax=100.)
            LSSizer.Add(Factr,0,WACV)
        if G2frame.Sngl:
            userReject = data['UsrReject']
            usrRej = {'minF/sig':[' Min obs/sig (0-5): ',[0.,5.], ],'MinExt':[' Min extinct. (0-.9): ',[0.,.9],],
                'MaxDF/F':[' Max delt-F/sig (3-1000): ',[3.,1000.],],'MaxD':[' Max d-spacing (3-500): ',[3.,500.],],
                'MinD':[' Min d-spacing (0.1-2.0): ',[0.1,2.0],]}

            fsqRef = wx.CheckBox(G2frame.dataWindow,-1,label='Refine HKLF as F^2? ')
            fsqRef.SetValue(data['F**2'])
            fsqRef.Bind(wx.EVT_CHECKBOX,OnFsqRef)
            LSSizer.Add(fsqRef,0,WACV)
            LSSizer.Add((1,0),)
            for item in usrRej:
                LSSizer.Add(wx.StaticText(G2frame.dataWindow,-1,label=usrRej[item][0]),0,WACV)
                usrrej = G2G.ValidatedTxtCtrl(G2frame.dataWindow,userReject,item,nDig=(10,2),
                    xmin=usrRej[item][1][0],xmax=usrRej[item][1][1])
                LSSizer.Add(usrrej,0,WACV)
        return LSSizer
        
    def AuthSizer():
        def OnAuthor(event):
            event.Skip()
            data['Author'] = auth.GetValue()

        Author = data['Author']
        authSizer = wx.BoxSizer(wx.HORIZONTAL)
        authSizer.Add(wx.StaticText(G2frame.dataWindow,label=' CIF Author (last, first):'),0,WACV)
        auth = wx.TextCtrl(G2frame.dataWindow,-1,value=Author,style=wx.TE_PROCESS_ENTER)
        auth.Bind(wx.EVT_TEXT_ENTER,OnAuthor)
        auth.Bind(wx.EVT_KILL_FOCUS,OnAuthor)
        authSizer.Add(auth,0,WACV)
        return authSizer

    def ClearFrozen(event):
        'Removes all frozen parameters by clearing the entire dict'
        Controls['parmFrozen'] = {}
        wx.CallAfter(UpdateControls,G2frame,data)
        
    # start of UpdateControls
    if 'SVD' in data['deriv type']:
        G2frame.GetStatusBar().SetStatusText('Hessian SVD not recommended for initial refinements; use analytic Hessian or Jacobian',1)
    else:
        G2frame.GetStatusBar().SetStatusText('',1)
    G2frame.dataWindow.ClearData()
    SetDataMenuBar(G2frame,G2frame.dataWindow.ControlsMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    subSizer = wx.BoxSizer(wx.HORIZONTAL)
    subSizer.Add((-1,-1),1,wx.EXPAND)
    subSizer.Add(wx.StaticText(G2frame.dataWindow,label='Refinement Controls'),0,WACV)    
    subSizer.Add((-1,-1),1,wx.EXPAND)
    mainSizer.Add(subSizer,0,wx.EXPAND)
    mainSizer.Add((5,5),0)
    mainSizer.Add(LSSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(SeqSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(AuthSizer())
    mainSizer.Add((5,5),0)
    Controls = data
    # count frozen variables (in appropriate place)
    for key in ('parmMinDict','parmMaxDict','parmFrozen'):
        if key not in Controls: Controls[key] = {}
    parmFrozen = Controls['parmFrozen']
    if G2frame.testSeqRefineMode():
        frozenList = set()
        for h in parmFrozen:
            if h == 'FrozenList': continue
            frozenList = frozenList.union(parmFrozen[h])
        count = len(frozenList)
    elif 'FrozenList' in parmFrozen:
        count = len(parmFrozen['FrozenList'])
    else:
        count = 0
    if count > 0:
        subSizer = wx.BoxSizer(wx.HORIZONTAL)
        subSizer.Add(wx.StaticText(G2frame.dataWindow,
            label='There are {} frozen variables (values refined outside limits)'.format(count)),0,WACV)
        subSizer.Add((5,-1))
        btn = wx.Button(G2frame.dataWindow, wx.ID_ANY,'Clear All Frozen')
        btn.Bind(wx.EVT_BUTTON,ClearFrozen)
        subSizer.Add(btn)
        mainSizer.Add(subSizer)
        
    mainSizer.Layout()
    mainSizer.FitInside(G2frame.dataWindow)    
    G2frame.dataWindow.SetSizer(mainSizer)
    G2frame.dataWindow.SetDataSize()
    G2frame.SendSizeEvent()
     
################################################################################
#####  Display of Sequential Results
################################################################################           
       
def UpdateSeqResults(G2frame,data,prevSize=None):
    """
    Called when the Sequential Results data tree entry is selected
    to show results from a sequential refinement.
    
    :param wx.Frame G2frame: main GSAS-II data tree windows

    :param dict data: a dictionary containing the following items:  

            * 'histNames' - list of histogram names in order as processed by Sequential Refinement
            * 'varyList' - list of variables - identical over all refinements in sequence
              note that this is the original list of variables, prior to processing
              constraints.
            * 'variableLabels' -- a dict of labels to be applied to each parameter
              (this is created as an empty dict if not present in data).
            * keyed by histName - dictionaries for all data sets processed, which contains:

              * 'variables'- result[0] from leastsq call
              * 'varyList' - list of variables passed to leastsq call (not same as above)
              * 'sig' - esds for variables
              * 'covMatrix' - covariance matrix from individual refinement
              * 'title' - histogram name; same as dict item name
              * 'newAtomDict' - new atom parameters after shifts applied
              * 'newCellDict' - refined cell parameters after shifts to A0-A5 from Dij terms applied'
    """
    def GetSampleParms():
        '''Make a dictionary of the sample parameters that are not the same over the
        refinement series. Controls here is local
        '''
        if 'IMG' in histNames[0]:
            sampleParmDict = {'Sample load':[],}
        else:
            sampleParmDict = {'Temperature':[],'Pressure':[],'Time':[],
                'FreePrm1':[],'FreePrm2':[],'FreePrm3':[],'Omega':[],
                'Chi':[],'Phi':[],'Azimuth':[],}
        Controls = G2frame.GPXtree.GetItemPyData(
            GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        sampleParm = {}
        for name in histNames:
            if 'IMG' in name:
                if name not in data:
                    continue
                for item in sampleParmDict:
                    sampleParmDict[item].append(data[name]['parmDict'].get(item,0))
            else:
                if 'PDF' in name:
                    name = 'PWDR' + name[4:]
                Id = GetGPXtreeItemId(G2frame,G2frame.root,name)
                if Id:
                    sampleData = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                    for item in sampleParmDict:
                        sampleParmDict[item].append(sampleData.get(item,0))
        for item in sampleParmDict:
            if sampleParmDict[item]:
                frstValue = sampleParmDict[item][0]
                if np.any(np.array(sampleParmDict[item])-frstValue):
                    if item.startswith('FreePrm'):
                        sampleParm[Controls[item]] = sampleParmDict[item]
                    else:
                        sampleParm[item] = sampleParmDict[item]
        return sampleParm

    def GetColumnInfo(col):
        '''returns column label, lists of values and errors (or None) for each column in the table
        for plotting. The column label is reformatted from Unicode to MatPlotLib encoding
        '''
        colName = G2frame.SeqTable.GetColLabelValue(col)
        plotName = variableLabels.get(colName,colName)
        plotName = plotSpCharFix(plotName)
        return plotName,G2frame.colList[col],G2frame.colSigs[col]
            
    def PlotSelectedColRow(calltyp=''):
        '''Called to plot a selected column or row by clicking or 
        double-clicking on a row or column label. N.B. This is called 
        after the event is processed so that the column or row has been
        selected.

        :param str calltyp: ='single'/'double', specifies if this was 
          a single- or double-click, where a single click on row 
          plots histogram; Double click on row plots V-C matrix; 
          Single or double click on column: plots values in column
        '''
        cols = G2frame.dataDisplay.GetSelectedCols()
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif rows and calltyp == 'single':
            name = histNames[rows[0]]       #only does 1st one selected
            if not name.startswith('PWDR'): return 
            pickId = G2frame.PickId
            G2frame.PickId = G2frame.PatternId = GetGPXtreeItemId(G2frame, G2frame.root, name)
            G2plt.PlotPatterns(G2frame,newPlot=False,plotType='PWDR')
            G2frame.PickId = pickId
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
        else:
            G2frame.ErrorDialog(
                'Select row or columns',
                'Nothing selected in table. Click on column or row label(s) to plot. N.B. Grid selection can be a bit funky.'
                )
        
    def PlotSSelect(event):
        'Called by a single click on a row or column label. '
        event.Skip()
        wx.CallAfter(PlotSelectedColRow,'single')
        
    def PlotSelect(event):
        'Called by a double-click on a row or column label'
        event.Skip()
        wx.CallAfter(PlotSelectedColRow,'double')
        
    def OnKeyUp(event):
        event.Skip()
        wx.CallAfter(PlotSelectedColRow,'single')
            
    def OnPlotSelSeq(event):
        'plot the selected columns or row from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
        else:
            G2frame.ErrorDialog('Select columns',
                'No columns or rows selected in table. Click on row or column labels to select fields for plotting.')
                
    def OnAveSelSeq(event):
        'average the selected columns from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        useCol =  ~np.array(G2frame.SeqTable.GetColValues(1),dtype=bool)
        if cols:
            for col in cols:
                items = GetColumnInfo(col)[1]
                noneMask = np.array([item is None for item in items])
                info = ma.array(items,mask=useCol+noneMask)
                ave = ma.mean(ma.compressed(info))
                sig = ma.std(ma.compressed(info))
                print (u' Average for '+G2frame.SeqTable.GetColLabelValue(col)+u': '+'%.6g'%(ave)+u' +/- '+u'%.6g'%(sig))
        else:
            G2frame.ErrorDialog('Select columns',
                'No columns selected in table. Click on column labels to select fields for averaging.')
            
    def OnSelectUse(event):
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select rows to use','Select rows',histNames)
        sellist = [i for i,item in enumerate(G2frame.colList[1]) if item]
        dlg.SetSelections(sellist)
        if dlg.ShowModal() == wx.ID_OK:
            sellist = dlg.GetSelections()
            for row in range(G2frame.SeqTable.GetNumberRows()):
                G2frame.SeqTable.SetValue(row,1,False)
                G2frame.colList[1][row] = False
            for row in sellist:
                G2frame.SeqTable.SetValue(row,1,True)
                G2frame.colList[1][row] = True
            G2frame.dataDisplay.ForceRefresh()
        dlg.Destroy()
                
    def OnRenameSelSeq(event):
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        colNames = [G2frame.SeqTable.GetColLabelValue(c) for c in cols]
        newNames = colNames[:]
        for i,name in enumerate(colNames):
            if name in variableLabels:
                newNames[i] = variableLabels[name]
        if not cols:
            G2frame.ErrorDialog('Select columns',
                'No columns selected in table. Click on column labels to select fields for rename.')
            return
        dlg = G2G.MultiStringDialog(G2frame.dataDisplay,'Set column names',colNames,newNames)
        if dlg.Show():
            newNames = dlg.GetValues()            
            variableLabels.update(dict(zip(colNames,newNames)))
        data['variableLabels'] = variableLabels 
        dlg.Destroy()
        UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
            
    def OnSaveSelSeqCSV(event):
        'export the selected columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True)
        
    def OnSaveSeqCSV(event):
        'export all columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True,allcols=True)
        
    def OnSaveSelSeq(event,csv=False,allcols=False):
        'export the selected columns to a .txt or .csv file from menu command'
        def WriteLine(line):
            if '2' in platform.python_version_tuple()[0]:
                SeqFile.write(G2obj.StripUnicode(line))
            else:
                SeqFile.write(line)
                
        def WriteCSV():
            def WriteList(headerItems):
                line = ''
                for lbl in headerItems:
                    if line: line += ','
                    line += '"'+lbl+'"'
                return line
            head = ['name']
            for col in cols:
                # Excel does not like unicode
                item = G2obj.StripUnicode(G2frame.SeqTable.GetColLabelValue(col))
                if col in havesig:
                    head += [item,'esd-'+item]
                else:
                    head += [item]
            WriteLine(WriteList(head)+'\n')
            for row,name in enumerate(saveNames):
                line = '"'+saveNames[row]+'"'
                for col in cols:
                    if saveData[col][row] is None:
                        if col in havesig:
#                            line += ',0.0,0.0'
                            line += ',,'
                        else:
#                            line += ',0.0'
                            line += ','
                    else:
                        if col in havesig:
                            line += ','+str(saveData[col][row])+','+str(saveSigs[col][row])
                        else:    
                            line += ','+str(saveData[col][row])
                WriteLine(line+'\n')
        def WriteSeq():
            lenName = len(saveNames[0])
            line = '  %s  '%('name'.center(lenName))
            for col in cols:
                item = G2frame.SeqTable.GetColLabelValue(col)
                if col in havesig:
                    line += ' %12s %12s '%(item.center(12),'esd'.center(12))
                else:
                    line += ' %12s '%(item.center(12))
            WriteLine(line+'\n')
            for row,name in enumerate(saveNames):
                line = " '%s' "%(saveNames[row])
                for col in cols:
                    if col in havesig:
                        try:
                            line += ' %12.6f %12.6f '%(saveData[col][row],saveSigs[col][row])
                        except TypeError:
                            line += '                           '
                    else:
                        try:
                            line += ' %12.6f '%saveData[col][row]
                        except TypeError:
                            line += '              '
                WriteLine(line+'\n')

        # start of OnSaveSelSeq code
        if allcols:
            cols = range(G2frame.SeqTable.GetNumberCols())
        else:
            cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        nrows = G2frame.SeqTable.GetNumberRows()
        if not cols:
            choices = [G2frame.SeqTable.GetColLabelValue(r) for r in range(G2frame.SeqTable.GetNumberCols())]
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select columns to write',
                'select columns',choices)
            #dlg.SetSelections()
            if dlg.ShowModal() == wx.ID_OK:
                cols = dlg.GetSelections()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
            #G2frame.ErrorDialog('Select columns',
            #                 'No columns selected in table. Click on column labels to select fields for output.')
            #return
        saveNames = [G2frame.SeqTable.GetRowLabelValue(r) for r in range(nrows)]
        saveData = {}
        saveSigs = {}
        havesig = []
        for col in cols:
            name,vals,sigs = GetColumnInfo(col)
            saveData[col] = vals
            if sigs:
                havesig.append(col)
                saveSigs[col] = sigs
        if csv:
            wild = 'CSV output file (*.csv)|*.csv'
        else:
            wild = 'Text output file (*.txt)|*.txt'
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(
            G2frame,
            'Choose text output file for your selection', pth, '', 
            wild,wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                SeqTextFile = dlg.GetPath()
                SeqTextFile = G2IO.FileDlgFixExt(dlg,SeqTextFile) 
                SeqFile = open(SeqTextFile,'w')
                if csv:
                    WriteCSV()
                else:
                    WriteSeq()
                SeqFile.close()
        finally:
            dlg.Destroy()
                
    def striphist(var,insChar=''):
        'strip a histogram number from a var name'
        sv = var.split(':')
        if len(sv) <= 1: return var
        if sv[1]:
            sv[1] = insChar
        return ':'.join(sv)
        
    def plotSpCharFix(lbl):
        'Change selected unicode characters to their matplotlib equivalent'
        for u,p in [
            (u'\u03B1',r'$\alpha$'),
            (u'\u03B2',r'$\beta$'),
            (u'\u03B3',r'$\gamma$'),
            (u'\u0394\u03C7',r'$\Delta\chi$'),
            ]:
            lbl = lbl.replace(u,p)
        return lbl
    
    def SelectXaxis():
        'returns a selected column number (or None) as the X-axis selection'
        ncols = G2frame.SeqTable.GetNumberCols()
        colNames = [G2frame.SeqTable.GetColLabelValue(r) for r in range(ncols)]
        dlg = G2G.G2SingleChoiceDialog(
            G2frame.dataDisplay,
            'Select x-axis parameter for plot or Cancel for sequence number',
            'Select X-axis',
            colNames)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                col = dlg.GetSelection()
            else:
                col = None
        finally:
            dlg.Destroy()
        return col
    
    def EnablePseudoVarMenus():
        'Enables or disables the PseudoVar menu items that require existing defs'
        if data['SeqPseudoVars']:
            val = True
        else:
            val = False
        G2frame.dataWindow.SequentialPvars.Enable(G2G.wxID_DELSEQVAR,val)
        G2frame.dataWindow.SequentialPvars.Enable(G2G.wxID_EDITSEQVAR,val)

    def DelPseudoVar(event):
        'Ask the user to select a pseudo var expression to delete'
        choices = list(data['SeqPseudoVars'].keys())
        selected = G2G.ItemSelector(
            choices,G2frame,
            multiple=True,
            title='Select expressions to remove',
            header='Delete expression')
        if selected is None: return
        for item in selected:
            del data['SeqPseudoVars'][choices[item]]
        if selected:
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def EditPseudoVar(event):
        'Edit an existing pseudo var expression'
        choices = list(data['SeqPseudoVars'].keys())
        if len(choices) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                choices,G2frame,
                multiple=False,
                title='Select an expression to edit',
                header='Edit expression')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,PSvarDict,
                data['SeqPseudoVars'][choices[selected]],
                header="Edit the PseudoVar expression",
                VarLabel="PseudoVar #"+str(selected+1),
                fit=False)
            newobj = dlg.Show(True)
            if newobj:
                calcobj = G2obj.ExpressionCalcObj(newobj)
                del data['SeqPseudoVars'][choices[selected]]
                data['SeqPseudoVars'][calcobj.eObj.expression] = newobj
                UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        
    def AddNewPseudoVar(event):
        'Create a new pseudo var expression'
        dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,PSvarDict,
            header='Enter an expression for a PseudoVar here',
            VarLabel = "New PseudoVar",fit=False)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            calcobj = G2obj.ExpressionCalcObj(obj)
            data['SeqPseudoVars'][calcobj.eObj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
            
    def AddNewDistPseudoVar(event):
        obj = None
        dlg = G2exG.BondDialog(
            G2frame.dataDisplay,Phases,PSvarDict,
            header='Select a Bond here',
            VarLabel = "New Bond")
        if dlg.ShowModal() == wx.ID_OK:
            pName,Oatom,Tatom = dlg.GetSelection()
            if Tatom:
                Phase = Phases[pName]
                General = Phase['General']
                cx,ct = General['AtomPtrs'][:2]
                pId = Phase['pId']
                SGData = General['SGData']
                sB = Tatom.find('(')+1
                symNo = 0
                if sB:
                    sF = Tatom.find(')')
                    symNo = int(Tatom[sB:sF])
                cellNo = [0,0,0]
                cB = Tatom.find('[')
                if cB>0:
                    cF = Tatom.find(']')+1
                    cellNo = eval(Tatom[cB:cF])
                Atoms = Phase['Atoms']
                aNames = [atom[ct-1] for atom in Atoms]
                oId = aNames.index(Oatom)
                tId = aNames.index(Tatom.split(' +')[0])
                # create an expression object
                obj = G2obj.ExpressionObj()
                obj.expression = 'Dist(%s,\n%s)'%(Oatom,Tatom.split(' d=')[0].replace(' ',''))
                obj.distance_dict = {'pId':pId,'SGData':SGData,'symNo':symNo,'cellNo':cellNo}
                obj.distance_atoms = [oId,tId]
        else: 
            dlg.Destroy()
            return
        dlg.Destroy()
        if obj:
            data['SeqPseudoVars'][obj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def AddNewAnglePseudoVar(event):
        obj = None
        dlg = G2exG.AngleDialog(
            G2frame.dataDisplay,Phases,PSvarDict,
            header='Enter an Angle here',
            VarLabel = "New Angle")
        if dlg.ShowModal() == wx.ID_OK:
            pName,Oatom,Tatoms = dlg.GetSelection()
            if Tatoms:
                Phase = Phases[pName]
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
                Tatoms = Tatoms.split(';')
                for Tatom in Tatoms:
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
                obj = G2obj.ExpressionObj()
                obj.expression = 'Angle(%s,%s,\n%s)'%(Tatoms[0],Oatom,Tatoms[1])
                obj.angle_dict = {'pId':pId,'SGData':SGData,'symNo':symNos,'cellNo':cellNos}
                obj.angle_atoms = [oId,tIds]
        else: 
            dlg.Destroy()
            return
        dlg.Destroy()
        if obj:
            data['SeqPseudoVars'][obj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
            
    def UpdateParmDict(parmDict):
        '''generate the atom positions and the direct & reciprocal cell values,
        because they might be needed to evaluate the pseudovar
        '''
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        delList = []
        phaselist = []
        for item in parmDict: 
            if ':' not in item: continue
            key = item.split(':')
            if len(key) < 3: continue
            # remove the dA[xyz] terms, they would only bring confusion
            if key[0] and key[0] not in phaselist: phaselist.append(key[0])
            if key[2].startswith('dA'):
                delList.append(item)
            # compute and update the corrected reciprocal cell terms using the Dij values
            elif key[2] in Ddict:
                akey = key[0]+'::'+Ddict[key[2]]
                parmDict[akey] += parmDict[item]
                delList.append(item)
        for item in delList:
            del parmDict[item]                
        for i in phaselist:
            pId = int(i)
            # apply cell symmetry
            A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],parmDict,zeroDict[pId])
            # convert to direct cell & add the unique terms to the dictionary
            for i,val in enumerate(G2lat.A2cell(A)):
                if i in uniqCellIndx[pId]:
                    lbl = str(pId)+'::'+cellUlbl[i]
                    parmDict[lbl] = val
            lbl = str(pId)+'::'+'Vol'
            parmDict[lbl] = G2lat.calc_V(A)
        return parmDict

    def EvalPSvarDeriv(calcobj,parmDict,sampleDict,var,ESD):
        '''Evaluate an expression derivative with respect to a
        GSAS-II variable name.

        Note this likely could be faster if the loop over calcobjs were done
        inside after the Dict was created. 
        '''
        if not ESD:
            return 0.
        step = ESD/10
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        results = []
        phaselist = []
        VparmDict = sampleDict.copy()
        for incr in step,-step:
            VparmDict.update(parmDict.copy())           
            # as saved, the parmDict has updated 'A[xyz]' values, but 'dA[xyz]'
            # values are not zeroed: fix that!
            VparmDict.update({item:0.0 for item in parmDict if 'dA' in item})
            VparmDict[var] += incr
            G2mv.Dict2Map(VparmDict,[]) # apply constraints
            # generate the atom positions and the direct & reciprocal cell values now, because they might
            # needed to evaluate the pseudovar
            for item in VparmDict:
                if item in sampleDict:
                    continue 
                if ':' not in item: continue
                key = item.split(':')
                if len(key) < 3: continue
                # apply any new shifts to atom positions
                if key[2].startswith('dA'):
                    VparmDict[''.join(item.split('d'))] += VparmDict[item]
                    VparmDict[item] = 0.0
                # compute and update the corrected reciprocal cell terms using the Dij values
                if key[2] in Ddict:
                    if key[0] not in phaselist: phaselist.append(key[0])
                    akey = key[0]+'::'+Ddict[key[2]]
                    VparmDict[akey] += VparmDict[item]
            for i in phaselist:
                pId = int(i)
                # apply cell symmetry
                A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],VparmDict,zeroDict[pId])
                # convert to direct cell & add the unique terms to the dictionary
                for i,val in enumerate(G2lat.A2cell(A)):
                    if i in uniqCellIndx[pId]:
                        lbl = str(pId)+'::'+cellUlbl[i]
                        VparmDict[lbl] = val
                lbl = str(pId)+'::'+'Vol'
                VparmDict[lbl] = G2lat.calc_V(A)
            # dict should be fully updated, use it & calculate
            calcobj.SetupCalc(VparmDict)
            results.append(calcobj.EvalExpression())
        if None in results:
            return None
        return (results[0] - results[1]) / (2.*step)
        
    def EnableParFitEqMenus():
        'Enables or disables the Parametric Fit menu items that require existing defs'
        if data['SeqParFitEqList']:
            val = True
        else:
            val = False
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_DELPARFIT,val)
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_EDITPARFIT,val)
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_DOPARFIT,val)

    def ParEqEval(Values,calcObjList,varyList):
        '''Evaluate the parametric expression(s)
        :param list Values: a list of values for each variable parameter
        :param list calcObjList: a list of :class:`GSASIIobj.ExpressionCalcObj`
          expression objects to evaluate
        :param list varyList: a list of variable names for each value in Values
        '''
        result = []
        for calcobj in calcObjList:
            calcobj.UpdateVars(varyList,Values)
            if calcobj.depSig:
                result.append((calcobj.depVal-calcobj.EvalExpression())/calcobj.depSig)
            else:
                result.append(calcobj.depVal-calcobj.EvalExpression())
        return result

    def DoParEqFit(event,eqObj=None):
        'Parametric fit minimizer'
        varyValueDict = {} # dict of variables and their initial values
        calcObjList = [] # expression objects, ready to go for each data point
        if eqObj is not None:
            eqObjList = [eqObj,]
        else:
            eqObjList = data['SeqParFitEqList']
        UseFlags = G2frame.SeqTable.GetColValues(1)
        for obj in eqObjList:
            # assemble refined vars for this equation
            varyValueDict.update({var:val for var,val in obj.GetVariedVarVal()})
            # lookup dependent var position
            depVar = obj.GetDepVar()
            if depVar in colLabels:
                indx = colLabels.index(depVar)
            else:
                raise Exception('Dependent variable '+depVar+' not found')
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()
            # loop over each datapoint
            for j,row in enumerate(zip(*G2frame.colList)):
                if not UseFlags[j]: continue
                # assemble equations to fit
                calcobj = G2obj.ExpressionCalcObj(obj)
                # prepare a dict of needed independent vars for this expression
                indepVarDict = {var:row[i] for i,var in enumerate(colLabels) if var in indepVars}
                calcobj.SetupCalc(indepVarDict)                
                # values and sigs for current value of dependent var
                if row[indx] is None: continue
                calcobj.depVal = row[indx]
                calcobj.depSig = G2frame.colSigs[indx][j]
                calcObjList.append(calcobj)
        # varied parameters
        varyList = varyValueDict.keys()
        values = varyValues = [varyValueDict[key] for key in varyList]
        if not varyList:
            print ('no variables to refine!')
            return
        try:
            result = so.leastsq(ParEqEval,varyValues,full_output=True,   #ftol=Ftol,
                args=(calcObjList,varyList))
            values = result[0]
            covar = result[1]
            if covar is None:
                raise Exception
            chisq = np.sum(result[2]['fvec']**2)
            GOF = np.sqrt(chisq/(len(calcObjList)-len(varyList)))
            esdDict = {}
            for i,avar in enumerate(varyList):
                esdDict[avar] = np.sqrt(covar[i,i])
        except:
            print('====> Fit failed')
            return
        print('==== Fit Results ====')
        print ('  chisq =  %.2f, GOF = %.2f'%(chisq,GOF))
        for obj in eqObjList:
            obj.UpdateVariedVars(varyList,values)
            ind = '      '
            print(u'  '+obj.GetDepVar()+u' = '+obj.expression)
            for var in obj.assgnVars:
                print(ind+var+u' = '+obj.assgnVars[var])
            for var in obj.freeVars:
                avar = "::"+obj.freeVars[var][0]
                val = obj.freeVars[var][1]
                if obj.freeVars[var][2]:
                    print(ind+var+u' = '+avar + " = " + G2mth.ValEsd(val,esdDict[avar]))
                else:
                    print(ind+var+u' = '+avar + u" =" + G2mth.ValEsd(val,0))
        # create a plot for each parametric variable
        for fitnum,obj in enumerate(eqObjList):
            calcobj = G2obj.ExpressionCalcObj(obj)
            # lookup dependent var position
            indx = colLabels.index(obj.GetDepVar())
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()            
            # loop over each datapoint
            fitvals = []
            for j,row in enumerate(zip(*G2frame.colList)):
                calcobj.SetupCalc({var:row[i] for i,var in enumerate(colLabels) if var in indepVars})
                fitvals.append(calcobj.EvalExpression())
            G2plt.PlotSelectedSequence(G2frame,[indx],GetColumnInfo,SelectXaxis,fitnum,fitvals)

    def SingleParEqFit(eqObj):
        DoParEqFit(None,eqObj)

    def DelParFitEq(event):
        'Ask the user to select function to delete'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        selected = G2G.ItemSelector(
            txtlst,G2frame,
            multiple=True,
            title='Select a parametric equation(s) to remove',
            header='Delete equation')
        if selected is None: return
        data['SeqParFitEqList'] = [obj for i,obj in enumerate(data['SeqParFitEqList']) if i not in selected]
        EnableParFitEqMenus()
        if data['SeqParFitEqList']: DoParEqFit(event)
        
    def EditParFitEq(event):
        'Edit an existing parametric equation'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame,
                multiple=False,
                title='Select a parametric equation to edit',
                header='Edit equation')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,VarDict,
                data['SeqParFitEqList'][selected],depVarDict=VarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit],wildCard=False)
            newobj = dlg.Show(True)
            if newobj:
                data['SeqParFitEqList'][selected] = newobj
                EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)

    def AddNewParFitEq(event):
        'Create a new parametric equation to be fit to sequential results'

        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in data['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])

        dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,VarDict,depVarDict=VarDict,
            header='Define an equation to minimize in the parametric fit',
            ExtraButton=['Fit',SingleParEqFit],usedVars=usedvarlist,wildCard=False)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            data['SeqParFitEqList'].append(obj)
            EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)
                
    def CopyParFitEq(event):
        'Copy an existing parametric equation to be fit to sequential results'
        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in data['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame,
                multiple=False,
                title='Select a parametric equation to copy',
                header='Copy equation')
        if selected is not None:
            newEqn = copy.deepcopy(data['SeqParFitEqList'][selected])
            for var in newEqn.freeVars:
                newEqn.freeVars[var][0] = G2obj.MakeUniqueLabel(newEqn.freeVars[var][0],usedvarlist)
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,VarDict,newEqn,depVarDict=VarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit],wildCard=False)
            newobj = dlg.Show(True)
            if newobj:
                data['SeqParFitEqList'].append(newobj)
                EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)
                                            
    def GridSetToolTip(row,col):
        '''Routine to show standard uncertainties for each element in table
        as a tooltip
        '''
        if G2frame.colSigs[col]:
            if G2frame.colSigs[col][row] == -0.1: return 'frozen'
            return u'\u03c3 = '+str(G2frame.colSigs[col][row])
        return ''
        
    def GridColLblToolTip(col):
        '''Define a tooltip for a column. This will be the user-entered value
        (from data['variableLabels']) or the default name
        '''
        if col < 0 or col > len(colLabels):
            print ('Illegal column #%d'%col)
            return
        var = colLabels[col]
        return variableLabels.get(var,G2obj.fmtVarDescr(var))
        
    def SetLabelString(event):
        '''Define or edit the label for a column in the table, to be used
        as a tooltip and for plotting
        '''
        col = event.GetCol()
        if col < 0 or col > len(colLabels):
            return
        var = colLabels[col]
        lbl = variableLabels.get(var,G2obj.fmtVarDescr(var))
        head = u'Set a new name for variable {} (column {})'.format(var,col)
        dlg = G2G.SingleStringDialog(G2frame,'Set variable label',
                                 head,lbl,size=(400,-1))
        if dlg.Show():
            variableLabels[var] = dlg.GetValue()
            dlg.Destroy()
            wx.CallAfter(UpdateSeqResults,G2frame,data) # redisplay variables
        else:
            dlg.Destroy()

    def DoSequentialExport(event):
        '''Event handler for all Sequential Export menu items
        '''
        vals = G2frame.dataWindow.SeqExportLookup.get(event.GetId())
        if vals is None:
            print('Error: Id not found. This should not happen!')
            return
        G2IO.ExportSequential(G2frame,data,*vals)

    def onSelectSeqVars(event):
        '''Select which variables will be shown in table'''
        hides = [saveColLabels[1:].index(item) for item in G2frame.SeqTblHideList if
                     item in saveColLabels[1:]]
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select columns to hide',
                'Hide columns',saveColLabels[1:])
        dlg.SetSelections(hides)
        if dlg.ShowModal() == wx.ID_OK:
            G2frame.SeqTblHideList = [saveColLabels[1:][sel] for sel in dlg.GetSelections()]
            dlg.Destroy()
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        else:
            dlg.Destroy()
            
    def OnCellChange(event):
        r = event.GetRow()
        val = G2frame.SeqTable.GetValue(r,0)
#        print (r,val)
        G2frame.SeqTable.SetValue(r,0, val)
        
    def OnSelectUpdate(event):
        '''Update all phase parameters from a selected column in the Sequential Table. 
        If no histogram is selected (or more than one), ask the user to make a selection.

        Loosely based on :func:`GSASIIstrIO.SetPhaseData`
        '''
        rows = G2frame.dataDisplay.GetSelectedRows()
        if len(rows) == 1:
            sel = rows[0]
        else:
            dlg = G2G.G2SingleChoiceDialog(G2frame, 'Select histogram to use for update',
                                           'Select row',histNames)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
        parmDict = data[histNames[sel]]['parmDict']
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        for phase in Phases:
            print('Updating {} from Seq. Ref. row {}'.format(phase,histNames[sel]))
            Phase = Phases[phase]
            General = Phase['General']
            SGData = General['SGData']
            Atoms = Phase['Atoms']
            cell = General['Cell']
            pId = Phase['pId']
            pfx = str(pId)+'::'
            # there should not be any changes to the cell because those terms are not refined
            A,sigA = G2stIO.cellFill(pfx,SGData,parmDict,{})
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            textureData = General['SH Texture']    
            if textureData['Order']:
#                SHtextureSig = {}
                for name in ['omega','chi','phi']:
                    aname = pfx+'SH '+name
                    textureData['Sample '+name][1] = parmDict[aname]
                for name in textureData['SH Coeff'][1]:
                    aname = pfx+name
                    textureData['SH Coeff'][1][name] = parmDict[aname]
            ik = 6  #for Pawley stuff below
            if General.get('Modulated',False):
                ik = 7
            # how are these updated?
            #General['SuperVec']
            #RBModels = Phase['RBModels']
            if Phase['General'].get('doPawley'):
                pawleyRef = Phase['Pawley ref']
                for i,refl in enumerate(pawleyRef):
                    key = pfx+'PWLref:'+str(i)
                    refl[ik] = parmDict[key]
#                    if key in sigDict:  #TODO: error here sigDict not defined. What was intended
#                        refl[ik+1] = sigDict[key]
#                    else:
#                        refl[ik+1] = 0
                continue
            General['Mass'] = 0.
            cx,ct,cs,cia = General['AtomPtrs']
            for i,at in enumerate(Atoms):
                names = {cx:pfx+'Ax:'+str(i),cx+1:pfx+'Ay:'+str(i),cx+2:pfx+'Az:'+str(i),cx+3:pfx+'Afrac:'+str(i),
                    cia+1:pfx+'AUiso:'+str(i),cia+2:pfx+'AU11:'+str(i),cia+3:pfx+'AU22:'+str(i),cia+4:pfx+'AU33:'+str(i),
                    cia+5:pfx+'AU12:'+str(i),cia+6:pfx+'AU13:'+str(i),cia+7:pfx+'AU23:'+str(i),
                    cx+4:pfx+'AMx:'+str(i),cx+5:pfx+'AMy:'+str(i),cx+6:pfx+'AMz:'+str(i)}
                for ind in range(cx,cx+4):
                    at[ind] = parmDict[names[ind]]
                if at[cia] == 'I':
                    at[cia+1] = parmDict[names[cia+1]]
                else:
                    for ind in range(cia+2,cia+8):
                        at[ind] = parmDict[names[ind]]
                if General['Type'] == 'magnetic':
                    for ind in range(cx+4,cx+7):
                        at[ind] = parmDict[names[ind]]
                ind = General['AtomTypes'].index(at[ct])
                General['Mass'] += General['AtomMass'][ind]*at[cx+3]*at[cx+5]
                if General.get('Modulated',False):
                    AtomSS = at[-1]['SS1']
                    waveType = AtomSS['waveType']
                    for Stype in ['Sfrac','Spos','Sadp','Smag']:
                        Waves = AtomSS[Stype]
                        for iw,wave in enumerate(Waves):
                            stiw = str(i)+':'+str(iw)
                            if Stype == 'Spos':
                                if waveType in ['ZigZag','Block',] and not iw:
                                    names = ['Tmin:'+stiw,'Tmax:'+stiw,'Xmax:'+stiw,'Ymax:'+stiw,'Zmax:'+stiw]
                                else:
                                    names = ['Xsin:'+stiw,'Ysin:'+stiw,'Zsin:'+stiw,
                                        'Xcos:'+stiw,'Ycos:'+stiw,'Zcos:'+stiw]
                            elif Stype == 'Sadp':
                                names = ['U11sin:'+stiw,'U22sin:'+stiw,'U33sin:'+stiw,
                                    'U12sin:'+stiw,'U13sin:'+stiw,'U23sin:'+stiw,
                                    'U11cos:'+stiw,'U22cos:'+stiw,'U33cos:'+stiw,
                                    'U12cos:'+stiw,'U13cos:'+stiw,'U23cos:'+stiw]
                            elif Stype == 'Sfrac':
                                if 'Crenel' in waveType and not iw:
                                    names = ['Fzero:'+stiw,'Fwid:'+stiw]
                                else:
                                    names = ['Fsin:'+stiw,'Fcos:'+stiw]
                            elif Stype == 'Smag':
                                names = ['MXsin:'+stiw,'MYsin:'+stiw,'MZsin:'+stiw,
                                    'MXcos:'+stiw,'MYcos:'+stiw,'MZcos:'+stiw]
                            for iname,name in enumerate(names):
                                AtomSS[Stype][iw][0][iname] = parmDict[pfx+name]
    
    # lookup table for unique cell parameters by symmetry
    cellGUIlist = [
        [['m3','m3m'],(0,)],
        [['3R','3mR'],(0,3)],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],(0,2)],
        [['mmm'],(0,1,2)],
        [['2/m'+'a'],(0,1,2,3)],
        [['2/m'+'b'],(0,1,2,4)],
        [['2/m'+'c'],(0,1,2,5)],
        [['-1'],(0,1,2,3,4,5)],
        ]
    # cell labels
    cellUlbl = ('a','b','c',u'\u03B1',u'\u03B2',u'\u03B3') # unicode a,b,c,alpha,beta,gamma

    #======================================================================
    # start processing sequential results here (UpdateSeqResults)
    #======================================================================
    if not data:
        print ('No sequential refinement results')
        return
    variableLabels = data.get('variableLabels',{})
    data['variableLabels'] = variableLabels
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    Controls = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,G2frame.root,'Controls'))
    # create a place to store Pseudo Vars & Parametric Fit functions, if not present
    if 'SeqPseudoVars' not in data: data['SeqPseudoVars'] = {}
    if 'SeqParFitEqList' not in data: data['SeqParFitEqList'] = []
    histNames = data['histNames']
    foundNames = [name for name in histNames if name in data]
    histNames = foundNames
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.GetStatusBar().SetStatusText("Select column to export; Double click on column to plot data; on row for Covariance",1)
    sampleParms = GetSampleParms()

    # make dict of varied atom coords keyed by absolute position
    newAtomDict = data[histNames[0]].get('newAtomDict',{}) # dict with atom positions; relative & absolute
    # Possible error: the next might need to be data[histNames[0]]['varyList']
    # error will arise if there constraints on coordinates?
    atomLookup = {newAtomDict[item][0]:item for item in newAtomDict if item in data['varyList']}
    
    # make dict of varied cell parameters equivalents
    ESDlookup = {} # provides the Dij term for each Ak term (where terms are refined)
    Dlookup = {} # provides the Ak term for each Dij term (where terms are refined)
    # N.B. These Dij vars are missing a histogram #
    newCellDict = {}
    for name in histNames:
        if name in data and 'newCellDict' in data[name]:
            newCellDict.update(data[name]['newCellDict'])
    cellAlist = []
    for item in newCellDict:
        cellAlist.append(newCellDict[item][0])
        if item in data.get('varyList',[]):
            ESDlookup[newCellDict[item][0]] = item
            Dlookup[item] = newCellDict[item][0]
    # add coordinate equivalents to lookup table
    for parm in atomLookup:
        Dlookup[atomLookup[parm]] = parm
        ESDlookup[parm] = atomLookup[parm]

    # get unit cell & symmetry for all phases & initial stuff for later use
    RecpCellTerms = {}
    SGdata = {}
    uniqCellIndx = {}
    initialCell = {}
    RcellLbls = {}
    zeroDict = {}
    for phase in Phases:
        phasedict = Phases[phase]
        pId = phasedict['pId']
        pfx = str(pId)+'::' # prefix for A values from phase
        RcellLbls[pId] = [pfx+'A'+str(i) for i in range(6)]
        RecpCellTerms[pId] = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        zeroDict[pId] = dict(zip(RcellLbls[pId],6*[0.,]))
        SGdata[pId] = phasedict['General']['SGData']
        laue = SGdata[pId]['SGLaue']
        if laue == '2/m':
            laue += SGdata[pId]['SGUniq']
        for symlist,celllist in cellGUIlist:
            if laue in symlist:
                uniqCellIndx[pId] = celllist
                break
        else: # should not happen
            uniqCellIndx[pId] = list(range(6))
        for i in uniqCellIndx[pId]:
            initialCell[str(pId)+'::A'+str(i)] =  RecpCellTerms[pId][i]

    SetDataMenuBar(G2frame,G2frame.dataWindow.SequentialMenu)
    G2frame.Bind(wx.EVT_MENU, OnSelectUse, id=G2G.wxID_SELECTUSE)
    G2frame.Bind(wx.EVT_MENU, OnRenameSelSeq, id=G2G.wxID_RENAMESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=G2G.wxID_SAVESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSaveSelSeqCSV, id=G2G.wxID_SAVESEQSELCSV)
    G2frame.Bind(wx.EVT_MENU, OnSaveSeqCSV, id=G2G.wxID_SAVESEQCSV)
    G2frame.Bind(wx.EVT_MENU, OnPlotSelSeq, id=G2G.wxID_PLOTSEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnAveSelSeq, id=G2G.wxID_AVESEQSEL)
    #G2frame.Bind(wx.EVT_MENU, OnReOrgSelSeq, id=G2G.wxID_ORGSEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSelectUpdate, id=G2G.wxID_UPDATESEQSEL)
    G2frame.Bind(wx.EVT_MENU, onSelectSeqVars, id=G2G.wxID_ORGSEQINC)
    G2frame.Bind(wx.EVT_MENU, AddNewPseudoVar, id=G2G.wxID_ADDSEQVAR)
    G2frame.Bind(wx.EVT_MENU, AddNewDistPseudoVar, id=G2G.wxID_ADDSEQDIST)
    G2frame.Bind(wx.EVT_MENU, AddNewAnglePseudoVar, id=G2G.wxID_ADDSEQANGLE)
    G2frame.Bind(wx.EVT_MENU, DelPseudoVar, id=G2G.wxID_DELSEQVAR)
    G2frame.Bind(wx.EVT_MENU, EditPseudoVar, id=G2G.wxID_EDITSEQVAR)
    G2frame.Bind(wx.EVT_MENU, AddNewParFitEq, id=G2G.wxID_ADDPARFIT)
    G2frame.Bind(wx.EVT_MENU, CopyParFitEq, id=G2G.wxID_COPYPARFIT)
    G2frame.Bind(wx.EVT_MENU, DelParFitEq, id=G2G.wxID_DELPARFIT)
    G2frame.Bind(wx.EVT_MENU, EditParFitEq, id=G2G.wxID_EDITPARFIT)
    G2frame.Bind(wx.EVT_MENU, DoParEqFit, id=G2G.wxID_DOPARFIT)

    for id in G2frame.dataWindow.SeqExportLookup:        
        G2frame.Bind(wx.EVT_MENU, DoSequentialExport, id=id)

    EnablePseudoVarMenus()
    EnableParFitEqMenus()

    # scan for locations where the variables change
    VaryListChanges = [] # histograms where there is a change
    combinedVaryList = []
    firstValueDict = {}
    vallookup = {}
    posdict = {}
    prevVaryList = []
    foundNames = []
    missing = 0
    for i,name in enumerate(histNames):
        if name not in data:
            if missing < 5:
                print(" Warning: "+name+" not found")
            elif missing == 5:
                print (' Warning: more are missing')
            missing += 1
            continue
        foundNames.append(name)
        maxPWL = 5
        for var,val,sig in zip(data[name]['varyList'],data[name]['variables'],data[name]['sig']):
            svar = striphist(var,'*') # wild-carded
            if 'PWL' in svar:
                if int(svar.split(':')[-1]) > maxPWL:
                    continue
            if svar not in combinedVaryList:
                # add variables to list as they appear
                combinedVaryList.append(svar)
                firstValueDict[svar] = (val,sig)
        if prevVaryList != data[name]['varyList']: # this refinement has a different refinement list from previous
            prevVaryList = data[name]['varyList']
            vallookup[name] = dict(zip(data[name]['varyList'],data[name]['variables']))
            posdict[name] = {}
            for var in data[name]['varyList']:
                svar = striphist(var,'*')
                if 'PWL' in svar:
                    if int(svar.split(':')[-1]) > maxPWL:
                        continue
                posdict[name][combinedVaryList.index(svar)] = svar
            VaryListChanges.append(name)
    if missing:
        print (' Warning: Total of %d data sets missing from sequential results'%(missing))
    #if len(VaryListChanges) > 1:
    #    G2frame.dataWindow.SequentialFile.Enable(G2G.wxID_ORGSEQSEL,True)
    #else:
    #    G2frame.dataWindow.SequentialFile.Enable(G2G.wxID_ORGSEQSEL,False)
    #-----------------------------------------------------------------------------------
    # build up the data table by columns -----------------------------------------------
    histNames = foundNames
    nRows = len(histNames)
    G2frame.colList = [list(range(nRows)),nRows*[True]]
    G2frame.colSigs = [None,None,]
    colLabels = ['No.','Use',]
    Types = [wg.GRID_VALUE_LONG,wg.GRID_VALUE_BOOL,]
    # start with Rwp values
    if 'IMG ' not in histNames[0][:4]:
        G2frame.colList += [[data[name]['Rvals']['Rwp'] for name in histNames]]
        G2frame.colSigs += [None]
        colLabels += ['Rwp']
        Types += [wg.GRID_VALUE_FLOAT+':10,3',]
    # add % change in Chi^2 in last cycle
    if histNames[0][:4] not in ['SASD','IMG ','REFD'] and Controls.get('ShowCell'):
        G2frame.colList += [[100.*data[name]['Rvals'].get('DelChi2',-1) for name in histNames]]
        G2frame.colSigs += [None]
        colLabels += [u'\u0394\u03C7\u00B2 (%)']
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
    deltaChiCol = len(colLabels)-1
    # frozen variables?
    if 'parmFrozen' in Controls:
        f = [len(Controls['parmFrozen'].get(h,[])) for h in histNames]
        if any(f):
            G2frame.colList += [f]
            G2frame.colSigs += [None]
            colLabels += ['frozen']
            Types += [wg.GRID_VALUE_LONG]
    # add changing sample parameters to table
    for key in sampleParms:
        G2frame.colList += [sampleParms[key]]
        G2frame.colSigs += [None]
        colLabels += [key]
        Types += [wg.GRID_VALUE_FLOAT,]
    sampleDict = {}
    for i,name in enumerate(histNames):
        sampleDict[name] = dict(zip(sampleParms.keys(),[sampleParms[key][i] for key in sampleParms.keys()])) 
    # add unique cell parameters  
    if Controls.get('ShowCell',False) and len(newCellDict):
        phaseLookup = {Phases[phase]['pId']:phase for phase in Phases}
        for pId in sorted(RecpCellTerms):
            pfx = str(pId)+'::' # prefix for A values from phase
            cells = []
            cellESDs = []
            colLabels += [pfx+cellUlbl[i] for i in uniqCellIndx[pId]]
            colLabels += [pfx+'Vol']
            Types += (len(uniqCellIndx[pId]))*[wg.GRID_VALUE_FLOAT+':10,5',]
            Types += [wg.GRID_VALUE_FLOAT+':10,3',]
            Albls = [pfx+'A'+str(i) for i in range(6)]
            for name in histNames:
                hId = Histograms[name]['hId']
                phfx = '%d:%d:'%(pId,hId)
                esdLookUp = {}
                dLookup = {}
                for item in data[name]['newCellDict']:
                    if phfx+item.split('::')[1] in data[name]['varyList']:
                        esdLookUp[newCellDict[item][0]] = item
                        dLookup[item] = newCellDict[item][0]
                covData = {'varyList': [dLookup.get(striphist(v),v) for v in data[name]['varyList']],
                    'covMatrix': data[name]['covMatrix']}
                A = RecpCellTerms[pId][:] # make copy of starting A values
                # update with refined values
                for i,j in enumerate(('D11','D22','D33','D12','D13','D23')):
                    var = str(pId)+'::A'+str(i)
                    Dvar = str(pId)+':'+str(hId)+':'+j
                    # apply Dij value if non-zero
                    if Dvar in data[name]['parmDict']:
                        parmDict = data[name]['parmDict']
                        if parmDict[Dvar] != 0.0:
                            A[i] += parmDict[Dvar]
                    # override with fit result if is Dij varied
                    if var in cellAlist:
                        try:
                            A[i] = data[name]['newCellDict'][esdLookUp[var]][1] # get refined value 
                        except KeyError:
                            pass
                # apply symmetry
                cellDict = dict(zip(Albls,A))
                try:    # convert to direct cell
                    A,zeros = G2stIO.cellFill(pfx,SGdata[pId],cellDict,zeroDict[pId])
                    c = G2lat.A2cell(A)
                    vol = G2lat.calc_V(A)
                    cE = G2stIO.getCellEsd(pfx,SGdata[pId],A,covData)
                except:
                    c = 6*[None]
                    cE = 6*[None]
                    vol = None
                # add only unique values to table
                if name in Phases[phaseLookup[pId]]['Histograms']:
                    cells += [[c[i] for i in uniqCellIndx[pId]]+[vol]]
                    cellESDs += [[cE[i] for i in uniqCellIndx[pId]]+[cE[-1]]]
                else:
                    cells += [[None for i in uniqCellIndx[pId]]+[None]]
                    cellESDs += [[None for i in uniqCellIndx[pId]]+[None]]
            G2frame.colList += zip(*cells)
            G2frame.colSigs += zip(*cellESDs)
    # sort out the variables in their selected order
    varcols = 0
    for d in posdict.values():
        varcols = max(varcols,max(d.keys())+1)
    # get labels for each column
    for i in range(varcols):
        lbl = ''
        for h in VaryListChanges:
            if posdict[h].get(i):
                if posdict[h].get(i) in lbl: continue
                if lbl != "": lbl += '/'
                lbl += posdict[h].get(i)
        colLabels.append(lbl)
    Types += varcols*[wg.GRID_VALUE_FLOAT,]
    vals = []
    esds = []
    varsellist = None        # will be a list of variable names in the order they are selected to appear
    # tabulate values for each hist, leaving None for blank columns
    for name in histNames:
        if name in posdict:
            varsellist = [posdict[name].get(i) for i in range(varcols)]
            # translate variable names to how they will be used in the headings
            vs = [striphist(v,'*') for v in data[name]['varyList']]
            # determine the index for each column (or None) in the data[]['variables'] and ['sig'] lists
            sellist = [vs.index(v) if v is not None else None for v in varsellist]
            #sellist = [i if striphist(v,'*') in varsellist else None for i,v in enumerate(data[name]['varyList'])]
        if not varsellist: raise Exception()
        vals.append([data[name]['variables'][s] if s is not None else None for s in sellist])
        esds.append([data[name]['sig'][s] if s is not None else None for s in sellist])
    G2frame.colList += zip(*vals)
    G2frame.colSigs += zip(*esds)
    # tabulate constrained variables, removing histogram numbers if needed
    # from parameter label
    depValDict = {}
    depSigDict = {}
    for name in histNames:
        for var in data[name].get('depParmDict',{}):
            val,sig = data[name]['depParmDict'][var]
            svar = striphist(var,'*')
            if svar not in depValDict:
               depValDict[svar] = [val]
               depSigDict[svar] = [sig]
            else:
               depValDict[svar].append(val)
               depSigDict[svar].append(sig)
    
    # add the dependent constrained variables to the table
    for var in sorted(depValDict):
        if len(depValDict[var]) != len(histNames): continue
        colLabels.append(var)
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
        G2frame.colSigs += [depSigDict[var]]
        G2frame.colList += [depValDict[var]]

    # add refined atom parameters to table
    colLabels += sorted(atomLookup.keys())
    for parm in sorted(atomLookup):
        G2frame.colList += [[data[name]['newAtomDict'][atomLookup[parm]][1] for name in histNames]]
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
        if atomLookup[parm] in data[histNames[0]]['varyList']:
            col = data[histNames[0]]['varyList'].index(atomLookup[parm])
            G2frame.colSigs += [[data[name]['sig'][col] for name in histNames]]
        else:
            G2frame.colSigs += [None]
            
    # compute and add weight fractions to table if varied
    for phase in Phases:
        var = str(Phases[phase]['pId'])+':*:Scale'
        if var not in combinedVaryList+list(depValDict.keys()): continue
        wtFrList = []
        sigwtFrList = []
        for i,name in enumerate(histNames):
            if name not in Phases[phase]['Histograms']:
                wtFrList.append(None)
                sigwtFrList.append(0.0)
                continue
            elif not Phases[phase]['Histograms'][name]['Use']:
                wtFrList.append(None)
                sigwtFrList.append(0.0)
                continue
            wtFrSum = 0.
            for phase1 in Phases:
                if name not in Phases[phase1]['Histograms']: continue
                if not Phases[phase1]['Histograms'][name]['Use']: continue
                wtFrSum += Phases[phase1]['Histograms'][name]['Scale'][0]*Phases[phase1]['General']['Mass']
            var = str(Phases[phase]['pId'])+':'+str(i)+':Scale'
            wtFr = Phases[phase]['Histograms'][name]['Scale'][0]*Phases[phase]['General']['Mass']/wtFrSum
            wtFrList.append(wtFr)
            if var in data[name]['varyList']:
                sig = data[name]['sig'][data[name]['varyList'].index(var)]*wtFr/Phases[phase]['Histograms'][name]['Scale'][0]
            else:
                sig = 0.0
            sigwtFrList.append(sig)
        colLabels.append(str(Phases[phase]['pId'])+':*:WgtFrac')
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
        G2frame.colList += [wtFrList]
        G2frame.colSigs += [sigwtFrList]
                
    # evaluate Pseudovars, their ESDs and add them to grid
    for expr in data['SeqPseudoVars']:
        obj = data['SeqPseudoVars'][expr]
        calcobj = G2obj.ExpressionCalcObj(obj)
        valList = []
        esdList = []
        for seqnum,name in enumerate(histNames):
            sigs = data[name]['sig']
            G2mv.InitVars()
            parmDict = data[name].get('parmDict')
            constraintInfo = data[name].get('constraintInfo',[[],[],{},[],seqnum])
            groups,parmlist,constrDict,fixedList,ihst = constraintInfo
            varyList = data[name]['varyList']
            parmDict = data[name]['parmDict']
            msg = G2mv.EvaluateMultipliers(constrDict,parmDict)
            if msg:
                print('Unable to interpret multiplier(s) for',name,':',msg)
                continue
            G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict,SeqHist=ihst)
            if 'Dist' in expr:
                derivs = G2mth.CalcDistDeriv(obj.distance_dict,obj.distance_atoms, parmDict)
                pId = obj.distance_dict['pId']
                aId,bId = obj.distance_atoms
                varyNames = ['%d::dA%s:%d'%(pId,ip,aId) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId) for ip in ['x','y','z']]
                VCoV = G2mth.getVCov(varyNames,varyList,data[name]['covMatrix'])
                esdList.append(np.sqrt(np.inner(derivs,np.inner(VCoV,derivs.T)) ))
#                GSASIIpath.IPyBreak()
            elif 'Angle' in expr:
                derivs = G2mth.CalcAngleDeriv(obj.angle_dict,obj.angle_atoms, parmDict)
                pId = obj.angle_dict['pId']
                aId,bId = obj.angle_atoms
                varyNames = ['%d::dA%s:%d'%(pId,ip,aId) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId[0]) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId[1]) for ip in ['x','y','z']]
                VCoV = G2mth.getVCov(varyNames,varyList,data[name]['covMatrix'])
                esdList.append(np.sqrt(np.inner(derivs,np.inner(VCoV,derivs.T)) ))
            else:
                derivs = np.array(
                    [EvalPSvarDeriv(calcobj,parmDict.copy(),sampleDict[name],var,ESD)
                     for var,ESD in zip(varyList,sigs)])
                if None in list(derivs):
                    esdList.append(None)
                else:
                    esdList.append(np.sqrt(
                        np.inner(derivs,np.inner(data[name]['covMatrix'],derivs.T)) ))
            PSvarDict = parmDict.copy()
            PSvarDict.update(sampleDict[name])
            UpdateParmDict(PSvarDict)
            calcobj.UpdateDict(PSvarDict)
            valList.append(calcobj.EvalExpression())
#            if calcobj.su is not None: esdList[-1] = calcobj.su
        if not esdList:
            esdList = None
        G2frame.colList += [valList]
        G2frame.colSigs += [esdList]
        colLabels += [expr]
        Types += [wg.GRID_VALUE_FLOAT+':10,5']
    #---- table build done -------------------------------------------------------------

    # Make dict needed for creating & editing pseudovars (PSvarDict).
    
    name = histNames[0]
    parmDict = data[name].get('parmDict',{})
    PSvarDict = parmDict.copy()
    PSvarDict.update(sampleParms)
    UpdateParmDict(PSvarDict)
    # Also dicts of variables 
    # for Parametric fitting from the data table
    parmDict = dict(zip(colLabels,list(zip(*G2frame.colList))[0])) # scratch dict w/all values in table
    parmDict.update({var:val for var,val in newCellDict.values()}) #  add varied reciprocal cell terms
    del parmDict['Use']
    name = histNames[0]

    # remove selected items from table
    saveColLabels = colLabels[:]
    if G2frame.SeqTblHideList is None:      #set default hides
        G2frame.SeqTblHideList = [item for item in saveColLabels if 'Back' in item]
        G2frame.SeqTblHideList += [item for item in saveColLabels if 'dA' in item]
        G2frame.SeqTblHideList += [item for item in saveColLabels if ':*:D' in item]
    #******************************************************************************
    # create a set of values for example evaluation of pseudovars and 
    # this does not work for refinements that have differing numbers of variables.
    VarDict = {}
    for i,var in enumerate(colLabels):
        if var in ['Use','Rwp',u'\u0394\u03C7\u00B2 (%)']: continue
        if G2frame.colList[i][0] is None:
            val,sig = firstValueDict.get(var,[None,None])
        elif G2frame.colSigs[i]:
            val,sig = G2frame.colList[i][0],G2frame.colSigs[i][0]
        else:
            val,sig = G2frame.colList[i][0],None
        if striphist(var) not in Dlookup:
            VarDict[var] = val
    # add recip cell coeff. values
    VarDict.update({var:val for var,val in newCellDict.values()})

    # remove items to be hidden from table
    for l in reversed(range(len(colLabels))):
        if colLabels[l] in G2frame.SeqTblHideList:
            del colLabels[l]
            del G2frame.colList[l]
            del G2frame.colSigs[l]

    # make a copy of the column labels substituting alternate labels when defined
    displayLabels = colLabels[:]
    for i,l in enumerate(colLabels):
        if l in variableLabels:
            displayLabels[i] = variableLabels[l]
            
    G2frame.dataWindow.ClearData()
    G2frame.dataWindow.currentGrids = []
    G2frame.dataDisplay = G2G.GSGrid(parent=G2frame.dataWindow)
    G2frame.dataDisplay.SetScrollRate(10,10)
    G2frame.dataWindow.GetSizer().Add(G2frame.dataDisplay,1,wx.ALL|wx.EXPAND)
    if histNames[0].startswith('PWDR'):
        #rowLabels = [str(i)+': '+l[5:30] for i,l in enumerate(histNames)]
        rowLabels = [l[5:] for i,l in enumerate(histNames)]
    else:
        rowLabels = histNames
    G2frame.SeqTable = G2G.Table([list(cl) for cl in zip(*G2frame.colList)],     # convert from columns to rows
        colLabels=displayLabels,rowLabels=rowLabels,types=Types)
    G2frame.dataDisplay.SetTable(G2frame.SeqTable, True)
    # make all Use editable all others ReadOnly
    for c in range(len(colLabels)):
        for r in range(nRows):
            if c == 1:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    if 'phoenix' in wx.version():
        G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
    else:
        G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
#    G2frame.dataDisplay.Bind(wx.EVT_KEY_UP,OnKeyUp)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, PlotSSelect)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, PlotSelect)
#    G2frame.dataDisplay.Bind(wg.EVT_GRID_SELECT_CELL,PlotSSelect)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, SetLabelString)
    G2frame.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    # highlight unconverged shifts 
    if histNames[0][:4] not in ['SASD','IMG ','REFD',]:
        for row,name in enumerate(histNames):
            deltaChi = G2frame.SeqTable.GetValue(row,deltaChiCol)
            if deltaChi > 10.:
                G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,0,0))
            elif deltaChi > 1.0:
                G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,255,0))
    G2frame.dataDisplay.InstallGridToolTip(GridSetToolTip,GridColLblToolTip)
    #G2frame.dataDisplay.SendSizeEvent() # resize needed on mac
    #G2frame.dataDisplay.Refresh() # shows colored text on mac
    G2frame.dataWindow.SetDataSize()

################################################################################
#####  Main PWDR panel
################################################################################           
       
def UpdatePWHKPlot(G2frame,kind,item):
    '''Called when the histogram main tree entry is called. Displays the
    histogram weight factor, refinement statistics for the histogram
    and the range of data for a simulation.

    Also invokes a plot of the histogram.
    '''
    def onEditSimRange(event):
        'Edit simulation range'
        if 'TOF' in data[0].get('simType','CW'):
            inp = [
                min(data[1][0])/1000.,
                max(data[1][0])/1000.,
                None
                ]
            inp[2] = (np.log(inp[1]) - np.log(inp[0]))/(len(data[1][0])-1.)
            names = ('start TOF(ms)', 'end TOF(ms)', 'DT/T')
            dlg = G2G.ScrolledMultiEditor(
                G2frame,[inp] * len(inp), range(len(inp)), names,
                header='Edit simulation range',
                minvals=(0.5,1.0,0.0001),
                maxvals=(200.,200.,.001),
                )            
        else:
            inp = [
                min(data[1][0]),
                max(data[1][0]),
                None
                ]
            inp[2] = (inp[1] - inp[0])/(len(data[1][0])-1.)
            names = ('start angle', 'end angle', 'step size')
            dlg = G2G.ScrolledMultiEditor(
                G2frame,[inp] * len(inp), range(len(inp)), names,
                header='Edit simulation range',
                minvals=(0.001,0.001,0.0001),
                maxvals=(180.,180.,.1),
                )
        dlg.CenterOnParent()
        val = dlg.ShowModal()
        dlg.Destroy()
        if val != wx.ID_OK: return
        if inp[0] > inp[1]:
            end,start,step = inp
        else:                
            start,end,step = inp
        step = abs(step)
        if 'TOF' in data[0].get('simType','CW'):
            N = (np.log(end)-np.log(start))/step
            newdata = np.exp((np.arange(0,N))*step+np.log(start*1000.))
        else:
            N = int((end-start)/step)+1
            newdata = np.linspace(start,end,N,True)
            if len(newdata) < 2: return # too small a range - reject
        data[1] = [newdata,np.zeros_like(newdata),np.ones_like(newdata),
            np.zeros_like(newdata),np.zeros_like(newdata),np.zeros_like(newdata)]
        Tmin = newdata[0]
        Tmax = newdata[-1]
        G2frame.GPXtree.SetItemPyData(GetGPXtreeItemId(G2frame,item,'Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        wx.CallAfter(UpdatePWHKPlot,G2frame,kind,item) # redisplay data screen
        
    def OnPlot1DHKL(event):
        refList = data[1]['RefList']
        G2plt.Plot1DSngl(G2frame,newPlot=True,hklRef=refList,Super=Super,Title=phaseName)

    def OnPlot3DHKL(event):
        refList = data[1]['RefList']
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        Vpoint = np.array([int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))])
        controls = {'Type' : 'Fosq','Iscale' : False,'HKLmax' : Hmax,'HKLmin' : Hmin,'Zone':False,'viewKey':'L',
            'FoMax' : FoMax,'Scale' : 1.0,'Drawing':{'viewPoint':[Vpoint,[]],'default':Vpoint[:],
            'backColor':[0,0,0],'depthFog':False,'Zclip':10.0,'cameraPos':10.,'Zstep':0.05,'viewUp':[0,1,0],
            'Scale':1.0,'oldxy':[],'viewDir':[0,0,1]},'Super':Super,'SuperVec':SuperVec}
        G2plt.Plot3DSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
        
    def OnPlotAll3DHKL(event):
        choices = GetGPXtreeDataNames(G2frame,['HKLF',])
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select reflection sets to plot',
            'Use data',choices)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                refNames = [choices[i] for i in dlg.GetSelections()]
            else:
                return
        finally:
            dlg.Destroy()
        refList = np.zeros(0)
        for name in refNames:
            Id = GetGPXtreeItemId(G2frame,G2frame.root, name)
            reflData = G2frame.GPXtree.GetItemPyData(Id)[1]
            if len(refList):
                refList = np.concatenate((refList,reflData['RefList']))
            else:
                refList = reflData['RefList']
            
        FoMax = np.max(refList.T[8+Super])
        Hmin = np.array([int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))])
        Hmax = np.array([int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))])
        Vpoint = [int(np.mean(refList.T[0])),int(np.mean(refList.T[1])),int(np.mean(refList.T[2]))]
        controls = {'Type' : 'Fosq','Iscale' : False,'HKLmax' : Hmax,'HKLmin' : Hmin,'Zone':False,'viewKey':'L',
            'FoMax' : FoMax,'Scale' : 1.0,'Drawing':{'viewPoint':[Vpoint,[]],'default':Vpoint[:],
            'backColor':[0,0,0],'depthFog':False,'Zclip':10.0,'cameraPos':10.,'Zstep':0.05,'viewUp':[0,1,0],
            'Scale':1.0,'oldxy':[],'viewDir':[1,0,0]},'Super':Super,'SuperVec':SuperVec}
        G2plt.Plot3DSngl(G2frame,newPlot=True,Data=controls,hklRef=refList,Title=phaseName)
                  
    def OnMergeHKL(event):
        Name = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        Inst = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,
            G2frame.PatternId,'Instrument Parameters'))
        CId = GetGPXtreeItemId(G2frame,G2frame.PatternId,'Comments')
        if CId:
            Comments = G2frame.GPXtree.GetItemPyData(CId)
        else:
            Comments = []
        refList = np.copy(data[1]['RefList'])
        Comments.append(' Merging %d reflections from %s'%(len(refList),Name))
        dlg = MergeDialog(G2frame,data)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                Trans,Cent,Laue = dlg.GetSelection()
            else:
                return
        finally:
            dlg.Destroy()
        Super = data[1]['Super']
        isup = 0
        if Super:
            isup = 1
        refList,badRefs = G2lat.transposeHKLF(Trans,isup,refList)
        if len(badRefs):    #do I want to list badRefs?
            G2frame.ErrorDialog('Failed transformation','Matrix yields fractional hkl indices')
            return
        Comments.append(" Transformation M*H = H' applied; M=")
        Comments.append(str(Trans))
        refList = G2lat.LaueUnique(Laue,refList)
        dlg = wx.ProgressDialog('Build HKL dictonary','',len(refList)+1, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        HKLdict = {}
        for ih,hkl in enumerate(refList):
            if str(hkl[:3+Super]) not in HKLdict:
                HKLdict[str(hkl[:3+Super])] = [hkl[:3+Super],[hkl[3+Super:],]]
            else:
                HKLdict[str(hkl[:3+Super])][1].append(hkl[3+Super:])
            dlg.Update(ih)
        dlg.Destroy()
        mergeRef = []
        dlg = wx.ProgressDialog('Processing merge','',len(HKLdict)+1, 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
        sumDf = 0.
        sumFo = 0.
        for ih,hkl in enumerate(HKLdict):
            HKL = HKLdict[hkl]
            newHKL = list(HKL[0])+list(HKL[1][0])
            if len(HKL[1]) > 1:
                fos = np.array(HKL[1])
                wFo = 1/fos[:,3]**2
                Fo = np.average(fos[:,2],weights=wFo)
                std = np.std(fos[:,2])
                sig = np.sqrt(np.mean(fos[:,3])**2+std**2)
                sumFo += np.sum(fos[:,2])
                sumDf += np.sum(np.abs(fos[:,2]-Fo))
                dlg.Update(ih)
                newHKL[5+Super] = Fo
                newHKL[6+Super] = sig
                newHKL[8+Super] = Fo
            if newHKL[5+Super] > 0.:
                mergeRef.append(list(newHKL)) 
        dlg.Destroy()
        if Super:
            mergeRef = G2mth.sortArray(G2mth.sortArray(G2mth.sortArray(G2mth.sortArray(mergeRef,3),2),1),0)
        else:
            mergeRef = G2mth.sortArray(G2mth.sortArray(G2mth.sortArray(mergeRef,2),1),0)
        mergeRef = np.array(mergeRef)
        if sumFo:
            mtext = ' merge R = %6.2f%s for %d reflections in %s'%(100.*sumDf/sumFo,'%',mergeRef.shape[0],Laue)
            print (mtext)
            Comments.append(mtext)
        else:
            print( 'nothing to merge for %s reflections'%(mergeRef.shape[0]))
        HKLFlist = []
        newName = Name+u' '+Laue
        if G2frame.GPXtree.GetCount():
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.GPXtree.GetItemText(item)
                if name.startswith('HKLF ') and name not in HKLFlist:
                    HKLFlist.append(name)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        newName = G2obj.MakeUniqueLabel(newName,HKLFlist)
        newData = copy.deepcopy(data)
        newData[0]['ranId'] = ran.randint(0,sys.maxsize)
        newData[1]['RefList'] = mergeRef
        Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=newName)
        G2frame.GPXtree.SetItemPyData(
            G2frame.GPXtree.AppendItem(Id,text='Comments'),Comments)
        G2frame.GPXtree.SetItemPyData(Id,newData)
        G2frame.GPXtree.SetItemPyData(
            G2frame.GPXtree.AppendItem(Id,text='Instrument Parameters'),Inst)
        G2frame.GPXtree.SetItemPyData(
            G2frame.GPXtree.AppendItem(Id,text='Reflection List'),{})  #dummy entry for GUI use
                   
    def OnErrorAnalysis(event):
        G2plt.PlotDeltSig(G2frame,kind)
        
#    def OnCompression(event):
#        data[0] = int(comp.GetValue())
        
    def onCopyPlotCtrls(event):
        '''Respond to menu item to copy multiple sections from a histogram.
        Need this here to pass on the G2frame object. 
        '''
        G2pdG.CopyPlotCtrls(G2frame)

    def onCopySelectedItems(event):
        '''Respond to menu item to copy multiple sections from a histogram.
        Need this here to pass on the G2frame object. 
        '''
        G2pdG.CopySelectedHistItems(G2frame)

    def OnAddMag(event):
        'Respond to the request to create a magnification region'
        if not data[0]['Magnification']: data[0]['Magnification'] = [[None,1.0]]
        data[0]['Magnification'] += [[(data[1][0][0]+data[1][0][-1])/2,2.]]
        wx.CallAfter(UpdatePWHKPlot,G2frame,kind,G2frame.PatternId)
    def OnDelMag(event):
        'Respond to the request to delete a magnification region'
        del data[0]['Magnification'][event.EventObject.row]
        if len(data[0]['Magnification']) == 1: data[0]['Magnification'] = []
        wx.CallAfter(UpdatePWHKPlot,G2frame,kind,G2frame.PatternId)
    def OnEditMag(**args):
        'Update to show edits to mag factors in window and plot'
        wx.CallAfter(UpdatePWHKPlot,G2frame,kind,G2frame.PatternId)

    # Start of UpdatePWHKPlot
    data = G2frame.GPXtree.GetItemPyData(item)
#patches
    if not data:
        return
    if 'wtFactor' not in data[0]:
        data[0] = {'wtFactor':1.0}
#    if kind == 'PWDR' and 'Compression' not in data[0]:
#        data[0]['Compression'] = 1
    #if isinstance(data[1],list) and kind == 'HKLF':
    if 'list' in str(type(data[1])) and kind == 'HKLF':
        RefData = {'RefList':[],'FF':[]}
        for ref in data[1]:
            RefData['RefList'].append(ref[:11]+[ref[13],])
            RefData['FF'].append(ref[14])
        data[1] = RefData
        G2frame.GPXtree.SetItemPyData(item,data)
#end patches
    if kind in ['PWDR','SASD','REFD']:
        SetDataMenuBar(G2frame,G2frame.dataWindow.PWDRMenu)
        G2frame.Bind(wx.EVT_MENU, OnErrorAnalysis, id=G2G.wxID_PWDANALYSIS)
        G2frame.Bind(wx.EVT_MENU, onCopySelectedItems, id=G2G.wxID_PWDCOPY)
        G2frame.Bind(wx.EVT_MENU, onCopyPlotCtrls, id=G2G.wxID_PLOTCTRLCOPY)
    elif kind in ['HKLF',]:
        SetDataMenuBar(G2frame,G2frame.dataWindow.HKLFMenu)
        G2frame.Bind(wx.EVT_MENU, OnErrorAnalysis, id=G2G.wxID_PWDANALYSIS)
        G2frame.Bind(wx.EVT_MENU, OnMergeHKL, id=G2G.wxID_MERGEHKL)
        G2frame.Bind(wx.EVT_MENU, OnPlot1DHKL, id=G2G.wxID_1DHKLSTICKPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlot3DHKL, id=G2G.wxID_PWD3DHKLPLOT)
        G2frame.Bind(wx.EVT_MENU, OnPlotAll3DHKL, id=G2G.wxID_3DALLHKLPLOT)
    
    if G2frame.dataWindow:
        G2frame.dataWindow.ClearData() 
    mainSizer = G2frame.dataWindow.GetSizer()
    mainSizer.Add((5,5),)
    wtSizer = wx.BoxSizer(wx.HORIZONTAL)
    wtSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Weight factor: '),0,WACV)
    wtSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[0],'wtFactor',nDig=(10,3),xmin=1.e-9),0,WACV)
#    if kind == 'PWDR':         #possible future compression feature; NB above patch as well
#        wtSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Compression factor: '),0,WACV)
#        choice = ['1','2','3','4','5','6']
#        comp = wx.ComboBox(parent=G2frame.dataWindow,choices=choice,
#            style=wx.CB_READONLY|wx.CB_DROPDOWN)
#        comp.SetValue(str(data[0]['Compression']))
#        comp.Bind(wx.EVT_COMBOBOX, OnCompression)
#        wtSizer.Add(comp,0,WACV)
    mainSizer.Add(wtSizer)
    wtSizer = wx.BoxSizer(wx.HORIZONTAL)
    wtSizer.Add(wx.StaticText(G2frame.dataWindow,-1,' Histogram label: '),0,WACV)
    if 'histTitle' not in data[0]: data[0]['histTitle'] = ''
    wtSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,data[0],'histTitle',typeHint=str,
        notBlank=False,size=(300,-1)),1,WACV)
    mainSizer.Add(wtSizer,0)
    if data[0].get('Dummy',False):
        simSizer = wx.BoxSizer(wx.HORIZONTAL)
        Tmin = min(data[1][0])
        Tmax = max(data[1][0])
        num = len(data[1][0])
        if 'TOF' in data[0].get('simType','CW'):
            step = (np.log(Tmax) - np.log(Tmin))/(num-1.)
            t = u'\u00b5s'
            lbl =  u'Simulation range: {:.2f} to {:.2f} {:s} with {:.4f} resolution ({:d} points)'
        else:
            step = (Tmax - Tmin)/(num-1)
            t = u'2\u03b8' # 2theta
            lbl =  u'Simulation range: {:.2f} to {:.2f} {:s} with {:.4f} steps ({:d} points)'
        lbl += u'\n(Edit range resets observed intensities).'
        lbl = lbl.format(Tmin,Tmax,t,step,num)
        simSizer.Add(wx.StaticText(G2frame.dataWindow,wx.ID_ANY,lbl),0,WACV)
        but = wx.Button(G2frame.dataWindow,wx.ID_ANY,"Edit range")
        but.Bind(wx.EVT_BUTTON,onEditSimRange)
        mainSizer.Add(simSizer)
        mainSizer.Add(but,0,WACV)
    if 'Nobs' in data[0]:
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,
            ' Data residual wR: %.3f%% on %d observations'%(data[0]['wR'],data[0]['Nobs'])))
        if kind == 'PWDR':
            mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,
                ' Durbin-Watson statistic: %.3f'%(data[0].get('Durbin-Watson',0.))))
        for value in data[0]:
            if 'Nref' in value:
                pfx = value.split('Nref')[0]
                name = data[0].get(pfx.split(':')[0]+'::Name','?')
                if 'SS' in value:
                    mainSizer.Add((5,5),)
                    mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' For incommensurate phase '+name+u':'))
                    for m,(Rf2,Rf,Nobs) in enumerate(zip(data[0][pfx+'Rf^2'],data[0][pfx+'Rf'],data[0][value])):
                        mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,
                            u' m = +/- %d: RF\u00b2: %.3f%%, RF: %.3f%% on %d reflections  '% \
                            (m,Rf2,Rf,Nobs)))
                else:
                    mainSizer.Add((5,5),)
                    mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,u' For phase '+name+u':'))
                    mainSizer.Add(wx.StaticText(G2frame.dataWindow,-1,
                        u' Unweighted phase residuals RF\u00b2: %.3f%%, RF: %.3f%% on %d reflections  '% \
                        (data[0][pfx+'Rf^2'],data[0][pfx+'Rf'],data[0][value])))

    # Draw edit box for Magnification factors/positions
    if kind == 'PWDR':
        if 'Magnification' not in data[0]:
            data[0]['Magnification'] = []
        mainSizer.Add((-1,10))
        lenmag = len(data[0]['Magnification'])
        data[0]['Magnification'][1:] = sorted(data[0]['Magnification'][1:],key=lambda x: x[0])
        if lenmag > 1:
            panel = wx.StaticBox(G2frame.dataWindow, wx.ID_ANY, 'Magnification regions',
                                 style=wx.ALIGN_CENTER)
            mSizer = wx.StaticBoxSizer(panel,wx.VERTICAL)
            magSizer = wx.FlexGridSizer(lenmag+1,3,0,0)
            Name = G2frame.GPXtree.GetItemText(G2frame.PatternId)
            Inst = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,
                G2frame.PatternId,'Instrument Parameters'))
            if 'C' in Inst[0]['Type'][0]:
                magSizer.Add(wx.StaticText(panel,-1,'2Theta'),1,wx.ALIGN_CENTER,1)
            else:
                magSizer.Add(wx.StaticText(panel,-1,'TOF'),1,wx.ALIGN_CENTER,1)
            magSizer.Add(wx.StaticText(panel,-1,'Magnification\nfactor',
                style=wx.ALIGN_CENTRE_HORIZONTAL),1,wx.ALIGN_CENTER,1)
            magSizer.Add(wx.StaticText(panel,-1,'Delete\nbutton',
                style=wx.ALIGN_CENTRE_HORIZONTAL),1,wx.ALIGN_CENTER,1)
            magSizer.Add(wx.StaticText(panel,-1,'(start)'),1,wx.ALIGN_CENTER,1)
            edit = G2G.ValidatedTxtCtrl(panel,data[0]['Magnification'][0],1,
                nDig=(7,2),xmin=0.01,xmax=1000.,OnLeave=OnEditMag,size=(65,-1))
            magSizer.Add(edit,1,wx.ALIGN_CENTER,5)
            magSizer.Add((1,1))
            for i in range(1,lenmag):
                edit = G2G.ValidatedTxtCtrl(panel,data[0]['Magnification'][i],0,nDig=(10,3), 
                    xmin=data[1][0][0],xmax=data[1][0][-1],OnLeave=OnEditMag)
                magSizer.Add(edit)
                edit = G2G.ValidatedTxtCtrl(panel,data[0]['Magnification'][i],1,
                    nDig=(7,2),xmin=0.01,xmax=1000.,OnLeave=OnEditMag,size=(65,-1))
                magSizer.Add(edit,1,wx.ALIGN_CENTER,5)
                delmag = wx.Button(panel,wx.ID_ANY,label='Del',size=(40,-1))
                delmag.Bind(wx.EVT_BUTTON,OnDelMag)
                delmag.row = i
                magSizer.Add(delmag,1,wx.ALIGN_CENTER,5)
            mSizer.Add(magSizer)
        else:
            panel = G2frame.dataWindow
            mSizer = wx.BoxSizer(wx.VERTICAL)
        addmag = wx.Button(panel,wx.ID_ANY,label='Add a magnification region')
        addmag.Bind(wx.EVT_BUTTON,OnAddMag)
        mSizer.Add(addmag,1,wx.ALIGN_CENTER,1)
        mainSizer.Add(mSizer)
        
    G2frame.GPXtree.SetItemPyData(item,data)
    G2frame.PatternId = item
    if kind in ['PWDR','SASD','REFD',]:
        NewPlot = True
        if 'xylim' in dir(G2frame):
            NewPlot = False
        G2plt.PlotPatterns(G2frame,plotType=kind,newPlot=NewPlot)
    elif kind == 'HKLF':
        Name = G2frame.GPXtree.GetItemText(item)
        phaseName = G2pdG.IsHistogramInAnyPhase(G2frame,Name)
        if phaseName:
            pId = GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
            phaseId =  GetGPXtreeItemId(G2frame,pId,phaseName)
            General = G2frame.GPXtree.GetItemPyData(phaseId)['General']
            Super = General.get('Super',0)
            SuperVec = General.get('SuperVec',[])
        else:
            Super = 0
            SuperVec = []       
        refList = data[1]['RefList']
        FoMax = np.max(refList.T[5+data[1].get('Super',0)])
        page = G2frame.G2plotNB.nb.GetSelection()
        tab = ''
        if page >= 0:
            tab = G2frame.G2plotNB.nb.GetPageText(page)
        if '3D' in tab:
            Page = G2frame.G2plotNB.nb.GetPage(page)
            controls = Page.controls
            G2plt.Plot3DSngl(G2frame,newPlot=False,Data=controls,hklRef=refList,Title=phaseName)
        else:
            controls = {'Type' : 'Fo','ifFc' : True,     
                'HKLmax' : [int(np.max(refList.T[0])),int(np.max(refList.T[1])),int(np.max(refList.T[2]))],
                'HKLmin' : [int(np.min(refList.T[0])),int(np.min(refList.T[1])),int(np.min(refList.T[2]))],
                'FoMax' : FoMax,'Zone' : '001','Layer' : 0,'Scale' : 1.0,'Super':Super,'SuperVec':SuperVec}
            G2plt.PlotSngl(G2frame,newPlot=True,Data=controls,hklRef=refList)
    G2frame.dataWindow.SetDataSize()
                 
################################################################################
#####  Data (GPX) tree routines
################################################################################           
       
def GetGPXtreeDataNames(G2frame,dataTypes):
    '''Finds all items in tree that match a 4 character prefix
    
    :param wx.Frame G2frame: Data tree frame object
    :param list dataTypes: Contains one or more data tree item types to be matched
      such as ['IMG '] or ['PWDR','HKLF']
    :returns: a list of tree item names for the matching items  
    '''
    names = []
    item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)        
    while item:
        name = G2frame.GPXtree.GetItemText(item)
        if name[:4] in dataTypes:
            names.append(name)
        item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
    return names
                          
def GetGPXtreeItemId(G2frame, parentId, itemText):
    '''Find the tree item that matches the text in itemText starting with parentId

    :param wx.Frame G2frame: Data tree frame object
    :param wx.TreeItemId parentId: tree item to start search with
    :param str itemText: text for tree item
    '''
    try:
        item, cookie = G2frame.GPXtree.GetFirstChild(parentId)
        while item:
            if G2frame.GPXtree.GetItemText(item) == itemText:
                return item
            item, cookie = G2frame.GPXtree.GetNextChild(parentId, cookie)
        return 0
    except:         #trap C++ error?
        return 0

def SelectDataTreeItem(G2frame,item,oldFocus=None):
    '''Called from :meth:`GSASIIdataGUI.GSASII.OnDataTreeSelChanged` when a item is selected on the tree.
    Also called from GSASII.OnGPXtreeEndDrag, OnAddPhase -- might be better to select item, triggering
    the the bind to SelectDataTreeItem

    Also Called in GSASIIphsGUI.UpdatePhaseData by OnTransform callback. 
    '''
    if G2frame.PickIdText == G2frame.GetTreeItemsList(item): # don't redo the current data tree item 
        if GSASIIpath.GetConfigValue('debug'): print('Skipping SelectDataTreeItem as G2frame.PickIdText unchanged')
        return

    # save or finish processing of outstanding events
    for grid in G2frame.dataWindow.currentGrids:  # complete any open wx.Grid edits
        #if GSASIIpath.GetConfigValue('debug'): print 'Testing grid edit in',grid
        try: 
            if grid.IsCellEditControlEnabled(): # complete any grid edits in progress
                if GSASIIpath.GetConfigValue('debug'): print ('Completing grid edit in%s'%str(grid))
                grid.SaveEditControlValue()
                # not sure if the next two things do anything
                grid.HideCellEditControl()
                grid.DisableCellEditControl()
        except:
            pass
    if G2frame.dataWindow.GetLabel() == 'Comments': # save any recently entered comments 
        try:
            data = [G2frame.dataDisplay.GetValue()]
            G2frame.dataDisplay.ClearData() 
            Id = GetGPXtreeItemId(G2frame,G2frame.root, 'Comments')
            if Id: G2frame.GPXtree.SetItemPyData(Id,data)
        except:     #clumsy but avoids dead window problem when opening another project
            pass
    elif G2frame.dataWindow.GetLabel() == 'Notebook': # save any recent notebook entries
        try:
            data = [G2frame.dataDisplay.GetValue()]
            G2frame.dataDisplay.ClearData() 
            Id = GetGPXtreeItemId(G2frame,G2frame.root, 'Notebook')
            if Id: G2frame.GPXtree.SetItemPyData(Id,data)
        except:     #clumsy but avoids dead window problem when opening another project
            pass
#    elif 'Phase Data for' in G2frame.dataWindow.GetLabel():
#        if G2frame.dataDisplay: 
#            oldPage = G2frame.dataDisplay.GetSelection()
        
    SetDataMenuBar(G2frame)
    G2frame.SetTitleByGPX()
    G2frame.PickId = item
    G2frame.PickIdText = None
    parentID = G2frame.root
    if item == G2frame.root:
        G2frame.dataWindow.ClearData()
        G2frame.helpKey = "Data tree"
        G2frame.dataWindow.GetSizer().Add(
            wx.StaticText(G2frame.dataWindow, wx.ID_ANY,
                    'Select an item from the tree to see/edit parameters'))
        G2frame.dataWindow.SetDataSize()
        return
    else:
        # Set up the help entry for the current selected tree item
        try:    #don't know why we get here when opening new project
            parentID = G2frame.GPXtree.GetItemParent(item)
            # save name of calling tree item for help. N.B. may want to override this later
            prfx = G2frame.GPXtree.GetItemText(item)
            if prfx:
                prfx = prfx.split()[0].upper()
            else:   #just deleted item - escape!!
                return
            prfx1 = G2frame.GPXtree.GetItemText(parentID).split()[0]
            if prfx in ('IMG','PKS','PWDR','SASD','HKLF','PDF','REFD',):
                G2frame.dataWindow.helpKey = prfx
            elif prfx1 == 'Phases':
                G2frame.dataWindow.helpKey = 'Phases'
            elif prfx1 in ('IMG','PKS','PWDR','SASD','HKLF','PDF','REFD',):
                suffix = G2frame.GPXtree.GetItemText(item)
                suffix1 = suffix.split()[0]
                if '(Q)' in suffix1 or '(R)' in suffix1: suffix = suffix1
                G2frame.dataWindow.helpKey = prfx1 + '_' + suffix
            else:
                G2frame.dataWindow.helpKey = G2frame.GPXtree.GetItemText(item) # save name of calling tree item for help
        except IndexError:
            G2frame.dataWindow.helpKey = ''
            if GSASIIpath.GetConfigValue('debug'):
                print ('bug: why here? prfx=%s prfx1=%s'%(prfx,prfx1))
                G2obj.HowDidIgetHere()
        
    # clear out the old panel contents
    if G2frame.dataWindow:
        G2frame.dataWindow.ClearData() 
#    G2frame.dataWindow.ClearData()
    # process first-level entries in tree
    if G2frame.GPXtree.GetItemParent(item) == G2frame.root:
        G2frame.PatternId = 0
        if G2frame.GPXtree.GetItemText(item) == 'Notebook':
            SetDataMenuBar(G2frame,G2frame.dataWindow.DataNotebookMenu)
            data = G2frame.GPXtree.GetItemPyData(item)
            UpdateNotebook(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item) == 'Controls':
            data = G2frame.GPXtree.GetItemPyData(item)
            if not data:           #fill in defaults
                data = copy.copy(G2obj.DefaultControls)    #least squares controls
                G2frame.GPXtree.SetItemPyData(item,data)                             
            for i in G2frame.Refine: i.Enable(True)
            UpdateControls(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item).startswith('Sequential '):
            G2frame.dataWindow.helpKey = 'Sequential'  # for now all sequential refinements are documented in one place
            data = G2frame.GPXtree.GetItemPyData(item)
            UpdateSeqResults(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item) == 'Covariance':
            data = G2frame.GPXtree.GetItemPyData(item)
            text = ''
            if 'Rvals' in data:
                Nvars = len(data['varyList'])
                Rvals = data['Rvals']
                text = ('Residuals after last refinement:\n'+
                        '\twR = {:.3f}\n\tchi**2 = {:.1f}\n\tGOF = {:.2f}').format(
                        Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
                text += '\n\tNobs = {}\n\tNvals = {}\n\tSVD zeros = {}'.format(
                    Rvals['Nobs'],Nvars,Rvals.get('SVD0',0.))
                text += '\n\tmax shift/esd = {:.3f}'.format(Rvals.get('Max shft/sig',0.0))
                if 'lamMax' in Rvals:
                    text += '\n\tlog10 MaxLambda = {:.1f}'.format(np.log10(Rvals['lamMax']))
                G2frame.dataWindow.GetSizer().Add(
                    wx.StaticText(G2frame.dataWindow,wx.ID_ANY,text)
                )
            G2plt.PlotCovariance(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item) == 'Constraints':
            data = G2frame.GPXtree.GetItemPyData(item)
            #import imp
            #imp.reload(G2cnstG)  # for testing changes to GSASIIconstrGUI
            G2cnstG.UpdateConstraints(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item) == 'Rigid bodies':
            data = G2frame.GPXtree.GetItemPyData(item)
            # import imp
            # imp.reload(G2cnstG)  # for testing changes to GSASIIconstrGUI
            # print('debug: reloaded',G2cnstG)
            G2cnstG.UpdateRigidBodies(G2frame,data)
        elif G2frame.GPXtree.GetItemText(item).startswith('IMG '):
            G2frame.Image = item
            data = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(
                G2frame,item,'Image Controls'))
            G2imG.UpdateImageData(G2frame,data)
            G2plt.PlotImage(G2frame,newPlot=True)
        elif G2frame.GPXtree.GetItemText(item).startswith('PKS '):
            G2frame.PatternId = item
            G2plt.PlotPowderLines(G2frame)
        elif G2frame.GPXtree.GetItemText(item).startswith('PWDR '):
            G2frame.PatternId = item
            #for i in G2frame.ExportPattern: i.Enable(True)
            if G2frame.EnablePlot:
                UpdatePWHKPlot(G2frame,'PWDR',item)
        elif G2frame.GPXtree.GetItemText(item).startswith('SASD '):
            G2frame.PatternId = item
            #for i in G2frame.ExportPattern: i.Enable(True)
            if G2frame.EnablePlot:
                UpdatePWHKPlot(G2frame,'SASD',item)
        elif G2frame.GPXtree.GetItemText(item).startswith('REFD '):
            G2frame.PatternId = item
            #for i in G2frame.ExportPattern: i.Enable(True)
            if G2frame.EnablePlot:
                UpdatePWHKPlot(G2frame,'REFD',item)
        elif G2frame.GPXtree.GetItemText(item).startswith('HKLF '):
            G2frame.Sngl = True
            UpdatePWHKPlot(G2frame,'HKLF',item)
        elif G2frame.GPXtree.GetItemText(item).startswith('PDF '):
            G2frame.PatternId = item
            for i in G2frame.ExportPDF: i.Enable(True) # this should be done on .gpx load; is done on OnMakePDFs (GSASII.py)
            data = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,item,'PDF Controls'))
            G2pdG.UpdatePDFGrid(G2frame,data)
            if len(data['G(R)']):
                G2plt.PlotISFG(G2frame,data,plotType='G(R)')
        elif G2frame.GPXtree.GetItemText(item) == 'Phases':
            G2frame.dataWindow.GetSizer().Add(
                wx.StaticText(G2frame.dataWindow,wx.ID_ANY,'Select one phase to see its parameters'))
        elif G2frame.GPXtree.GetItemText(item) == 'Restraints':
            data = G2frame.GPXtree.GetItemPyData(item)
#patch - put phases in restraint tree
            names = G2frame.GetPhaseNames()
            for name in names:
                if not GetGPXtreeItemId(G2frame,item,name):
                    G2frame.GPXtree.AppendItem(parent=item,text=name)
                if name not in data:
                    data[name] = {}
#end patch            
            G2frame.dataWindow.GetSizer().Add(
                wx.StaticText(G2frame.dataWindow,wx.ID_ANY,'Select one phase to see its restraints'))
    ############################################################################
    # process second-level entries in tree            
    elif G2frame.GPXtree.GetItemText(item) == 'PDF Peaks':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        peaks = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,G2frame.PatternId,'PDF Peaks'))
        data = G2frame.GPXtree.GetItemPyData(GetGPXtreeItemId(G2frame,G2frame.PatternId,'PDF Controls'))
        G2pdG.UpdatePDFPeaks(G2frame,peaks,data)
        if len(data['G(R)']):
            G2plt.PlotISFG(G2frame,data,plotType='G(R)',newPlot=True,peaks=peaks)            
    elif G2frame.GPXtree.GetItemText(item) == 'PDF Controls':
        for i in G2frame.ExportPDF: i.Enable(True) # this should be done on .gpx load; is done on OnMakePDFs (GSASII.py)
        G2frame.dataWindow.helpKey = G2frame.GPXtree.GetItemText(item) # special treatment to avoid PDF_PDF Controls
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(G2frame,data)
        if len(data['G(R)']):
            if 'I(Q)' in data:  G2plt.PlotISFG(G2frame,data,plotType='I(Q)')
            if 'S(Q)' in data:  G2plt.PlotISFG(G2frame,data,plotType='S(Q)')
            if 'F(Q)' in data:  G2plt.PlotISFG(G2frame,data,plotType='F(Q)')
            G2plt.PlotISFG(G2frame,data,plotType='G(R)')
    elif G2frame.GPXtree.GetItemText(parentID) == 'Phases':
        data = G2frame.GPXtree.GetItemPyData(item)
        # debug stuff
        # if GSASIIpath.GetConfigValue('debug'):
        #     print('Debug: reloading G2phG')
        #     import imp
        #     imp.reload(G2phG)
        # end debug stuff            
        G2phG.UpdatePhaseData(G2frame,item,data)
    elif G2frame.GPXtree.GetItemText(parentID) == 'Restraints':
        data = G2frame.GPXtree.GetItemPyData(parentID)
        phaseName = G2frame.GPXtree.GetItemText(item)
        if phaseName not in data:
            data[phaseName] = {}
        G2restG.UpdateRestraints(G2frame,data,phaseName)
    elif G2frame.GPXtree.GetItemText(item) == 'Comments':
        SetDataMenuBar(G2frame,G2frame.dataWindow.DataCommentsMenu)
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        UpdateComments(G2frame,data)
    elif G2frame.GPXtree.GetItemText(item) == 'Image Controls':
        G2frame.Image = G2frame.GPXtree.GetItemParent(item)
        masks = G2frame.GPXtree.GetItemPyData(
            GetGPXtreeItemId(G2frame,G2frame.Image, 'Masks'))
        data = G2frame.GPXtree.GetItemPyData(item)
        G2frame.ImageZ = G2imG.GetImageZ(G2frame,data)
        G2imG.UpdateImageControls(G2frame,data,masks)
        G2plt.PlotImage(G2frame,newPlot=False)
    elif G2frame.GPXtree.GetItemText(item) == 'Masks':
        G2frame.Image = G2frame.GPXtree.GetItemParent(item)
        masks = G2frame.GPXtree.GetItemPyData(item)
        data = G2frame.GPXtree.GetItemPyData(
            GetGPXtreeItemId(G2frame,G2frame.Image, 'Image Controls'))
        G2frame.ImageZ = G2imG.GetImageZ(G2frame,data)
        G2imG.UpdateMasks(G2frame,masks)
        G2plt.PlotImage(G2frame,newPlot=False)
    elif G2frame.GPXtree.GetItemText(item) == 'Stress/Strain':
        G2frame.Image = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(
            GetGPXtreeItemId(G2frame,G2frame.Image, 'Image Controls'))
        G2frame.ImageZ = G2imG.GetImageZ(G2frame,data,newRange=False)
        strsta = G2frame.GPXtree.GetItemPyData(item)
        G2plt.PlotStrain(G2frame,strsta,newPlot=True)
        G2plt.PlotImage(G2frame,newPlot=False)
        G2imG.UpdateStressStrain(G2frame,strsta)
    elif G2frame.GPXtree.GetItemText(item) == 'Peak List':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        data = G2frame.GPXtree.GetItemPyData(item)
#patch
        if 'list' in str(type(data)):
            data = {'peaks':data,'sigDict':{}}
            G2frame.GPXtree.SetItemPyData(item,data)
#end patch
        G2pdG.UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Background':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        G2pdG.UpdateBackground(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Limits':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        datatype = G2frame.GPXtree.GetItemText(G2frame.PatternId)[:4]
        data = G2frame.GPXtree.GetItemPyData(item)
        G2pdG.UpdateLimitsGrid(G2frame,data,datatype)
        G2plt.PlotPatterns(G2frame,plotType=datatype,newPlot=True)
    elif G2frame.GPXtree.GetItemText(item) == 'Instrument Parameters':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)[0]
        G2pdG.UpdateInstrumentGrid(G2frame,data)
        if 'P' in data['Type'][0]:          #powder data only
            G2plt.PlotPeakWidths(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Models':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        if prfx1 == 'SASD':
            G2pdG.UpdateModelsGrid(G2frame,data)
        elif prfx1 == 'REFD':
            G2pdG.UpdateREFDModelsGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType=prfx1)
        if prfx1 == 'SASD':
            if len(data['Size']['Distribution']):
                G2plt.PlotSASDSizeDist(G2frame)
            if len(data['Pair']['Distribution']):
                G2plt.PlotSASDPairDist(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Substances':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        G2pdG.UpdateSubstanceGrid(G2frame,data)
    elif G2frame.GPXtree.GetItemText(item) == 'Sample Parameters':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        datatype = G2frame.GPXtree.GetItemText(G2frame.PatternId)[:4]

        if 'Temperature' not in data:           #temp fix for old gpx files
            data = {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],'DisplaceX':[0.0,False],
                'DisplaceY':[0.0,False],'Temperature':300.,'Pressure':1.0,
                    'FreePrm1':0.,'FreePrm2':0.,'FreePrm3':0.,
                    'Gonio. radius':200.0}
            G2frame.GPXtree.SetItemPyData(item,data)
    
        G2pdG.UpdateSampleGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame,plotType=datatype)
    elif G2frame.GPXtree.GetItemText(item) == 'Index Peak List':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        data = G2frame.GPXtree.GetItemPyData(item)
#patch
        if len(data) != 2:
            data = [data,[]]
            G2frame.GPXtree.SetItemPyData(item,data)
#end patch
        G2pdG.UpdateIndexPeaksGrid(G2frame,data)
        if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Unit Cells List':
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        if not data:
            data.append([0,0.0,4,25.0,0,'P1',1.,1.,1.,90.,90.,90.,1.,'P 1']) #zero error flag, zero value, max Nc/No, start volume
            data.append([0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            data.append({})                                 #empty superlattice stuff
            data.append([])                                 #empty mag cells list
            G2frame.GPXtree.SetItemPyData(item,data)                             
#patch
        if len(data) < 5:
            data.append({'Use':False,'ModVec':[0,0,0.1],'maxH':1,'ssSymb':''})                                 #empty superlattice stuff
            G2frame.GPXtree.SetItemPyData(item,data)  
#end patch
        # debug stuff
        #if GSASIIpath.GetConfigValue('debug'):
        #    import imp
        #    imp.reload(G2pdG)  # for testing changes
        #    print('debug: reloaded',G2pdG)
        # end debug stuff
        G2pdG.UpdateUnitCellsGrid(G2frame,data)
        if 'PKS' in G2frame.GPXtree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Reflection Lists':   #powder reflections
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        data = G2frame.GPXtree.GetItemPyData(item)
        G2frame.RefList = ''
        if len(data):
            G2frame.RefList = list(data.keys())[0]
        G2pdG.UpdateReflectionGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.GPXtree.GetItemText(item) == 'Reflection List':    #HKLF reflections
        G2frame.PatternId = G2frame.GPXtree.GetItemParent(item)
        name = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
        G2pdG.UpdateReflectionGrid(G2frame,data,HKLF=True,Name=name)

    if G2frame.PickId:
        G2frame.PickIdText = G2frame.GetTreeItemsList(G2frame.PickId)
    # window has been filled, now resize scroll bars
    G2frame.dataWindow.SetDataSize()

    G2frame.Raise()
    if oldFocus:
        oldFocus.GetTopLevelParent().Raise()
        oldFocus.SetFocus()
    #debug code
    # print 'Got here!'
    # def FillWindow(panel,size=1.):
    #         sizer = panel.GetSizer()
    #         sizer.Add(wx.StaticText(panel, wx.ID_ANY, "Panel Two: long line "+int(size*40)*'*', (5,5)))
    #         for i in range(int(size*15)): sizer.Add(wx.StaticText(panel, wx.ID_ANY, "Line "+str(2+i), (5,5)))
    #         panel.Layout()
    #         panel.SendSizeEvent()
    # G2frame.dataWindow.ClearData()
    # FillWindow(G2frame.dataWindow)
        
def SetDataMenuBar(G2frame,menu=None):
    '''Set the menu for the data frame.

    Note that data frame items do not have menus, for these (menu=None)
    display the standard main menu for the data tree window.
    '''
    if menu is None:
        G2frame.SetMenuBar(G2frame.GSASIIMenu)
    else:
        G2frame.SetMenuBar(menu)
