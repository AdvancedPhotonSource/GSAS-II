'''
*makeBat: Create GSAS-II Batch File*
====================================

This script creates a file named ``RunGSASII.bat`` and a desktop shortcut to that file.
It registers the filetype .gpx so that the GSAS-II project files exhibit the
GSAS-II icon and so that double-clicking on them opens them in GSAS-II. 

Run this script with no arguments; the path to the ``GSASII.py`` file
is assumed to be the the same as the path to the ``makeBat.py`` file
and the path to Python is determined from the version of Python
used to run this script. 

'''
from __future__ import division, print_function
version = "$Id$"
# creates Windows files to aid in running GSAS-II
#   creates RunGSASII.bat and a desktop shortcut to that file
#   registers the filetype .gpx so that the GSAS-II project files exhibit the
#     GSAS-II icon and so that double-clicking on them opens them in GSAS-II
#
import os, sys
import datetime
import wx

Script = '''@echo ========================================================================
@echo                General Structure Analysis System-II
@echo              by Robert B. Von Dreele and Brian H. Toby
@echo                Argonne National Laboratory(C), 2010
@echo  This product includes software developed by the UChicago Argonne, LLC,
@echo             as Operator of Argonne National Laboratory.
@echo                            Please cite:
@echo      B.H. Toby and R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013)
@echo                   for small angle use also cite:
@echo      R.B. Von Dreele, J. Appl. Cryst. 47, 1784-9 (2014)
@echo                   for DIFFaX use also cite:
@echo      M.M.J. Treacy, J.M. Newsam and M.W. Deem, 
@echo                   Proc. Roy. Soc. Lond. 433A, 499-520 (1991)
@echo ========================================================================
@
{:s}{:s} {:s} "%~1"
@REM To keep the window from disappearing with any error messages
pause

'''

if __name__ == '__main__':
    try:
        import _winreg as winreg
    except ImportError:
        import winreg
    app = None # delay starting wx until we need it. Likely not needed. 
    gsaspath = os.path.split(sys.argv[0])[0]
    if not gsaspath: gsaspath = os.path.curdir
    gsaspath = os.path.abspath(os.path.expanduser(gsaspath))
    G2script = os.path.join(gsaspath,'GSASII.py')
    G2bat = os.path.join(gsaspath,'RunGSASII.bat')
    G2icon = os.path.join(gsaspath,'gsas2.ico')
    pythonexe = os.path.realpath(sys.executable)
    print('Python installed at',pythonexe)
    print('GSAS-II installed at',gsaspath)
    # Bob reports a problem using pythonw.exe w/Canopy on Windows, so change that if used
    if pythonexe.lower().endswith('pythonw.exe'):
        print("  using python.exe rather than "+pythonexe)
        pythonexe = os.path.join(os.path.split(pythonexe)[0],'python.exe')
        print("  now pythonexe="+pythonexe)
    # create a GSAS-II script
    fp = open(os.path.join(G2bat),'w')
    fp.write("@REM created by run of bootstrap.py on {:%d %b %Y %H:%M}\n".format(
        datetime.datetime.now()))
    activate = os.path.join(os.path.split(pythonexe)[0],'Scripts','activate')
    print("Looking for",activate)
    if os.path.exists(activate):
        activate = os.path.realpath(activate)
        if ' ' in activate:
            activate = 'call "'+ activate + '"\n'
        else:
            activate = 'call '+ activate + '\n'
        print('adding activate to .bat file ({})'.format(activate))
    else:
        print('Anaconda activate not found')
        activate = ''
    pexe = pythonexe
    if ' ' in pythonexe: pexe = '"'+pythonexe+'"'
    G2s = G2script
    if ' ' in G2s: G2s = '"'+G2script+'"'
    # is mingw-w64\bin present? If so add it to path.
    #d = os.path.split(pexe)[0]
    #mdir = os.path.join(d,'Library','mingw-w64','bin')
    #if os.path.exists(mdir):
    #    fp.write('@path={};%path%\n'.format(mdir))
    fp.write(Script.format(activate,pexe,G2s))
    fp.close()
    print('\nCreated GSAS-II batch file RunGSASII.bat in '+gsaspath)
    
    new = False
    oldBat = ''
    try: # patch for FileNotFoundError not in Python 2.7
        FileNotFoundError
    except NameError:
        FileNotFoundError = Exception
    # this code does not appear to work properly when paths have spaces
    try:
        oldgpx = winreg.OpenKey(winreg.HKEY_CURRENT_USER,r'Software\CLASSES\GSAS-II.project') # throws FileNotFoundError
        oldopen = winreg.OpenKey(oldgpx,r'shell\open\command')
        # get previous value & strip %1 off end
        oldBat = winreg.QueryValue(oldopen,None).strip()
        pos = oldBat.rfind(' ')
        if pos > 1:
            oldBat = oldBat[:pos]
        os.stat(oldBat)     #check if it is still around
    except FileNotFoundError:
        if oldBat:
            print('old GPX assignment',oldBat, 'not found; registry entry will be made for new one')
        new = True
    if not new:
        try:
            if oldBat != G2bat:
                if app is None:
                    app = wx.App()
                    app.MainLoop()
                dlg = wx.MessageDialog(None,'gpx files already assigned in registry to: \n'+oldBat+'\n Replace with: '+G2bat+'?','GSAS-II gpx in use', 
                        wx.YES_NO | wx.ICON_QUESTION | wx.STAY_ON_TOP)
                dlg.Raise()
                if dlg.ShowModal() == wx.ID_YES:
                    new = True
                dlg.Destroy()
        finally:
            pass
    if new:
        # Associate a script and icon with .gpx files
        gpxkey = winreg.CreateKey(winreg.HKEY_CURRENT_USER,r'Software\CLASSES\.gpx')
        winreg.SetValue(gpxkey, None, winreg.REG_SZ, 'GSAS-II.project')
        winreg.CloseKey(gpxkey)
        gpxkey = winreg.CreateKey(winreg.HKEY_CURRENT_USER,r'Software\CLASSES\GSAS-II.project')
        winreg.SetValue(gpxkey, None, winreg.REG_SZ, 'GSAS-II project')
        iconkey = winreg.CreateKey(gpxkey, 'DefaultIcon')
        winreg.SetValue(iconkey, None, winreg.REG_SZ, G2icon)
        openkey = winreg.CreateKey(gpxkey, r'shell\open\command')
        winreg.SetValue(openkey, None, winreg.REG_SZ, G2bat+' "%1"')
        winreg.CloseKey(iconkey)
        winreg.CloseKey(openkey)
        winreg.CloseKey(gpxkey)
        print('Assigned icon and batch file to .gpx files in registery')
    else:
        print('old assignment of icon and batch file in registery is retained')

    try:
        import win32com.shell.shell, win32com.shell.shellcon
        win32com.shell.shell.SHChangeNotify(
            win32com.shell.shellcon.SHCNE_ASSOCCHANGED, 0, None, None)
    except ImportError:
        print('Module pywin32 not present, login again to see file types properly')
    except:
        print('Unexpected error on explorer refresh. Please report:')
        import traceback
        print(traceback.format_exc())

    # make a desktop shortcut to GSAS-II
    try:
        import win32com.shell.shell, win32com.shell.shellcon, win32com.client
        desktop = win32com.shell.shell.SHGetFolderPath(
            0, win32com.shell.shellcon.CSIDL_DESKTOP, None, 0)
        shortbase = "GSAS-II.lnk"
        shortcut = os.path.join(desktop, shortbase)
        save = True
        if win32com.shell.shell.SHGetFileInfo(shortcut,0,0)[0]:
            print('GSAS-II shortcut exists!')
            if app is None:
                app = wx.App()
                app.MainLoop()
            dlg = wx.FileDialog(None, 'Choose new GSAS-II shortcut name',  desktop, shortbase,
                wildcard='GSAS-II shortcut (*.lnk)|*.lnk',style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            dlg.Raise()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    shortcut = dlg.GetPath()
                else:
                    save = False
            finally:
                dlg.Destroy()
        if save:
            shell = win32com.client.Dispatch('WScript.Shell')
            shobj = shell.CreateShortCut(shortcut)
            shobj.Targetpath = G2bat
            #shobj.WorkingDirectory = wDir # could specify a default project location here
            shobj.IconLocation = G2icon
            shobj.save()
            print('Created shortcut to start GSAS-II on desktop')
        else:
            print('No shortcut for this GSAS-II created on desktop')
    except ImportError:
        print('Module pywin32 not present, will not make desktop shortcut')
    except:
        print('Unexpected error making desktop shortcut. Please report:')
        import traceback
        print(traceback.format_exc())
    
