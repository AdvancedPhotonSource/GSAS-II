# creates Windows files to aid in running GSAS-II
'''
This script creates a file named ``RunGSASII.bat`` and a desktop shortcut to that file.
It registers the filetype .gpx so that the GSAS-II project files exhibit the
GSAS-II icon and so that double-clicking on them opens them in GSAS-II. 

Run this script with no arguments; the path to the ``GSASII.py`` file
is assumed to be in the parent directory to the one where this file 
(``makeBat.py``) is found.

The contents of this file may also be run from inside the gitstrap.py
installation script. In that case, the following variables are already 
defined: 

  * path2GSAS2   is the directory with all GSAS-II Python code
  * G2script     has the location of the GSASII.py file
  * path2repo    is the location of the GSAS-II git repository

The path to Python is determined from the version of Python used to 
run this script.
'''
#
import os, sys
import datetime

Script = '''@REM Script to start GSAS-II on Windows
@echo =========================================================================
@echo                General Structure Analysis System-II
@echo              by Robert B. Von Dreele and Brian H. Toby
@echo              Argonne National Laboratory(C), 2006-2014
@echo  This product includes software developed by the UChicago Argonne, LLC,
@echo             as Operator of Argonne National Laboratory.
@echo                            Please cite:
@echo      B.H. Toby and R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013)
@echo  + other papers for DIFFax, small angle, Bilbao, ISODISTORT,... as shown
@echo =========================================================================
@
{:s}{:s} {:s} "%~1"
@REM To keep the window from disappearing with any error messages
pause

'''
    
app = None # delay starting wx until we need it. Likely not needed. 
if __name__ == '__main__':
    try:
        import winreg
    except ImportError:
        try:
            import _winreg as winreg
        except ImportError:
            print('winreg not found')

    if __file__.lower().endswith("makebat.py"):
        print(f"running from file {__file__!r}")
        invokedDirectly = True
        path2GSAS2 = os.path.dirname(os.path.dirname(__file__))
        path2repo = os.path.dirname(path2GSAS2)
        G2script = os.path.join(path2GSAS2,'GSASII.py')
    else:
        print(f"running makeBat.py indirectly inside {__file__!r}")
        invokedDirectly = False

    G2icon = os.path.join(path2GSAS2,'icons','gsas2.ico')
    # create the .BAT file below the git directory, where gitstrap is located
    G2bat = os.path.normpath(os.path.join(path2repo,'..','RunGSASII.bat'))
    pythonexe = os.path.realpath(sys.executable)

    print('Python installed at ',pythonexe)
    print('GSAS-II installed at',path2GSAS2)
    print('GSASII.py at        ',G2script)
    print('GSASII icon at      ',G2icon)
    print('.bat file to be at  ',G2bat)

    # create a GSAS-II script
    fp = open(G2bat,'w')
    fp.write("@REM created by run of makeBat.py on {:%d %b %Y %H:%M}\n".format(
        datetime.datetime.now()))
    activate = os.path.join(os.path.split(pythonexe)[0],'Scripts','activate')
    print("Looking for",activate)
    # for a non-base conda install, it might be better to use the activate in
    # the base, but for now let's use the one we find relative to our python
    if os.path.exists(activate):
        activate = os.path.realpath(activate)
        if ' ' in activate:
            activate = 'call "'+ activate + '"\n'
        else:
            activate = 'call '+ activate + '\n'
        print(f'adding activate to .bat file ({activate})')
    else:
        print('conda activate not found')
        activate = ''
    pexe = pythonexe
    if ' ' in pythonexe: pexe = '"'+pythonexe+'"'
    G2s = G2script
    if ' ' in G2s: G2s = '"'+G2script+'"'
    fp.write(Script.format(activate,pexe,G2s))
    fp.close()
    print(f'\nCreated GSAS-II batch file {G2bat}')
    
    # create a reset script
    gitstrap = os.path.abspath(
        os.path.normpath(os.path.join(path2repo,'..','gitstrap.py')))
    if not os.path.exists:
        print(f'the installation script was not found: {gitstrap!r}')
    else:
        G2reset = os.path.normpath(os.path.join(path2repo,'..','Reset2FreshGSASII.bat'))
        fp = open(G2reset,'w')
        fp.write("@REM created by run of makeBat.py on {:%d %b %Y %H:%M}\n".format(
            datetime.datetime.now()))
        fp.write("REM This script will reset GSAS-II to the latest version, even if the program can't be started\n")
        pexe = pythonexe
        if ' ' in pythonexe: pexe = '"'+pythonexe+'"'
        G2s = gitstrap
        if ' ' in G2s: G2s = '"'+gitstrap+'"'
        if activate: fp.write(f"{activate}")
        fp.write('choice /c yn /n /m "Reset any local changes and install latest GSAS-II version? (y/n)"\n')
        fp.write(f"goto %ERRORLEVEL%\n")
        fp.write(f":1\n")
        fp.write(f"{pexe} {G2s} --reset\n")
        fp.write(f":2\n")
        fp.write(f"pause\n")
        fp.close()
        print(f'\nCreated GSAS-II reset script {G2reset}')

    new = False
    oldBat = ''
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
            print(f'old GPX assignment {oldBat} not found; registry entry will be made for new one')
        new = True
    except NameError:
        pass
    if invokedDirectly and not new:
        try:
            if oldBat != G2bat:
                if app is None:
                    import wx
                    app = wx.App()
                dlg = wx.MessageDialog(None,'gpx files already assigned in registry to: \n'+oldBat+'\n Replace with: '+G2bat+'?','GSAS-II gpx in use', 
                        wx.YES_NO | wx.ICON_QUESTION | wx.STAY_ON_TOP)
                dlg.Raise()
                if dlg.ShowModal() == wx.ID_YES:
                    new = True
                dlg.Destroy()
        finally:
            pass
    elif not invokedDirectly:  # force a choice, since we can't ask
        new = True
    if new:
        # Associate a script and icon with .gpx files
        try:
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
            print('Assigned icon and batch file to .gpx files in registry')
        except:
            print('Error assigning icon and batch file to .gpx files')
    else:
        print('old assignment of icon and batch file in registery is retained')

    try:
        import win32com.shell.shell, win32com.shell.shellcon
        win32com.shell.shell.SHChangeNotify(
            win32com.shell.shellcon.SHCNE_ASSOCCHANGED, 0, None, None)
    except ImportError:
        print('Module pywin32 not present, .gpx files will display properly after logging in again')
    except:
        print('Unexpected error on explorer refresh.')
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
            if invokedDirectly:
                if app is None:
                    import wx
                    app = wx.App()
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
            else:
                # set an installation location
                distdir = os.path.split(
                    os.path.dirname(os.path.dirname(path2GSAS2))
                    )[1]
                if distdir == '\\' or distdir == '': distdir = '/'
                shortbase = f"GSAS-II from {distdir}.lnk"
                shortcut = os.path.join(desktop, shortbase)
        if save:
            shell = win32com.client.Dispatch('WScript.Shell')
            shobj = shell.CreateShortCut(shortcut)
            shobj.Targetpath = G2bat
            #shobj.WorkingDirectory = wDir # could specify a default project location here
            shobj.IconLocation = G2icon
            shobj.save()
            print(f'Created shortcut {shortbase!r} to start GSAS-II on desktop')
        else:
            print('No shortcut for this GSAS-II created on desktop')
    except ImportError:
        print('Module pywin32 not present, will not make desktop shortcut')
    except:
        print('Unexpected error making desktop shortcut.')
        import traceback
        print(traceback.format_exc())
