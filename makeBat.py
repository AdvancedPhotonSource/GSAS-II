# creates Windows files to aid in running GSAS-II
#   creates RunGSASII.bat and a desktop shortcut to that file
#   registers the filetype .gpx so that the GSAS-II project files exhibit the
#     GSAS-II icon and so that double-clicking on them opens them in GSAS-II
import os, sys
import datetime
try:
    import _winreg as winreg
except ImportError:
    import winreg

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
@echo ========================================================================
@
{:s} {:s} %1
@REM To keep the window from disappearing with any error messages
pause

'''
gsaspath = os.path.split(sys.argv[0])[0]
if not gsaspath: gsaspath = os.path.curdir
gsaspath = os.path.abspath(os.path.expanduser(gsaspath))
G2script = os.path.join(gsaspath,'GSASII.py')
G2bat = os.path.join(gsaspath,'RunGSASII.bat')
G2icon = os.path.join(gsaspath,'gsas2.ico')
pythonexe = os.path.realpath(sys.executable)
# create a GSAS-II script
fp = open(os.path.join(G2bat),'w')
fp.write("@REM created by run of bootstrap.py on {:%d %b %Y %H:%M}\n".format(
    datetime.datetime.now()))
pexe = pythonexe
if ' ' in pythonexe: pexe = '"'+pythonexe+'"'
G2s = G2script
if ' ' in G2s: G2script = '"'+G2script+'"'
fp.write(Script.format(pexe,G2s))
fp.close()
print '\nCreated GSAS-II batch file RunGSASII.bat in '+gsaspath

# Associate a script and icon with .gpx files
#gpxkey = winreg.CreateKey(winreg.HKEY_CLASSES_ROOT, '.gpx')
gpxkey = winreg.CreateKey(winreg.HKEY_CURRENT_USER,r'Software\CLASSES\.gpx')
winreg.SetValue(gpxkey, None, winreg.REG_SZ, 'GSAS-II.project')
winreg.CloseKey(gpxkey)
gpxkey = winreg.CreateKey(winreg.HKEY_CURRENT_USER,r'Software\CLASSES\GSAS-II.project')
winreg.SetValue(gpxkey, None, winreg.REG_SZ, 'GSAS-II project')
iconkey = winreg.CreateKey(gpxkey, 'DefaultIcon')
winreg.SetValue(iconkey, None, winreg.REG_SZ, G2icon)
openkey = winreg.CreateKey(gpxkey, r'shell\open\command')
winreg.SetValue(openkey, None, winreg.REG_SZ, G2bat+" %1")
winreg.CloseKey(iconkey)
winreg.CloseKey(openkey)
winreg.CloseKey(gpxkey)
print 'Assigned icon and batch file to .gpx files'

try:
    import win32com.shell.shell, win32com.shell.shellcon
    win32com.shell.shell.SHChangeNotify(
        win32com.shell.shellcon.SHCNE_ASSOCCHANGED, 0, None, None)
except ImportError:
    print 'Module pywin32 not present, login again to see file types properly'
except:
    print 'Unexpected error on explorer refresh. Please report:'
    import traceback
    print traceback.format_exc()

# make a desktop shortcut to GSAS-II
try:
    import win32com.shell.shell, win32com.shell.shellcon, win32com.client
    desktop = win32com.shell.shell.SHGetFolderPath(
        0, win32com.shell.shellcon.CSIDL_DESKTOP, None, 0)
    shortcut = os.path.join(desktop, "GSAS-II.lnk")
    shell = win32com.client.Dispatch('WScript.Shell')
    shobj = shell.CreateShortCut(shortcut)
    shobj.Targetpath = G2bat
    #shobj.WorkingDirectory = wDir # could specify a default project location here
    shobj.IconLocation = G2icon
    shobj.save()
    print 'Created shortcut to start GSAS-II on desktop'
except ImportError:
    print 'Module pywin32 not present, will not make desktop shortcut'
except:
    print 'Unexpected error making desktop shortcut. Please report:'
    import traceback
    print traceback.format_exc()
    
