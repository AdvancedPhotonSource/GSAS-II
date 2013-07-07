#!/usr/bin/env python
'''This script creates an AppleScript app to launch GSAS-II. The app is
created in the directory where the GSAS-II script is located.
A softlink to Python is created, but is named GSAS-II, so that 
GSAS-II shows up as the name of the app rather than Python in the
menu bar, etc. Note that this requires finding an app version of Python
(expected name .../Resources/Python.app/Contents/MacOS/Python in
directory tree of the calling python interpreter).

Run this script with one optional argument, the path to the GSASII.py
The script path may be specified relative to the current path or given
an absolute path, but will be accessed via an absolute path. 
If no arguments are supplied, the GSASII.py script is assumed to be in the
same directory as this file.

'''
import sys, os, os.path, stat, shutil, subprocess, plistlib
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>]\n")
    sys.exit()

project="GSAS-II"

# find the main GSAS-II script if not on command line
if len(sys.argv) == 1:
    script = os.path.abspath(os.path.join(
        os.path.split(__file__)[0],
        "GSASII.py"
        ))
elif len(sys.argv) == 2:
    script = os.path.abspath(sys.argv[1])
else:
    Usage()

# make sure we found it
if not os.path.exists(script):
    print("\nFile "+script+" not found")
    Usage()
if os.path.splitext(script)[1].lower() != '.py':
    print("\nScript "+script+" does not have extension .py")
    Usage()
# where the app will be created
scriptdir = os.path.split(script)[0]
# name of app
apppath = os.path.abspath(os.path.join(scriptdir,project+".app"))
iconfile = os.path.join(scriptdir,'gsas2.icns') # optional icon file

# find the python application; must be an OS X app
pythonpath,top = os.path.split(os.path.realpath(sys.executable))
while top:
    if 'Resources' in pythonpath:
        pass
    elif os.path.exists(os.path.join(pythonpath,'Resources')):
        break
    pythonpath,top = os.path.split(pythonpath)
else:
    print("\nSorry, failed to find a Resources directory associated with "+str(sys.executable))
    sys.exit()
pythonapp = os.path.join(pythonpath,'Resources','Python.app','Contents','MacOS','Python')
if not os.path.exists(pythonapp): 
    print("\nSorry, failed to find a Python app in "+str(pythonapp))
    sys.exit()

# new name to call python
newpython =  os.path.join(apppath,"Contents","MacOS",project)

if os.path.exists(apppath): # cleanup
    print("\nRemoving old "+project+" app ("+str(apppath)+")")
    shutil.rmtree(apppath)

# create an AppleScript that launches python with the requested app
shell = os.path.join("/tmp/","appscrpt.script")
f = open(shell, "w")
f.write('''(*   GSAS-II AppleScript by B. Toby (brian.toby@anl.gov)
     It can launch GSAS-II by double clicking or by dropping a data file
     or folder over the app.
     It runs GSAS-II in a terminal window.
*)

(* test if a file is present and exit with an error message if it is not  *)
on TestFilePresent(appwithpath)
	tell application "System Events"
		if (file appwithpath exists) then
		else
			display dialog "Error: file " & appwithpath & " not found. If you have moved this file recreate the AppleScript with bootstrap.py." with icon caution buttons {{"Quit"}}
			return
		end if
	end tell
end TestFilePresent

(* 
------------------------------------------------------------------------
this section responds to a double-click. No file is supplied to GSAS-II
------------------------------------------------------------------------ 
*)
on run
	set python to "{:s}"
	set appwithpath to "{:s}"
	TestFilePresent(appwithpath)
	TestFilePresent(python)
	tell application "Terminal"
		activate
		do script python & " " & appwithpath & "; exit"
	end tell
end run

(*
-----------------------------------------------------------------------------------------------
this section handles starting with files dragged into the AppleScript
 o it goes through the list of file(s) dragged in
 o then it converts the colon-delimited macintosh file location to a POSIX filename
 o for every non-directory file dragged into the icon, it starts GSAS-II, passing the file name
------------------------------------------------------------------------------------------------
*)

on open names
	set python to "{:s}"
	set appwithpath to "{:s}"
	TestFilePresent(appwithpath)
	repeat with filename in names
		set filestr to (filename as string)
		if filestr ends with ":" then
                        (* should not happen, skip directories *)
		else
			(* if this is an input file, open it *)
			set filename to the quoted form of the POSIX path of filename
			tell application "Terminal"
				activate
				do script python & " " & appwithpath & " " & filename & "; exit"
			end tell
		end if
	end repeat
end open
'''.format(newpython,script,newpython,script))
f.close()

try: 
    subprocess.check_output(["osacompile","-o",apppath,shell],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError, msg:
    print '''Error compiling AppleScript.
    Report the next message along with details about your Mac to toby@anl.gov'''
    print msg.output
    sys.exit()

# create a link to the python app, but named to match the project
os.symlink(pythonapp,newpython)

# change the icon
oldicon = os.path.join(apppath,"Contents","Resources","droplet.icns")
if os.path.exists(iconfile) and os.path.exists(oldicon):
    shutil.copyfile(iconfile,oldicon)

# Edit the app plist file to restrict the type of files that can be dropped
d = plistlib.readPlist(os.path.join(apppath,"Contents",'Info.plist'))
d['CFBundleDocumentTypes'] = [{
    'CFBundleTypeExtensions': ['gpx'],
    'CFBundleTypeName': 'GSAS-II project',
    'CFBundleTypeRole': 'Editor'}]
plistlib.writePlist(d,os.path.join(apppath,"Contents",'Info.plist'))

print("\nCreated "+project+" app ("+str(apppath)+
      ").\nViewing app in Finder so you can drag it to the dock if, you wish.")
subprocess.call(["open","-R",apppath])
