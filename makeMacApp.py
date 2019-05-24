#!/usr/bin/env python
'''
*makeMacApp: Create Mac Applet*
===============================

This script creates an AppleScript app bundle to launch GSAS-II. The app is
created in the directory where the GSAS-II script (GSASII.py) is located.
A softlink to Python is created, but is named GSAS-II, and placed inside
the bundle so that GSAS-II shows up as the name of the app rather than 
Python in the menu bar, etc. A soft link named GSAS-II.py, referencing 
GSASII.py, is created so that appropriate menu items also are labeled with 
GSAS-II (but not the right capitalization, alas). 

This has been tested with several versions of Python interpreters 
from Anaconda and does not require pythonw (Python.app). It tests to 
make sure that a wx python script will run inside Python and if not, 
it searches for a pythonw image and tries that. 

Run this script with no arguments or with one optional argument, 
a reference to the GSASII.py, with a relative or absolute path. Either way 
the path will be converted via an absolute path. If no arguments are 
supplied, the GSASII.py script must be located in the same directory as 
this file.

'''
from __future__ import division, print_function
import sys, os, os.path, stat, shutil, subprocess, plistlib
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>]\n")
    sys.exit()

def RunPython(image,cmd):
    'Run a command in a python image'
    try:
        err=None
        p = subprocess.Popen([image,'-c',cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out = p.stdout.read()
        err = p.stderr.read()
        p.communicate()
        return out,err
    except Exception(err):
        return '','Exception = '+err

AppleScript = ''
'''Contains an AppleScript to start GSAS-II, launching Python and the
GSAS-II python script
'''

if __name__ == '__main__':
    project="GSAS-II"

    # set scriptdir: find the main GSAS-II script if not on command line
    if len(sys.argv) == 1:
        script = os.path.abspath(os.path.join(
                os.path.split(__file__)[0],
                "GSASII.py"
                ))
    elif len(sys.argv) == 2:
            script = os.path.abspath(sys.argv[1])
    else:
        Usage()
        raise Exception
    # make sure we found it
    if not os.path.exists(script):
        print("\nFile "+script+" not found")
        Usage()
    if os.path.splitext(script)[1].lower() != '.py':
        print("\nScript "+script+" does not have extension .py")
        Usage()
    # where the app will be created
    scriptdir = os.path.split(script)[0]
    
    appPath = os.path.abspath(os.path.join(scriptdir,project+".app"))
    iconfile = os.path.join(scriptdir,'gsas2.icns') # optional icon file
    
    AppleScript = '''(*   GSAS-II AppleScript by B. Toby (brian.toby@anl.gov)
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
	set env to "{:s}"
	TestFilePresent(appwithpath)
	TestFilePresent(python)
	tell application "Terminal"
		activate
		do script env & python & " " & appwithpath & "; exit"
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
	set env to "{:s}"
 
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
				do script env & python & " " & appwithpath & " " & filename & "; exit"
			end tell
		end if
	end repeat
end open
'''
    # create a link named GSAS-II.py to the script
    newScript = os.path.join(scriptdir,project+'.py')
    if os.path.exists(newScript): # cleanup
        print("\nRemoving sym link",newScript)
        os.remove(newScript)
    os.symlink(os.path.split(script)[1],newScript)
    script=newScript

    # find Python used to run GSAS-II and set a new to use to call it
    # inside the app that will be created
    pythonExe = os.path.realpath(sys.executable)
    newpython =  os.path.join(appPath,"Contents","MacOS",project)

    # create an app using this python and if that fails to run wx, look for a
    # pythonw and try that with wx
    for i in 1,2,3:
        if os.path.exists(appPath): # cleanup
            print("\nRemoving old "+project+" app ("+str(appPath)+")")
            shutil.rmtree(appPath)
        
        shell = os.path.join("/tmp/","appscrpt.script")
        f = open(shell, "w")
        f.write(AppleScript.format(newpython,script,'',newpython,script,''))
        f.close()

        try: 
            subprocess.check_output(["osacompile","-o",appPath,shell],stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as msg:
            print('''Error compiling AppleScript.
            Report the next message along with details about your Mac to toby@anl.gov''')
            print(msg.output)
            sys.exit()
        # create a link to the python inside the app, if named to match the project
        if pythonExe != newpython: os.symlink(pythonExe,newpython)

        # test if newpython can run wx
        testout,errout = RunPython(newpython,'import numpy; import wx; wx.App(); print("-OK-")')
        if isinstance(testout,bytes): testout = testout.decode()
        if "-OK-" in testout:
            break
        elif i == 1:
            print('Run of wx in',pythonExe,'failed, looking for pythonw')
            pythonExe = os.path.join(os.path.split(sys.executable)[0],'pythonw')
            if not os.path.exists(pythonExe):
                print('Warning no pythonw found with ',sys.executable,
                      '\ncontinuing, hoping for the best')
        elif i == 2:
            print('Warning could not run wx with',pythonExe,
                      'will try with that external to app')
            newpython = pythonExe
        else:
            print('Warning still could not run wx with',pythonExe,
                      '\ncontinuing, hoping for the best')

    # change the icon
    oldicon = os.path.join(appPath,"Contents","Resources","droplet.icns")
    #if os.path.exists(iconfile) and os.path.exists(oldicon):
    if os.path.exists(iconfile):
        shutil.copyfile(iconfile,oldicon)

    # Edit the app plist file to restrict the type of files that can be dropped
    if hasattr(plistlib,'load'):
        fp = open(os.path.join(appPath,"Contents",'Info.plist'),'rb')
        d = plistlib.load(fp)
        fp.close()
    else:
        d = plistlib.readPlist(os.path.join(appPath,"Contents",'Info.plist'))
    d['CFBundleDocumentTypes'] = [{
        'CFBundleTypeExtensions': ['gpx'],
        'CFBundleTypeName': 'GSAS-II project',
        'CFBundleTypeRole': 'Editor'}]
    
    if hasattr(plistlib,'dump'):
        fp = open(os.path.join(appPath,"Contents",'Info.plist'),'wb')
        plistlib.dump(d,fp)
        fp.close()
    else:
        plistlib.writePlist(d,os.path.join(appPath,"Contents",'Info.plist'))

    print("\nCreated "+project+" app ("+str(appPath)+
          ").\nViewing app in Finder so you can drag it to the dock if, you wish.")
    subprocess.call(["open","-R",appPath])
