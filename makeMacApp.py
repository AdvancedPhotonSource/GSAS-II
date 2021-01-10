#!/usr/bin/env python
'''
*makeMacApp: Create Mac Applet*
===============================

This script creates an AppleScript app bundle to launch GSAS-II. The app is
usually created in the directory where the GSAS-II script (.../GSASII/GSASII.py) 
is located. A softlink to Python is created inside that app bundle, 
but the softlink name is GSAS-II so that "GSAS-II" shows up as the name 
of the app in the menu bar, etc. rather than "Python". A soft link named 
GSAS-II.py, referencing the GSASII.py script, is created so that some file 
menu items also are labeled with GSAS-II (but not the right capitalization, 
alas). 

This can be used two different ways. 

 1. In the usual way, for conda-type installations
    where Python is in <condaroot>/bin and GSAS-II is in <condaroot>/GSASII, a premade 
    app is restored from a tar file. This works best for 11.0 (Big Sur) where there are security 
    constraints in place. 

 2. If python is not in that location or a name/location is specified
    for the app that will be created, this script creates an app (AppleScript) with the GSAS-II
    and the python locations hard coded. When an AppleScript is created,  
    this script tests to make sure that a wxpython script will run inside the 
    app and if not, it searches for a pythonw image and tries that. 

This has been tested with several versions of Python interpreters 
from Anaconda and does not require pythonw (Python.app). 

Run this script with no arguments or with one or two arguments.

The first argument, if supplied, is a reference to the GSASII.py script, 
which can have a relative or absolute path (the absolute path is determined).
If not supplied, the GSASII.py script will be used from the directory where 
this (makeMacApp.py) script is found. 

The second argument, if supplied, 
provides the name/location for the app to be created. This can be used to create 
multiple app copies using different Python versions (likely use for
development only). If the second argument is used, the AppleScript is created rather than 
restored from g2app.tar.gz
'''

from __future__ import division, print_function
import sys, os, os.path, stat, shutil, subprocess, plistlib
import platform
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>] [project]\n")
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
GSAS-II python script.
'''

if __name__ == '__main__':
    project="GSAS-II"
    # set GSAS-II location (path2GSAS): find the main GSAS-II script 
    # from this file, if not on command line
    if len(sys.argv) == 1:
        path2GSAS = os.path.split(__file__)[0]
        script = os.path.abspath(os.path.join(path2GSAS,"GSASII.py"))
    elif len(sys.argv) == 2:
        script = os.path.abspath(sys.argv[1])
        path2GSAS = os.path.split(script)[0]
    elif len(sys.argv) == 3:
        script = os.path.abspath(sys.argv[1])
        path2GSAS = os.path.split(script)[0]
        project = sys.argv[2]
    else:
        Usage()
    if not os.path.exists(script):
        print("\nFile "+script+" not found")
        Usage()
    if os.path.splitext(script)[1].lower() != '.py':
        print("\nScript "+script+" does not have extension .py")
        Usage()
    projectname = os.path.split(project)[1]
    if os.path.split(project)[0] != '':
        appPath = os.path.abspath(project+".app")
    else:
        appPath = os.path.abspath(os.path.join(path2GSAS,project+".app"))

# new approach, if possible use previously created file 
if __name__ == '__main__' and sys.platform == "darwin" and os.path.exists(
                os.path.join(path2GSAS,"g2app.tar.gz")) and project =="GSAS-II":
    if os.path.exists(os.path.join(path2GSAS,'../bin/python')):
        print('found python, found g2app.tar.gz')
        subprocess.call(["rm","-rf",appPath])
        subprocess.call(["mkdir","-p",appPath])
        subprocess.call(["tar","xzvf",os.path.join(path2GSAS,"g2app.tar.gz"),'-C',appPath])
        # create a link named GSAS-II.py to the script
        newScript = os.path.join(path2GSAS,'GSAS-II.py')
        if os.path.exists(newScript): # cleanup
            print("\nRemoving sym link",newScript)
            os.remove(newScript)
        os.symlink(os.path.split(script)[1],newScript)
        print("\nCreated "+projectname+" app ("+str(appPath)+
          ").\nViewing app in Finder so you can drag it to the dock if, you wish.")
        subprocess.call(["open","-R",appPath])
        sys.exit()
    else:
        print('found g2app.tar.gz, but python not in expected location')        

if __name__ == '__main__' and sys.platform == "darwin":
    iconfile = os.path.join(path2GSAS,'gsas2.icns') # optional icon file
    
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
    newScript = os.path.join(path2GSAS,'GSAS-II.py')
    if os.path.exists(newScript): # cleanup
        print("\nRemoving sym link",newScript)
        os.remove(newScript)
    os.symlink(os.path.split(script)[1],newScript)
    script=newScript

    # find Python used to run GSAS-II and set a new to use to call it
    # inside the app that will be created
    pythonExe = os.path.realpath(sys.executable)
    newpython =  os.path.join(appPath,"Contents","MacOS",projectname)
    
    # create an app using this python and if that fails to run wx, look for a
    # pythonw and try that with wx
    for i in 1,2,3:
        if os.path.exists(appPath): # cleanup
            print("\nRemoving old "+projectname+" app ("+str(appPath)+")")
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
        testout,errout = RunPython(newpython,'import numpy; import wx; wx.App(); print("-"+"OK-")')
        if isinstance(testout,bytes): testout = testout.decode()
        if "-OK-" in testout:
            print('wxpython app ran',testout)
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

    # Big Sur: open & save the file in the editor to set authorization levels
    osascript = '''
    tell application "Script Editor"
       set MyName to open "{}"
       save MyName
       (* close MyName *)
       (* quit *)
    end tell
'''.format(appPath)
    # detect MacOS 11 (11.0 == 10.16!)
    if platform.mac_ver()[0].split('.')[0] == '11' or platform.mac_ver()[0][:5] == '10.16':
        print("\nFor Big Sur and later, save the app in Script Editor before using it\n")
        subprocess.Popen(["osascript","-e",osascript])
    print("\nCreated "+projectname+" app ("+str(appPath)+
          ").\nViewing app in Finder so you can drag it to the dock if, you wish.")
    subprocess.call(["open","-R",appPath])
