#!/usr/bin/env python
#------------------------------------------------------------
# this version intended for use with git installations
#------------------------------------------------------------
'''
This routine creates an app bundle named GSAS-II.app. Inside the 
bundle is a symbolic link to the Python executable named "GSAS-II" 
that will be used to run GSAS-II. Having this link named that
way causes the name of the app to shows in the menu bar as 
"GSAS-II" rather than "Python". Also used by the app, is another 
symbolic link named GSAS-II.py, which must be placed in the same 
directory as the app bundle. This file is linked to the GSASII.py 
script and the link is run using the link to Python. This also 
causes other items in the app to be labeled as GSAS-II (but not 
with the right capitalization, alas). 

The original contents of the app bundle was created interactively 
and, after some manual edits, the contents of that was placed into 
a tar file distributed with GSAS-II, and is expanded in this script. 
This method seems to be needed for MacOS 11.0+ (Big Sur and later) 
where Apple's security constraints seem to prevent creation of the 
app directly. Older code (not currently in use) created the 
app from "scratch" using the osacompile utility, but that no longer
seems to work.

Three different paths are needed to run this script:: 

    path2GSAS:  The location where the GSASII.py (and other GSAS-II 
       Python files) are found. 
    installLoc: The location where the GSAS-II.app app bundle and 
       the GSAS-II.py will be placed.
    pythonLoc:  The location of the Python executable. 

Under normal circumstances, the locations for all of these paths 
can be determined from the location of the makeMacApp.py file. 
Note that when GSAS-II is installed from git using gitstrap.py, 
the git repository is placed at <loc>/GSAS-II and the GSAS-II 
Python scripts are placed at the GSASII child directory, so that 
GSAS-II is started from the GSASII.py script at <loc>/GSAS-II/GSASII/
and the current script (makeMacApp.py) will be found in 
<loc>/GSAS-II/GSASII/install/. 

When the GSAS-II conda installers 
are used, the git repository is placed at $CONDA_HOME/GSAS-II so that 
<loc> above is $CONDA_HOME. Also, the Python executable will be found
in $CONDA_HOME/bin/Python. Thus, if this file is in makePath (typically 
<loc>/GSAS-II/GSASII/install), then 

    * path2GSAS will be makePath/.. and 
    * installLoc will be path2GSAS/.. and 
    * pythonLoc will be installLoc/../bin/python, 

but these locations can be overridden from the command-line arguments. 
If a Python location is not supplied and is not at the default location 
(installLoc/../bin/python) then the Python executable currently 
running this script (from sys.executable) is used. 
    
Run this script with no arguments or with one or two arguments.

The first argument, if supplied, provides the path to be used for the 
app bundle will be created. Note that GSAS-II.app and GSAS-II.py will
be created in this directory. 

The second argument, if supplied, is path2GSAS, a path to the 
location GSASII.py script, which can be a relative path 
(the absolute path is determined). If not supplied, the GSASII.py script 
is expected to be located in the directory above where this 
(makeMacApp.py) script is found. 

The third argument, if supplied, provides the full path for the Python 
installation to be used inside the app bundle that will be created. If not 
supplied, and Python exists at installLoc/../bin/python, that will be used. 
If that does not exist, then the location of the current Python executable
(from sys.executable) will be used. 

'''

from __future__ import division, print_function
import sys, os, os.path, stat, shutil, subprocess, plistlib
import platform
def Usage():
    print(f"\nUsage:\n\tpython {os.path.abspath(sys.argv[0])} [install_path] [<GSAS-II_script>] [Python_loc]\n")
    sys.exit()

AppleScript = ''
'''Will be set to contain an AppleScript to start GSAS-II by launching 
Python and the GSAS-II Python script. Not currently used.
'''

if __name__ == '__main__':
    # set defaults
    project="GSAS-II"  # name of app
    makePath = os.path.dirname(__file__)  # location of this script
    path2GSAS = os.path.dirname(makePath)
    installLoc = os.path.dirname(path2GSAS)
    pythonLoc = None
    # use command line args, if any
    if len(sys.argv) > 4: # too many
        Usage()
    if len(sys.argv) >= 2:
        installLoc = os.path.abspath(sys.argv[1])
    if len(sys.argv) >= 3:
        path2GSAS = os.path.abspath(sys.argv[2])
    if len(sys.argv) == 4:
        pythonLoc = os.path.abspath(sys.argv[3])
    # sanity checking
    G2script = os.path.abspath(os.path.join(path2GSAS,"GSASII.py"))
    if not os.path.exists(G2script):
        print(f"\nERROR: File {G2script!r} not found")
        Usage()
    if os.path.splitext(G2script)[1].lower() != '.py':
        print(f"\nScript {G2script!r} does not have extension .py")
        Usage()
    if not os.path.exists(installLoc):
        print(f"\nERROR: directory {installLoc!r} not found")
        Usage()
    tarLoc = os.path.join(path2GSAS,'install',"g2app.tar.gz")
    if not os.path.exists(tarLoc):
        print(f"\nERROR: file {tarLoc!r} not found")
        Usage()
    if pythonLoc is None:
        pythonLoc = os.path.join(installLoc,'../bin',"python")
        if not os.path.exists(pythonLoc):
            pythonLoc = sys.executable
    if not os.path.exists(pythonLoc):
        print(f"\nERROR: Python not found at {pythonLoc!r}")
        Usage()

    print(f'Using Python: {pythonLoc}')
    print(f'Using GSAS-II script: {G2script}')
    print(f'Install location: {installLoc}')

    # files to be created
    appName = os.path.abspath(os.path.join(installLoc,project+".app"))
    g2Name = os.path.abspath(os.path.join(installLoc,project+'.py'))

# new approach, use previously created tar (.tgz) file 
if __name__ == '__main__' and sys.platform == "darwin":
    if os.path.exists(appName):
        print(f"\nRemoving old Mac app {appName!r}")
        subprocess.call(["rm","-rf",appName])
    subprocess.call(["mkdir","-p",appName])
    subprocess.call(["tar","xzvf",tarLoc,'-C',appName])
    # create a link named GSAS-II.py to the script
    if os.path.exists(g2Name): # cleanup
        print(f"\nRemoving sym link {g2Name!r}")
        os.remove(g2Name)
    os.symlink(G2script,g2Name)
    if pythonLoc != os.path.join(installLoc,'../bin',"python"):
        link = os.path.join(appName,'Contents','MacOS','GSAS-II')
        try:
            os.remove(link)
            print(f"\nRemoved sym link {link!r}")
        except FileNotFoundError:
            pass
        print(f"\nOverriding {link!r} with Python location {pythonLoc!r}")
        os.symlink(pythonLoc,link)
    print(f"\nCreated app {appName!r} and {g2Name!r}" +
    "\nViewing app in Finder so you can drag it to the dock if, you wish.")
    subprocess.call(["open","-R",appName])
    
    sys.exit()

#if __name__ == '__main__' and sys.platform == "darwin":
if False:
    iconfile = os.path.join(path2GSAS,'icons','gsas2.icns') # optional icon file
    if not os.path.exists(iconfile): # patch 3/2024 for svn dir organization
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
    os.symlink(os.path.split(G2script)[1],newScript)
    G2script=newScript

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
        f.write(AppleScript.format(newpython,G2script,'',newpython,G2script,''))
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
