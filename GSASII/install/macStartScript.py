# create a MacOS applet to run GSAS-II. The purpose of the Applet is
# so that the Apple menus are named GSAS-II rather than Python. It also
# allows .gpx files to be dropped on the applet.

# this runs but has not been fully tested

import sys
import os
import subprocess
import shutil
import plistlib
import platform

makePath = os.path.dirname(__file__)  # location of this script
path2GSAS = os.path.dirname(makePath)
appPath = '/tmp/GSAS-II.app'
projectname = 'GSAS-II'
#G2script = os.path.abspath(os.path.join(path2GSAS,"G2.py"))

if __name__ == '__main__' and sys.platform == "darwin":
    iconfile = os.path.join(path2GSAS,'icons','gsas2.icns') # optional icon file
    if not os.path.exists(iconfile): # patch 3/2024 for svn dir organization
        iconfile = os.path.join(path2GSAS,'gsas2.icns') # optional icon file

#     AppleScript = '''(*   GSAS-II AppleScript by B. Toby (brian.toby@anl.gov)
#      It can launch GSAS-II by double clicking or by dropping a data file
#      or folder over the app.
#      It runs GSAS-II in a terminal window.
# *)

# (* test if a file is present and exit with an error message if it is not  *)
# on TestFilePresent(appwithpath)
# 	tell application "System Events"
# 		if (file appwithpath exists) then
# 		else
# 			display dialog "Error: file " & appwithpath & " not found. If you have moved this file recreate the AppleScript with bootstrap.py." with icon caution buttons {{"Quit"}}
# 			return
# 		end if
# 	end tell
# end TestFilePresent

# (*
# ------------------------------------------------------------------------
# this section responds to a double-click. No file is supplied to GSAS-II
# ------------------------------------------------------------------------
# *)
# on run
# 	set python to "{:s}"
# 	set appwithpath to "{:s}"
# 	set env to "{:s}"
# 	TestFilePresent(appwithpath)
# 	TestFilePresent(python)
# 	tell application "Terminal"
# 		do script env & python & " " & appwithpath & "; exit"
# 	end tell
# end run

# (*
# -----------------------------------------------------------------------------------------------
# this section handles starting with files dragged into the AppleScript
#  o it goes through the list of file(s) dragged in
#  o then it converts the colon-delimited macintosh file location to a POSIX filename
#  o for every non-directory file dragged into the icon, it starts GSAS-II, passing the file name
# ------------------------------------------------------------------------------------------------
# *)

# on open names
# 	set python to "{:s}"
# 	set appwithpath to "{:s}"
# 	set env to "{:s}"

# 	TestFilePresent(appwithpath)
# 	repeat with filename in names
# 		set filestr to (filename as string)
# 		if filestr ends with ":" then
#                         (* should not happen, skip directories *)
# 		else
# 			(* if this is an input file, open it *)
# 			set filename to the quoted form of the POSIX path of filename
# 			tell application "Terminal"
# 				activate
# 				do script env & python & " " & appwithpath & " " & filename & "; exit"
# 			end tell
# 		end if
# 	end repeat
# end open
# '''
    AppleScript = '''(*   GSAS-II AppleScript by B. Toby (brian.toby@anl.gov)
     It can launch GSAS-II by double clicking or by dropping a data file
     or folder over the app.
     It runs GSAS-II in a terminal window.
	 
	 This is intended for use with a conda-based GSAS-II distribution where
	 Python is linked from a file (./Contents/MacOS/GSAS-II) inside the current app, 
	 and where the GSAS-II .py files are located in the same directory as this
	 script. This can be renamed but cannot be moved.
*)

on GetPythonLocation()
	(* find python in this script's bundle *)
	tell application "Finder"
		set scriptpath to the POSIX path of (path to me)
	end tell
	set python to (scriptpath & "Contents/MacOS/GSAS-II")
	TestFilePresent(python)
	return python
end GetPythonLocation

on GetG2Location()
	(* find the GSAS-II.py script in the same directory as this script *)
	tell application "Finder"
		set scriptdir to the POSIX path of (container of (path to me) as alias)
	end tell
	set g2script to (scriptdir & "GSAS-II.py")
	TestFilePresent(g2script)
	return g2script
end GetG2Location

on TestFilePresent(filepath)
	(* test if a file is present and exit with an error message if it is not  *)
	tell application "System Events"
		if (file filepath exists) then
		else
			display dialog "Error: file " & filepath & " not found. Was this app moved from the GSASII directory?" with icon caution buttons {"Quit"}
			error number -128
		end if
	end tell
end TestFilePresent

(* 
----------------------------------------------------------------------------
this section responds to a double-click. No file is supplied to GSAS-II
---------------------------------------------------------------------------- 
*)
on run
	set python to GetPythonLocation()
	set appwithpath to GetG2Location()
	
	tell application "Terminal"
		activate
		do script (quoted form of python) & " " & (quoted form of appwithpath) & "; exit"
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
	set python to GetPythonLocation()
	set appwithpath to GetG2Location()
	
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
'''
    # # create a link named GSAS-II.py to the script
    # newScript = os.path.join(path2GSAS,'GSAS-II.py')
    # if os.path.exists(newScript): # cleanup
    #     print("\nRemoving sym link",newScript)
    #     os.remove(newScript)
    # os.symlink(os.path.split(G2script)[1],newScript)
    # G2script=newScript

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
        #f.write(AppleScript.format(newpython,G2script,'',newpython,G2script,''))
        f.write(AppleScript)
        f.close()

        try:
            subprocess.check_output(["osacompile","-o",appPath,shell],stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as msg:
            print('''Error compiling AppleScript.
            Report the next message along with details about your Mac to toby@anl.gov''')
            print(msg.output)
            sys.exit()
        # create a link to conda python relative to this app, named to match the project
        os.symlink('../../../../bin/python',newpython)

        # # test if newpython can run wx
        # def RunPython(image,cmd):
        #     'Run a command in a python image'
        #     try:
        #         err=None
        #         p = subprocess.Popen([image,'-c',cmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #         out = p.stdout.read()
        #         err = p.stderr.read()
        #         p.communicate()
        #         return out,err
        #     except Exception(err):
        #         return '','Exception = '+err

        # testout,errout = RunPython(newpython,'import numpy; import wx; wx.App(); print("-"+"OK-")')
        # if isinstance(testout,bytes): testout = testout.decode()
        # if "-OK-" in testout:
        #     print('wxpython app ran',testout)
        #     break
        # elif i == 1:
        #     print('Run of wx in',pythonExe,'failed, looking for pythonw')
        #     pythonExe = os.path.join(os.path.split(sys.executable)[0],'pythonw')
        #     if not os.path.exists(pythonExe):
        #         print('Warning no pythonw found with ',sys.executable,
        #               '\ncontinuing, hoping for the best')
        # elif i == 2:
        #     print('Warning could not run wx with',pythonExe,
        #               'will try with that external to app')
        #     newpython = pythonExe
        # else:
        #     print('Warning still could not run wx with',pythonExe,
        #               '\ncontinuing, hoping for the best')

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
