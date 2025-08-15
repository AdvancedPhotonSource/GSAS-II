#!/usr/bin/env python
#------------------------------------------------------------
# this version intended for use with git installations
#------------------------------------------------------------
# TODO: clean up use of args: installLoc will always be the parent of path2GSAS
'''
This routine creates an app bundle named GSAS-II.app. Inside the
bundle is a symbolic link to the Python executable named "GSAS-II"
that will be used to run GSAS-II. Having this link named that
way causes the name of the app to shows in the menu bar as
"GSAS-II" rather than "Python". Also used by the app, is another
symbolic link named GSAS-II.py, which must be placed in the same
directory as the app bundle. This file is linked to the G2.py
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

    path2GSAS:  The location where the GSAS-II Python files are found.
    installLoc: The location where the GSAS-II.app app bundle and
       the GSAS-II.py will be placed. This will be the parent of path2GSAS
    pythonLoc:  The location of the Python executable.

Under normal circumstances, the locations for all of these paths
can be determined from the location of the makeMacApp.py file.
Note that when GSAS-II is installed from git using gitstrap.py,
the git repository is placed at <loc>/GSAS-II and the GSAS-II
Python scripts are placed in the <loc>/GSAS-II/GSASII child directory. 
GSAS-II is started from the <loc>/GSAS-II/GSAS-II.py script created here
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
location G2.py script, which can be a relative path
(the absolute path is determined). If not supplied, the G2.py script
is expected to be located in the directory above where this
(makeMacApp.py) script is found.

The third argument, if supplied, provides the full path for the Python
installation to be used inside the app bundle that will be created. If not
supplied, and Python exists at installLoc/../bin/python, that will be used.
If that does not exist, then the location of the current Python executable
(from sys.executable) will be used.

'''

import os
import os.path
import subprocess
import sys


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
    # G2script = os.path.abspath(os.path.join(path2GSAS,"G2.py"))
    # if not os.path.exists(G2script):
    #     print(f"\nERROR: File {G2script!r} not found")
    #     Usage()
    # if os.path.splitext(G2script)[1].lower() != '.py':
    #     print(f"\nScript {G2script!r} does not have extension .py")
    #     Usage()
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
    #print(f'Using GSAS-II script: {G2script}')
    print(f'Install location: {installLoc}')

    # files to be created
    appName = os.path.abspath(os.path.join(installLoc,project+".app"))
    g2Name = os.path.abspath(os.path.join(installLoc,project+'.py'))

# new approach, use previously created tar (.tgz) file
# with pre-built startup app. See macStartScript.py
if __name__ == '__main__' and sys.platform == "darwin":
    if os.path.exists(appName):
        print(f"\nRemoving old Mac app {appName!r}")
        subprocess.call(["rm","-rf",appName])
    subprocess.call(["mkdir","-p",appName])
    subprocess.call(["tar","xzvf",tarLoc,'-C',appName])
    # create a script named GSAS-II.py to be run by the AppleScript
    if os.path.exists(g2Name): # cleanup
        print(f"\nRemoving {g2Name!r}")
        os.remove(g2Name)
    #os.symlink(G2script,g2Name)
    with open(g2Name,'w') as fp:
        fp.write('''# this script starts GSAS-II when not installed into Python
# it will be called from the GSAS-II.app AppleScript
# it should be not be renamed or moved
import sys,os
print('Hacking sys.path to provide access to GSAS-II')
sys.path.insert(0,os.path.dirname(__file__))
from GSASII.GSASIIGUI import main
main()
''')
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
