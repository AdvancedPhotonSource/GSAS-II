#!/usr/bin/env python
'''
This script creates a menu entry and dektop shortcut for Gnome
(and perhaps KDE) desktop managers. The most recent testing 
has been on Raspberry Pi OS.
My hope is to improve this further to work conveniently with a wider 
range of Linux desktop managers. 

Run this script with one optional argument, the location of the GSASII.py
file. That location may be specified relative to the current path or given
an absolute path, but will be accessed via an absolute path. 
If no arguments are supplied, the path to the ``GSASII.py`` file
is assumed to be in the parent directory to the one where this file 
(``makeLinux.py``) is found. 

The contents of this file may also be run from inside the gitstrap.py
installation script. In that case, the following variables are already 
defined: 

  * path2GSAS2   is the directory with all GSAS-II Python code
  * G2script     has the location of the GSASII.py file
  * path2repo    is the location of the GSAS-II git repository

The path to Python is determined from the version of Python used to 
run this script.
'''
import sys, os, os.path, stat, shutil, subprocess, plistlib
import datetime
desktop_template = """
[Desktop Entry]
Version=1.0
Type=Application
Terminal=false
Exec={} bash -c "{}"
Name=GSAS-II
Icon={}
Categories=Science;
"""
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>]\n")
    sys.exit()

if __name__ == '__main__' and sys.platform.startswith('linux'):
    if __file__.endswith("makeLinux.py") and len(sys.argv) == 1:
        # find the main GSAS-II script if not on command line
        path2GSAS2 = os.path.dirname(os.path.dirname(__file__))
        path2repo = os.path.dirname(path2GSAS2)
        G2script = os.path.abspath(path2GSAS2,"GSASII.py")
    elif __file__.endswith("makeLinux.py") and len(sys.argv) == 2:
        G2script = os.path.abspath(sys.argv[1])
        path2GSAS2 = os.path.dirname(G2script)
        path2repo = os.path.dirname(path2GSAS2)
    elif __file__.endswith("makeLinux.py"):
        Usage()
    else:
        print(f"running makeLinux.py indirectly inside {__file__!r}")
    pythonexe = os.path.realpath(sys.executable)
    G2icon = os.path.join(path2GSAS2,'icons','gsas2.png')
    
    print('Python installed at ',pythonexe)
    print('GSAS-II installed at',path2GSAS2)
    print('GSASII.py at        ',G2script)
    print('GSASII icon at      ',G2icon)

    # make sure we have the stuff we need
    if not os.path.exists(G2script):
        print("\nFile "+G2script+" not found")
        Usage()
    if not os.path.exists(G2icon):
        print("\nWarning: File "+G2icon+" not found")
        #Usage()
    if os.path.splitext(G2script)[1] != '.py':
        print("\nScript "+G2script+" does not have extension .py")
        Usage()
        
    mfile = None   # menu file
    if "XDG_DATA_HOME" in os.environ and os.path.exists(
            os.path.join(os.environ.get("XDG_DATA_HOME"),"applications")):
        mfile = os.path.join(os.environ.get("XDG_DATA_HOME"),"applications",
                                 'GSASII.desktop')
    elif "HOME" in os.environ and os.path.exists(
            os.path.join(os.environ.get("HOME"),".local/share/applications")):
        mfile = os.path.join(os.environ.get("HOME"),".local/share/applications",
                                 'GSASII.desktop')
    dfile = None   # desktop file
    if os.path.exists(os.path.expanduser('~/Desktop/')):
        dfile = os.path.join(os.path.expanduser('~/Desktop'),'GSASII.desktop')
    os.chmod(G2script, # make the .py file executable and readable
             stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)

    # create a script that activates conda Python and then starts GSAS-II
    # that would be run from a terminal
    G2start = os.path.abspath(os.path.normpath(os.path.join(path2repo,'..','RunGSASII.sh')))
    if os.path.exists(G2start): os.unlink(G2start)
    fp = open(G2start,'w')
    fp.write('#!/bin/bash\n')
    fp.write("# created by run of makeLinux.py on {:%d %b %Y %H:%M}\n".format(
        datetime.datetime.now()))
    activate = os.path.join(os.path.dirname(pythonexe),'activate')
    if os.path.exists(activate): fp.write(f'source {activate}\n')
    fp.write(f'{pythonexe} {G2script} $*\n')
    fp.close()
    os.chmod(G2start, # make the .py file executable and readable
             stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)

    script = ''
    import shutil
    for term in ("lxterminal", "gnome-terminal", 'konsole', "xterm",
                         "terminator", "terminology", "tilix"):
        try:
            found = shutil.which(term)
            if not found: continue
        except AttributeError:
            print(f'shutil.which() failed (why?); assuming {term} present')
            found = True
        if term == "gnome-terminal":
            terminal = 'gnome-terminal -t "GSAS-II console" --'
            script = "echo Press Enter to close window; read line"
            break
        elif term == "lxterminal":
            terminal = 'lxterminal -t "GSAS-II console" -e'
            script = "echo Press Enter to close window; read line"
            break
        elif term == "xterm":
            terminal = 'xterm -title "GSAS-II console" -hold -e'
            break
        elif term == "terminator":
            terminal = term + ' -T "GSAS-II console" -x'
            script = "echo;echo Press Enter to close window; read line"
            break
        elif term == "konsole":
            terminal = term + ' -p tabtitle="GSAS-II console" --hold -e'
            script = "echo; echo This window can now be closed"
            break
        elif term == "tilix":
            terminal = term + ' -t "GSAS-II console" -e'
            script = "echo;echo Press Enter to close window; read line"
            break
        elif term == "terminology":
            terminal = term + ' -T="GSAS-II console" --hold -e'
            script = "echo; echo This window can now be closed"
            break                
        else:
            print("unknown terminal",term)
            sys.exit()
    else:
        print("No terminal found -- no shortcuts created")
        sys.exit()
    add2script = ''
    if script: add2script = '; ' + script 
    for f,t in zip((dfile,mfile),('Desktop','Menu')):
        if f is None: continue
        try:
            open(f,'w').write(desktop_template.format(
                terminal,
                G2start+add2script,G2icon))
            os.chmod(
                dfile,
                stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
            print("Created {} shortcut calling {} as file\n\t{}".format(
                t,term,f))
        except Exception as msg:
            print("creation of file failed: "+f)
            print(msg)

    # now create a script that opens GSAS-II in a window
    G2startterm = os.path.normpath(os.path.join(path2repo,'..','RunG2inTerm.sh'))
    if os.path.exists(G2startterm): os.unlink(G2startterm)
    fp = open(G2startterm,'w')
    fp.write('#!/bin/bash\n')
    fp.write("# created by run of makeLinux.py on {:%d %b %Y %H:%M}\n".format(
        datetime.datetime.now()))
    if os.path.exists(activate): fp.write(f'source {activate}\n')
    fp.write(f'{terminal} {G2start}\n')
    if os.path.exists(script): fp.write(f'{script}\n')
    fp.close()
    os.chmod(G2startterm, stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)

    # Finally, create a script to reset GSAS-II, dropping any local
    # changes and downloading the latest version
    gitstrap = os.path.abspath(
        os.path.normpath(os.path.join(path2repo,'..','gitstrap.py')))
    if not os.path.exists:
        print(f'the installation script was not found: {gitstrap!r}')
    else:
        G2reset = os.path.normpath(os.path.join(path2repo,'..','Reset2FreshGSASII.sh'))
        if os.path.exists(G2reset): os.unlink(G2reset)
        fp = open(G2reset,'w')
        fp.write('#!/bin/bash\n')
        fp.write("# created by run of makeLinux.py on {:%d %b %Y %H:%M}\n".format(
            datetime.datetime.now()))
        if os.path.exists(activate): fp.write(f'source {activate}\n')
        fp.write('read -p "Reset any local changes and install latest GSAS-II version? (Y/[N]): " confirm\n[[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1\n')
        fp.write(f'{pythonexe} {gitstrap} --reset\n')
        fp.close()
        os.chmod(G2reset, stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
        print(f'Created {G2reset!r} to reset GSAS-II installation when all else fails...')
