#!/usr/bin/env python
'''
*makeLinux: Create Linux Shortcuts*
===================================

This script creates a menu entry and dektop shortcut for Gnome
(and perhaps KDE) desktop managers. Recent testing on Raspbian.

This is a work in progress as I learn more about shortcuts in Linux.

Run this script with one optional argument, the path to the GSASII.py
The script path may be specified relative to the current path or given
an absolute path, but will be accessed via an absolute path. 
If no arguments are supplied, the GSASII.py script is assumed to be in the
same directory as this file.

'''
import sys, os, os.path, stat, shutil, subprocess, plistlib
desktop_template = """
[Desktop Entry]
Version=1.0
Type=Application
Terminal=false
Exec={} bash -c "{} {}"
Name=GSAS-II
Icon={}
Categories=Science;
"""
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>]\n")
    sys.exit()

if __name__ == '__main__' and sys.platform.startswith('linux'):
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
    icon = os.path.join(os.path.split(script)[0], 'gsas2.png')

    # make sure we found it
    if not os.path.exists(script):
        print("\nFile "+script+" not found")
        Usage()
    if not os.path.exists(icon):
        print("\nWarning: File "+icon+" not found")
        #Usage()
    if os.path.splitext(script)[1].lower() != '.py':
        print("\nScript "+script+" does not have extension .py")
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
    os.chmod(script, # make the .py file executable and readable
             stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
    import shutil
    for term in ("lxterminal", "gnome-terminal", "xterm"):
        try:
            found = shutil.which(term)
        except AttributeError:
            print(' shutil.which() failed (Python 2.7?); assuming',term)
            found = True
        if not found: continue
        if term == "gnome-terminal":
            terminal = 'gnome-terminal -t "GSAS-II console" --'
            script += "; echo Press Enter to close window; read line"
            break
        elif term == "lxterminal":
            terminal = 'lxterminal -t "GSAS-II console" -e'
            script += "; echo Press Enter to close window; read line"
            break
        elif term == "xterm":
            terminal = 'xterm -title "GSAS-II console" -hold -e'
            break
        else:
            print("unknown terminal",term)
            sys.exit()
    else:
        print("No terminal found -- no shortcuts created")
        sys.exit()
    for f,t in zip((dfile,mfile),('Desktop','Menu')):
        if f is None: continue
        try:
            open(f,'w').write(desktop_template.format(
                terminal,
                sys.executable,script,icon))
            os.chmod(
                dfile,
                stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
            print("Created {} shortcut calling {} as file\n\t{}".format(
                t,term,f))
        except Exception as msg:
            print("creation of file failed: "+f)
            print(msg)
