#!/usr/bin/env python
'''
*makeLinux: Create Linux Shortcuts*
===================================

This script creates a menu entry for Gnome & KDE desktop managers
and puts it in a place where it can be "indexed."

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
Exec=xterm -hold -e bash -c "{} {}"
Name=GSAS-II
Icon={}
Categories=Science;
"""
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" [<GSAS-II script>]\n")
    sys.exit()
if __name__ == '__main__':
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
    loc = os.path.expanduser('~/Desktop/')
    if "XDG_DATA_HOME" in os.environ and os.path.exists(os.path.join(os.environ.get("XDG_DATA_HOME"),"applications")):
        loc = os.path.join(os.environ.get("XDG_DATA_HOME"),"applications")
    elif "HOME" in os.environ and os.path.exists(os.path.join(os.environ.get("HOME"),".local/share/applications")):
        loc = os.path.join(os.environ.get("HOME"),".local/share/applications")
    elif not os.path.exists(loc):
        loc = os.path.expanduser('~/')
    dfile = os.path.join(loc,'GSASII.desktop')
    os.chmod(script, # make the .py file executable and readable
             stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
    try:
        open(dfile,'w').write(desktop_template.format(sys.executable,script,icon))
        os.chmod(
            dfile,
            stat.S_IWUSR | stat.S_IXUSR | stat.S_IRUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH)
        print("created GNOME/KDE desktop shortcut "+dfile)
    except:
        print("creation of file failed: "+dfile)
