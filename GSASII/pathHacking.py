# modules for use where GSAS-II binaries are not co-located with
# the main GSAS-II files and the path is modified (I can hear Tom
# saying "Shame, Shame!").

import glob
import os
import sys
import numpy as np
from . import GSASIIpath

def _path_discovery(printInfo=False):
    def appendIfExists(searchpathlist,loc,subdir):
        newpath = os.path.join(loc,subdir)
        if os.path.exists(newpath):
            if newpath in searchpathlist: return False
            searchpathlist.append(newpath)
    inpver = GSASIIpath.intver(np.__version__)
    
    if GSASIIpath.path2GSAS2 not in sys.path:
        sys.path.insert(0,GSASIIpath.path2GSAS2)  # make sure current path is used
    binpath = None
    binprfx = GSASIIpath.GetBinaryPrefix()
    # places to look for the GSAS-II binary directory
    binseapath = [os.path.abspath(sys.path[0])]  # where Python is installed
    binseapath += [os.path.abspath(os.path.dirname(__file__))]  # directory where this file is found
    binseapath += [os.path.dirname(binseapath[-1])]  # parent of above directory
    binseapath += [os.path.expanduser('~/.GSASII')]       # directory in user's home
    searched = []
    for loc in binseapath:
        if loc in searched: continue
        if not os.path.exists(loc): continue
        searched.append(loc)
        # Look at bin directory (created by a local compile) before looking for standard dist files
        searchpathlist = []
        appendIfExists(searchpathlist,loc,'bin')
        appendIfExists(searchpathlist,loc,'bindist')
        appendIfExists(searchpathlist,loc,'GSASII-bin')
        # also look for directories named by platform etc in loc/AllBinaries or loc
        versions = {}
        namedpath =  glob.glob(os.path.join(loc,'AllBinaries',binprfx+'*'))
        namedpath += glob.glob(os.path.join(loc,'GSASII-bin',binprfx+'*'))
        for d in namedpath:
            d = os.path.realpath(d)
            v = GSASIIpath.intver(d.rstrip('/').split('_')[-1].lstrip('n'))
            versions[v] = d
        vmin = None
        vmax = None
        # try to order the search in a way that makes sense
        for v in sorted(versions.keys()):
            if v <= inpver:
                vmin = v
            elif v > inpver:
                vmax = v
                break
        if vmin in versions and versions[vmin] not in searchpathlist:
            searchpathlist.append(versions[vmin])
        if vmax in versions and versions[vmax] not in searchpathlist:
            searchpathlist.append(versions[vmax])
        for fpth in searchpathlist:
            if not glob.glob(os.path.join(fpth,'pyspg.*')): continue
            if GSASIIpath.pathhack_TestSPG(fpth):
                binpath = fpth   # got one that works, look no farther!
                break
        else:
            continue
        break
    if binpath:  # were GSAS-II binaries found?
        GSASIIpath.binaryPath = binpath
        GSASIIpath.BinaryPathLoaded = True
    else:
        print('*** ERROR: Unable to find GSAS-II binaries. Much of GSAS-II cannot function')
        if GSASIIpath.GetConfigValue('debug'):
            print('Searched directories:')
            for i in searched: print(f'\t{i}')
            print(f'''for subdirectories named .../bin, .../bindist, .../GSASII-bin, 
   .../AllBinaries/{binprfx}* and .../GSASII-bin/{binprfx}*''')
        return True

    # add the data import and export directory to the search path
    if binpath not in sys.path: sys.path.insert(0,binpath)
    if printInfo: print(f'GSAS-II binary directory: {binpath}')
        
    # *** Thanks to work by Tom, imports and exports are now found directly
    # *** and the code below is no longer needed.
    # ***
    #newpath = os.path.join(path2GSAS2,'imports')
    #if newpath not in sys.path: sys.path.append(newpath)
    #newpath = os.path.join(path2GSAS2,'exports')
    #if newpath not in sys.path: sys.path.append(newpath)
