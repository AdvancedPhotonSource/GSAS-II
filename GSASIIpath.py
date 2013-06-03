# -*- coding: utf-8 -*-
'''
*GSASIIpath: locations & updates*
---------------------------------

Routines for dealing with file locations, etc.

Determines the location of the compiled (.pyd or .so) libraries.

Interfaces with subversion (svn): 
Determine the subversion release number by determining the highest version number
where :func:`SetVersionNumber` is called (best done in every GSASII file).
Other routines will update GSASII from the subversion server if svn can be
found.
'''

import os
import sys
import platform
# determine a binary path for the pyd files based on the host OS and the python version,  
# path is relative to location of the script that is called as well as this file
# this must be imported before anything that imports any .pyd/.so file for GSASII
bindir = None
if sys.platform == "win32":
    if platform.architecture()[0] == '64bit':
        bindir = 'binwin64-%d.%d' % sys.version_info[0:2]
    else:
        bindir = 'binwin%d.%d' % sys.version_info[0:2]
elif sys.platform == "darwin":
    bindir = 'binmac%d.%d' % sys.version_info[0:2]
    import platform
    if platform.mac_ver()[0].startswith('10.5.'):
        bindir += '_10.5'
elif sys.platform == "linux2":
    if platform.architecture()[0] == '64bit':
        bindir = 'binlinux64-%d.%d' % sys.version_info[0:2]
    else:
        bindir = 'binlinux%d.%d' % sys.version_info[0:2]
for loc in sys.path[0],os.path.split(__file__)[0]:
    if bindir:
        if os.path.exists(os.path.join(loc,bindir)) and os.path.join(loc,bindir) not in sys.path: 
            sys.path.insert(0,os.path.join(loc,bindir))
        # is there a bin directory? (created by a local compile), if so put
        # that at the top of the path
    if os.path.exists(os.path.join(loc,'bin')) and os.path.getsize(os.path.join(loc,'bin')):
        bindir = 'bin'
        if os.path.join(loc,'bin') not in sys.path: 
            sys.path.insert(0,os.path.join(loc,bindir))
print 'GSAS-II binary directory: ',os.path.join(loc,bindir)
if bindir == None:
    raise Exception,"**** ERROR GSAS-II binary libraries not found, GSAS-II fails ****"

# routines for looking a version numbers in files
version = -1
def SetVersionNumber(RevString):
    '''Set the subversion version number

    :param str RevString: something like "$Revision$"
      that is set by subversion when the file is retrieved from subversion.

    Place ``GSASIIpath.SetVersionNumber("$Revision$")`` in every python
    file.
    '''
    try:
        RevVersion = int(RevString.split(':')[1].split()[0])
        global version
        version = max(version,RevVersion)
    except:
        pass
        
def GetVersionNumber():
    '''Return the maximum version number seen in :func:`SetVersionNumber`
    '''
    return version

# routines to interface with subversion
def whichsvn():
    '''Returns a path to the subversion exe file, if any is found.
    Searches the current path as well as subdirectory "svn" and
    "svn/bin" in the location of the GSASII source files.

    :returns: None if svn is not found.
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    svnprog = 'svn'
    if sys.platform == "win32": svnprog += '.exe'
    pathlist = os.environ["PATH"].split(os.pathsep)
    pathlist.insert(0,os.path.join(os.path.split(__file__)[0],'svn'))
    pathlist.insert(1,os.path.join(os.path.split(__file__)[0],'svn','bin'))
    for path in pathlist:
        exe_file = os.path.join(path, svnprog)
        if is_exe(exe_file):
            return exe_file

def svnGetRev(fpath=os.path.split(__file__)[0],local=True):
    '''This obtains the version number for the either the latest local last update
    or contacts the subversion server to get the latest update version (# of Head).

    :param fpath: path to repository dictionary, defaults to directory where
       the current file is located
    :param local: determines the type of version number, where
       True (default): returns the latest installed update 
       False: returns the version number of Head on the server

    :Returns: the version number as an str or 
       None if there is a subversion error (likely because the path is
       not a repository or svn is not found)
    '''

    import subprocess
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    if local:
        cmd = [svn,'info',fpath,'--xml']
    else:
        cmd = [svn,'info',fpath,'--xml','-rHEAD']
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err:
        print 'out=',out
        print 'err=',err
        return None
    x = ET.fromstring(out)
    for i in x.iter('entry'):
        rev = i.attrib.get('revision')
        if rev: return rev

def svnFindLocalChanges(fpath=os.path.split(__file__)[0]):
    '''Returns a list of files that were changed locally. If no files are changed,
       the list has length 0

    :param fpath: path to repository dictionary, defaults to directory where
       the current file is located

    :returns: None if there is a subversion error (likely because the path is
       not a repository or svn is not found)

    '''
    import subprocess
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    s = subprocess.Popen([svn,'status',fpath,'--xml'],
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err: return None
    x = ET.fromstring(out)
    changed = []
    for i in x.iter('entry'):
        if i.find('wc-status').attrib.get('item') == 'modified': 
            changed.append(i.attrib.get('path'))
    return changed

def svnUpdateDir(fpath=os.path.split(__file__)[0]):
    '''This performs an update of the files in a local directory from a server. 

    :param fpath: path to repository dictionary, defaults to directory where
       the current file is located

    :returns: A dictionary with the files that have been changed/added and
          a code describing how they have been updated (see changetype) ro 
          None if there is a subversion error (likely because the path is
          not a repository or svn is not found)

    '''
    import subprocess
    changetype = {'A': 'Added', 'D': 'Deleted', 'U': 'Updated',
                  'C': 'Conflict', 'G': 'Merged', 'E': 'Replaced'}
    svn = whichsvn()
    if not svn: return
    cmd = [svn,'update',fpath,'-rHEAD',
           '--non-interactive',
           '--accept','theirs-conflict','--force']
    s = subprocess.Popen(cmd, 
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err: return
    l = out.split()
    updates = {}
    for i,j in zip(l[::2],l[1::2]):
        if i == 'Updated': break
        if i == 'At': break
        t = changetype.get(i[0])
        if not t: continue
        updates[j] = t
    return updates
