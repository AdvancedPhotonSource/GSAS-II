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
    import platform
    if platform.architecture()[0] == '64bit':
        bindir = 'binmac64-%d.%d' % sys.version_info[0:2]
    else:
        bindir = 'binmac%d.%d' % sys.version_info[0:2]
    if platform.mac_ver()[0].startswith('10.5.'):
        bindir += '_10.5'
elif sys.platform == "linux2":
    if platform.architecture()[0] == '64bit':
        bindir = 'binlinux64-%d.%d' % sys.version_info[0:2]
    else:
        bindir = 'binlinux%d.%d' % sys.version_info[0:2]
for loc in sys.path[0],os.path.abspath(os.path.split(__file__)[0]):
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
# add the data import and export directory to the search path
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # location of this file; save before any changes in pwd
newpath = os.path.join(path2GSAS2,'imports')
if newpath not in sys.path: sys.path.append(newpath)
newpath = os.path.join(path2GSAS2,'exports')
if newpath not in sys.path: sys.path.append(newpath)

# setup read of config.py, if present
try:
    import config
    configDict = config.__dict__
    import inspect
    vals = [True for i in inspect.getmembers(config) if '__' not in i[0]]
    print str(len(vals))+' values read from config file '+os.path.abspath(config.__file__)
except ImportError:
    configDict = {}
    
def GetConfigValue(key,default=None):
    '''Return the configuration file value for key or a default value if not present
    
    :param str key: a value to be found in the configuration (config.py) file
    :param default: a value to be supplied is none is in the config file or
      the config file is not found. Defaults to None
    :returns: the value found or the default.
    '''
    return configDict.get(key,default)

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

def LoadConfigFile(filename):
    '''Read a GSAS-II configuration file.
    Comments (starting with "%") are removed, as are empty lines
    
    :param str filename: base file name (such as 'file.dat'). Files with this name
      are located from the path and the contents of each are concatenated.
    :returns: a list containing each non-empty (after removal of comments) line
      found in every matching config file.
    '''
    info = []
    for path in sys.path:
        fil = os.path.join(path,filename)
        if not os.path.exists(fil): continue
        try:
            i = 0
            fp = open(fil,'r')
            for line in fp:
                expr = line.split('#')[0].strip()
                if expr:
                    info.append(expr)
                    i += 1
            print(str(i)+' lines read from config file '+fil)
        finally:
            fp.close()
    return info


# routines to interface with subversion
def whichsvn():
    '''Returns a path to the subversion exe file, if any is found.
    Searches the current path as well as subdirectory "svn" and
    "svn/bin" in the location of the GSASII source files.

    :returns: None if svn is not found or an absolute path to the subversion
      executable file.
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
            return os.path.abspath(exe_file)

def svnGetLog(fpath=os.path.split(__file__)[0],version=None):
    '''Get the revision log information for a specific version of the 

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located.
    :param int version: the version number to be looked up or None (default)
       for the latest version.

    :returns: a dictionary with keys (one hopes) 'author', 'date', 'msg', and 'revision'

    '''
    import subprocess
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    if version is not None:
        vstr = '-r'+str(version)
    else:
        vstr = '-rHEAD'

    cmd = [svn,'log',fpath,'--xml',vstr]
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err:
        print 'out=',out
        print 'err=',err
        return None
    x = ET.fromstring(out)
    d = {}
    for i in x.iter('logentry'):
        d = {'revision':i.attrib.get('revision','?')}
        for j in i:
            d[j.tag] = j.text
        break # only need the first
    return d

def svnGetRev(fpath=os.path.split(__file__)[0],local=True):
    '''Obtain the version number for the either the last update of the local version
    or contacts the subversion server to get the latest update version (# of Head).

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    :param bool local: determines the type of version number, where
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
    s = subprocess.Popen(cmd+['--non-interactive', '--trust-server-cert'],
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err:
        print 'svn failed, retry w/o --trust...\nout=',out
        print 'err=',err
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
    cmd = [svn,'status',fpath,'--xml']
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err: return None
    x = ET.fromstring(out)
    changed = []
    for i in x.iter('entry'):
        if i.find('wc-status').attrib.get('item') == 'modified': 
            changed.append(i.attrib.get('path'))
    return changed

def svnUpdateDir(fpath=os.path.split(__file__)[0],version=None):
    '''This performs an update of the files in a local directory from a server. 

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    :param version: the number of the version to be loaded. Used only
       cast as a string, but should be an integer or something that corresponds to a
       string representation of an integer value when cast. A value of None (default)
       causes the latest version on the server to be used.
    '''
    import subprocess
    svn = whichsvn()
    if not svn: return
    if version:
        verstr = '-r' + str(version)
    else:
        verstr = '-rHEAD'
    cmd = [svn,'update',fpath,verstr,
           '--non-interactive',
           '--accept','theirs-conflict','--force']
    s = subprocess.Popen(cmd, 
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = s.communicate()
    print out
    if err:
        print(60*"=")
        print ("****** An error was noted, see below *********")
        print(60*"=")
        print err
        sys.exit()

def svnUpdateProcess(version=None,projectfile=None):
    '''perform an update of GSAS-II in a separate python process'''
    import subprocess
    if not projectfile:
        projectfile = ''
    else:
        projectfile = os.path.realpath(projectfile)
        print 'restart using',projectfile
    if not version:
        version = ''
    else:
        version = str(version)
    # start the upgrade in a separate interpreter (avoids loading .pyd files)
    subprocess.Popen([sys.executable,__file__,projectfile,version])
    sys.exit()

if __name__ == '__main__':
    import subprocess
    import time
    time.sleep(1) # delay to give the main process a chance to exit
    # perform an update and restart GSAS-II
    project,version = sys.argv[1:3]
    loc = os.path.dirname(__file__)
    if version:
        print("Regress to version "+str(version))
        svnUpdateDir(loc,version=version)
    else:
        print("Update to current version")
        svnUpdateDir(loc)
    if project:
        print("Restart GSAS-II with project file "+str(project))
        subprocess.Popen([sys.executable,os.path.join(loc,'GSASII.py'),project])
    else:
        print("Restart GSAS-II without a project file ")
        subprocess.Popen([sys.executable,os.path.join(loc,'GSASII.py')])
    print 'exiting update process'
    sys.exit()
    
