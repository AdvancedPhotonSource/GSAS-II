# -*- coding: utf-8 -*-
#GSASIIpath - file location & update routines
'''
:mod:`GSASIIpath` Classes & routines follow
'''

from __future__ import division, print_function
import os
import sys
import platform
import glob
import subprocess
import datetime as dt
try:
    import numpy as np
except ImportError:
    print("skipping numpy in GSASIIpath")
try:
    import requests
except:
    print('Python requests package not installed (required for web access')

# fix up path before using git. Needed when using conda without
#   activate (happens on MacOS w/GSAS-II.app)
pyPath = os.path.dirname(os.path.realpath(sys.executable))
if sys.platform != "win32" and pyPath not in os.environ['PATH'].split(':'):
    os.environ['PATH'] = pyPath + ':' + os.environ['PATH']
    
try:
    import git
except ImportError as msg:
    if 'Failed to initialize' in msg.msg:
        print('The gitpython package is unable to locate a git installation.')
        print('See https://gsas-ii.readthedocs.io/en/latest/packages.html for more information.')
    elif 'No module' in msg.msg:
        print('Python gitpython module not installed')
    else:
        print(f'gitpython failed to import, but why? Error:\n{msg}')
    print('Note: git and gitpython are required for GSAS-II to self-update')
except Exception as msg:
    print(f'git import failed with unexpected error:\n{msg}')
    print('Note: git and gitpython are required for GSAS-II to self-update')

# hard-coded github repo locations
G2binURL = "https://api.github.com/repos/AdvancedPhotonSource/GSAS-II-buildtools"
g2URL = "https://github.com/AdvancedPhotonSource/GSAS-II.git"
# tutorial repo owner & Repo name
gitTutorialOwn,gitTutorialRepo = 'AdvancedPhotonSource', 'GSAS-II-Tutorials'
    
path2GSAS2 = os.path.dirname(os.path.abspath(os.path.expanduser(__file__))) # location of this file; save before any changes in pwd

# convert version numbers as '1.2.3' to integers (1002) and back (to 1.2)
fmtver = lambda v: str(v//1000)+'.'+str(v%1000)
intver = lambda vs: sum([int(i) for i in vs.split('.')[0:2]]*np.array((1000,1)))

def GetConfigValue(key,default=None):
    '''Return the configuration file value for key or a default value if not present
    
    :param str key: a value to be found in the configuration (config.py) file
    :param default: a value to be supplied if none is in the config file or
      the config file is not found. Defaults to None
    :returns: the value found or the default.
    '''
    try:
        return configDict.get(key,default)
    except NameError: # this happens when building docs
        return None

def SetConfigValue(parmdict):
    '''Set configuration variables from a dictionary where elements are lists
    First item in list is the default value and second is the value to use.
    '''
    global configDict
    for var in parmdict:
        if var in configDict:
            del configDict[var]
        if isinstance(parmdict[var],tuple):
            configDict[var] = parmdict[var]
        else:
            if parmdict[var][1] is None: continue
            if parmdict[var][1] == '': continue
            if parmdict[var][0] == parmdict[var][1]: continue
            configDict[var] = parmdict[var][1]

def addPrevGPX(fil,cDict):
    '''Add a GPX file to the list of previous files. 
    Move previous names to start of list. Keep most recent five files
    '''
    global configDict
    fil = os.path.abspath(os.path.expanduser(fil))
    if 'previous_GPX_files' not in cDict: 
        cDict['previous_GPX_files'] = [[],[],[],'Previous .gpx files'] # unexpected
    try:
        pos = cDict['previous_GPX_files'][1].index(fil) 
        if pos == 0: return
        cDict['previous_GPX_files'][1].pop(pos) # if present, remove if not 1st
    except ValueError:
        pass
    except AttributeError:
        cDict['previous_GPX_files'][1] = []
    files = list(cDict['previous_GPX_files'][1])
    files.insert(0,fil)
    cDict['previous_GPX_files'][1] = files[:5]
    configDict['previous_GPX_files'] = cDict['previous_GPX_files'][1]

# routines for looking a version numbers in files
version = -1
def SetVersionNumber(RevString):
    '''Set the subversion (svn) version number

    :param str RevString: something like "$Revision: 5796 $"
      that is set by subversion when the file is retrieved from subversion.

    Place ``GSASIIpath.SetVersionNumber("$Revision: 5796 $")`` in every python
    file.
    '''
    try:
        RevVersion = int(RevString.split(':')[1].split()[0])
        global version
        version = max(version,RevVersion)
    except:
        pass
        
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
        fil = os.path.join(path,'inputs',filename)
        if not os.path.exists(fil):  # patch 3/2024 for svn dir organization
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

def GetBinaryPrefix(pyver=None):
    '''Creates the first part of the binary directory name
    such as linux_64_p3.9 (where the full name will be 
    linux_64_p3.9_n1.21). 

    Note that any change made here is also needed in GetBinaryDir in 
    fsource/SConstruct
    '''
    if sys.platform == "win32":
        prefix = 'win'
    elif sys.platform == "darwin":
        prefix = 'mac'
    elif sys.platform.startswith("linux"):
        prefix = 'linux'
    else:
        print(u'Unknown platform: '+sys.platform)
        raise Exception('Unknown platform')
    if 'arm' in platform.machine() and sys.platform == "darwin":
        bits = 'arm'
    elif 'aarch' in platform.machine() and '64' in platform.architecture()[0]:
        bits = 'arm64'
    elif 'arm' in platform.machine():
        bits = 'arm32'
    elif '64' in platform.architecture()[0]:
        bits = '64'
    else:
        bits = '32'

    # format current python version
    if pyver:
        pyver = 'p'+pyver
    else:
        pyver = 'p{}.{}'.format(*sys.version_info[0:2])

    return '_'.join([prefix,bits,pyver])

#==============================================================================
#==============================================================================
# hybrid routines that use git & svn (to be revised to remove svn someday)
G2_installed_result = None
def HowIsG2Installed():
    '''Determines if GSAS-II was installed with git, svn or none of the above.
    Result is cached to avoid time needed for multiple calls of this.

    :returns: 
      * a string starting with 'git' from git, 
        if installed from the GSAS-II GitHub repository (defined in g2URL), 
        the string is 'github', if the post-3/2024 directory structure is
        in use '-rev' is appended.
      * or 'svn' is installed from an svn repository (assumed as defined in g2home)
      * or 'noVCS' if installed without a connection to a version control system
    '''
    global G2_installed_result
    if G2_installed_result is not None: return G2_installed_result
    try:
        g2repo = openGitRepo(path2GSAS2)
        if os.path.samefile(os.path.dirname(g2repo.common_dir),path2GSAS2):
            rev = ''
        else:
            rev = '-rev'
        if g2URL in g2repo.remote().urls:
            return 'github'+rev
        G2_installed_result = 'git'+rev
        return G2_installed_result
    except:
        pass
    if svnGetRev(verbose=False):
        G2_installed_result = 'svn'
        return G2_installed_result
    G2_installed_result = 'noVCS'
    return G2_installed_result
        
def GetVersionNumber():
    '''Obtain a version number for GSAS-II from git, svn or from the 
    files themselves, if no other choice. 

    This routine was used to get the GSAS-II version from strings 
    placed in files by svn with the version number being the latest 
    number found, gathered by :func:`SetVersionNumber` (not 100% accurate
    as the latest version might have files changed that are not tagged
    by svn or with a call to SetVersionNumber. Post-svn this info
    will not be reliable, and this mechanism is replaced by a something
    created with a git hook, file git_verinfo.py (see the git_filters.py file).

    Before resorting to the approaches above, try to ask git or svn 
    directly.

    :returns: an int value usually, but a value of 'unknown' might occur 
    '''
    if HowIsG2Installed().startswith('git'):
        g2repo = openGitRepo(path2GSAS2)
        for h in list(g2repo.iter_commits('HEAD'))[:50]: # (don't go too far back)
            tags = g2repo.git.tag('--points-at',h).split('\n')
            try:
                for item in tags:
                    if item.isnumeric(): return int(item)
            except:
                pass
        
    elif HowIsG2Installed() == 'svn':
        rev = svnGetRev()
        if rev is not None: return rev

    # No luck asking, look up version information from git_verinfo.py
    try:
        import git_verinfo as gv
        try:
            for item in gv.git_tags:
                if item.isnumeric(): return int(item)
        except:
            pass
        try:
            for item in gv.git_prevtags:
                if item.isnumeric(): return int(item)
        except:
            pass
    except:
        pass
        
    # all else failed, use the SetVersionNumber value
    if version > 5000:  # a small number must be wrong
        return version
    else:
        return "unknown"
    
def getG2VersionInfo():
    if HowIsG2Installed().startswith('git'):
        g2repo = openGitRepo(path2GSAS2)
        commit = g2repo.head.commit
        ctim = commit.committed_datetime.strftime('%d-%b-%Y %H:%M')
        now = dt.datetime.now().replace(
            tzinfo=commit.committed_datetime.tzinfo)
        delta = now - commit.committed_datetime
        age = delta.total_seconds()/(60*60*24.)
        tags = g2repo.git.tag('--points-at',commit).split('\n')
        tags = [i for i in tags if i.isnumeric()]
        # get sequential version # (tag)
        gversion = ""
        if len(tags) >= 1:
            gversion = f"Tag: #{tags[0]}"
        else:
            for h in list(g2repo.iter_commits(commit))[:50]: # (don't go too far back)
                for t in g2repo.git.tag('--points-at',h).split('\n'):
                    if t.isnumeric():
                        gversion = f"Last tag: #{t}"
                        break
                if gversion: break
            else:
                gversion = "No tag?!" # if all else fails
        msg = ''
        if g2repo.head.is_detached:
            msg = ("\n" +
            "**** You have reverted to a past version of GSAS-II. Please \n"
            +
            "contact the developers with what is preferred in this version ****"
                    )
        else:
            msg = ''
            rc,lc,_ = gitCheckForUpdates(False,g2repo)
            if rc is None:
                msg += f"\n\tNo history found. On development branch? ({g2repo.active_branch})"
            elif str(g2repo.active_branch) != 'master':
                msg += f'\n\tUsing development branch "{g2repo.active_branch}"'
            elif age > 60 and len(rc) > 0:
                msg += f"\n\t**** This version is really old. Please update. >= {len(rc)} updates have been posted ****"
            elif age > 5 and len(rc) > 0:
                msg += f"\n\t**** Please consider updating. >= {len(rc)} updates have been posted"
            elif len(rc) > 0:
                msg += f"\n\tThis GSAS-II version is ~{len(rc)} updates behind current."
        return f"  GSAS-II:    {commit.hexsha[:6]}, {ctim} ({age:.1f} days old). {gversion}{msg}"
    elif HowIsG2Installed() == 'svn':
        rev = svnGetRev()
        if rev is None: 
            "no SVN"
        else:
            rev = f"SVN version {rev}"

        # patch 11/2020: warn if GSASII path has not been updated past v4576.
        # For unknown reasons on Mac with gsas2full, there have been checksum
        # errors in the .so files that prevented svn from completing updates.
        # If GSASIIpath.svnChecksumPatch is not present, then the fix for that
        # has not been retrieved, so warn. Keep for a year or so. 
        try:
            svnChecksumPatch
        except:
            print('Warning GSAS-II incompletely updated. Please contact toby@anl.gov')
        # end patch
            
        return f"Latest GSAS-II revision: {GetVersionNumber()} (svn {rev})"
    else:
        try:
            import git_verinfo as gv
            if gv.git_tags:
                msg = f"{' '.join(gv.git_tags)}, Git: {gv.git_version[:6]}"
            else:
                msg = (f"{gv.git_version[:6]}; "+
                       f"Prev ver: {' '.join(gv.git_prevtags)}"+
                       f", {gv.git_prevtaggedversion[:6]}")
            return f"GSAS-II version: {msg} (Manual update)"
        except:
            pass
    # all else fails, use the old version number routine
    return f"GSAS-II installed manually, last revision: {GetVersionNumber()}"

#==============================================================================
#==============================================================================
# routines to interface with git.

BASE_HEADER = {'Accept': 'application/vnd.github+json',
               'X-GitHub-Api-Version': '2022-11-28'}

def openGitRepo(repo_path):
    try:  # patch 3/2024 for svn dir organization
        return git.Repo(path2GSAS2)
    except git.InvalidGitRepositoryError:
        pass
    return git.Repo(os.path.dirname(path2GSAS2))
        
def gitLookup(repo_path,gittag=None,githash=None):
    '''Return information on a particular checked-in version
    of GSAS-II. 

    :param str repo_path: location where GSAS-II has been installed
    :param str gittag: a tag value. 
    :param str githash: hex hash code (abbreviated to as few characters as 
       needed to keep it unique). If None (default), a tag must be supplied.
    :returns: either None if the tag/hash is not found or a tuple with 
       four values (hash, tag-list, message,date_time) where 

        * hash (str) is the git checking hash code; 
        * tag-list is a list of tags (typically there will 
          be one or two); 
        * message is the check-in message (str)
        * date_time is the check-in date as a datetime object
    '''
    g2repo = openGitRepo(repo_path)
    if gittag is not None and githash is not None:
        raise ValueError("Cannot specify a hash and a tag")
    if gittag is not None:
        try:
            commit = g2repo.tag(gittag).commit
        except ValueError:
            return None
    elif githash is not None:
        try:
            commit = g2repo.commit(githash)
        except git.BadName:
            return None
    else:
        raise ValueError("Must specify either a hash or a tag")
    tags = [i.name for i in g2repo.tags if i.commit == commit]
    return (commit.hexsha, tags, commit.message,commit.committed_datetime)
    
def gitHash2Tags(githash=None,g2repo=None):
    '''Find tags associated with a particular git commit. 
    Note that if `githash` cannot be located because it does not 
    exist or is not unique, a `git.BadName` exception is raised. 

    :param str githash: hex hash code (abbreviated to as few characters as 
       needed to keep it unique). If None (default), the HEAD is used.
    :param str g2repo: git.Rwpo connecton to GSAS-II installation. If 
       None (default) it will be opened. 
    :returns: a list of tags (each a string)
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    if githash is None:
        commit = g2repo.head.object
    else:
        commit = g2repo.commit(githash)
    #return [i.name for i in g2repo.tags if i.commit == commit] # slow with a big repo
    return g2repo.git.tag('--points-at',commit).split('\n')
    
def gitTag2Hash(gittag,g2repo=None):
    '''Provides the hash number for a git tag.
    Note that if `gittag` cannot be located because it does not 
    exist or is too old and is beyond the `depth` of the local 
    repository, a `ValueError` exception is raised. 

    :param str repo_path: location where GSAS-II has been installed.
    :param str gittag: a tag value.
    :param str g2repo: git.Rwpo connecton to GSAS-II installation. If 
       None (default) it will be opened. 
    :returns: a str value with the hex hash for the commit.
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    return g2repo.tag(gittag).commit.hexsha

def gitTestGSASII(verbose=True,g2repo=None):
    '''Test a the status of a GSAS-II installation

    :param bool verbose: if True (default), status messages are printed
    :param str g2repo: git.Rwpo connecton to GSAS-II installation. If 
       None (default) it will be opened. 
    :returns: istat, with the status of the repository, with one of the 
      following values:

       * -1: path is not found
       * -2: no git repository at path
       * -3: unable to access repository

       * value&1==1: repository has local changes (uncommitted/stashed)
       * value&2==2: repository has been regressed (detached head)
       * value&4==4: repository has staged files
       * value&8==8: repository has has been switched to non-master branch

       * value==0:   no problems noted
    '''
    if g2repo is None:
        if not os.path.exists(path2GSAS2): 
            if verbose: print(f'Warning: Directory {path2GSAS2} not found')
            return -1
        if os.path.exists(os.path.join(path2GSAS2,'..','.git')):
            path2repo = os.path.join(path2GSAS2,'..')  # expected location
        elif os.path.exists(os.path.join(path2GSAS2,'.git')):
            path2repo = path2GSAS2
        else:
            if verbose: print(f'Warning: Repository {path2GSAS2} not found')
            return -2
        try:
            g2repo = openGitRepo(path2repo)
        except Exception as msg:
            if verbose: print(f'Warning: Failed to open repository. Error: {msg}')
            return -3
    code = 0
    if g2repo.is_dirty():                     # has changed files
        code += 1
        count_modified_files = len(g2repo.index.diff(None))
    if g2repo.head.is_detached:
        code += 2                             # detached
    else:
        if g2repo.active_branch.name != 'master':
            code += 8                         # not on master branch
    if g2repo.index.diff("HEAD"): code += 4   # staged

    # test if there are local changes committed
    return code

def gitCheckForUpdates(fetch=True,g2repo=None):
    '''Provides a list of the commits made locally and those in the 
    local copy of the repo that have not been applied. Does not 
    provide useful information in the case of a detached Head (see 
    :func:`countDetachedCommits` for that.)

    :param bool fetch: if True (default), updates are copied over from
      the remote repository (git fetch), before checking for changes.
    :param str g2repo: git.Rwpo connecton to GSAS-II installation. If 
       None (default) it will be opened. 
    :returns: a list containing (remotecommits, localcommits, fetched) where 

       * remotecommits is a list of hex hash numbers of remote commits and 
       * localcommits is a list of hex hash numbers of local commits and
       * fetched is a bool that will be True if the update (fetch)
         step ran successfully

       Note that if the head is detached (GSAS-II has been reverted to an 
       older version) or the branch has been changed, the values for each 
       of the three items above will be None.
    '''
    fetched = False
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    if g2repo.head.is_detached:
        return (None,None,None)
    if fetch:
        try:
            g2repo.remote().fetch()
            fetched = True
        except git.GitCommandError as msg:
            print(f'Failed to get updates from {g2repo.remote().url}')
    try:
        head = g2repo.head.ref
        tracking = head.tracking_branch()
        localcommits = [i.hexsha for i in head.commit.iter_items(g2repo, f'{tracking.path}..{head.path}')]
        remotecommits = [i.hexsha for i in head.commit.iter_items(g2repo, f'{head.path}..{tracking.path}')]
        return remotecommits,localcommits,fetched
    except:
        return (None,None,None)
    
def countDetachedCommits(g2repo=None):
    '''Count the number of commits that have been made since
    a commit that is containined in the master branch

    returns the count and the commit object for the 
    parent commit that connects the current stranded
    branch to the master branch.

    None is returned if no connection is found
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    if not g2repo.head.is_detached:
        return 0,g2repo.commit()
    # is detached head in master branch?
    if g2repo.commit() in g2repo.iter_commits('master'):
        return 0,g2repo.commit()
    # count number of commits since leaving master branch
    masterList = list(g2repo.iter_commits('master'))
    for c,i in enumerate(g2repo.commit().iter_parents()):
        if i in masterList:
            return c+1,i
    else:
        return None,None

def gitCountRegressions(g2repo=None):
    '''Count the number of new check ins on the master branch since
    the head was detached as well as any checkins made on the detached
    head. 

    :returns: mastercount,detachedcount, where 

      * mastercount is the number of check ins made on the master branch 
        remote repository since the reverted check in was first made. 
      * detachedcount is the number of check ins made locally 
        starting from the detached head (hopefully 0)

      If the connection between the current head and the master branch 
      cannot be established, None is returned for both.
      If the connection from the reverted check in to the newest version
      (I don't see how this could happen) then only mastercount will be None.
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    # get parent of current head that is in master branch
    detachedcount,parent = countDetachedCommits(g2repo)
    if detachedcount is None: return None,None
    mastercount = 0
    for h in g2repo.iter_commits('master'):
        if h == parent:
            return mastercount,detachedcount
        mastercount += 1
    return None,detachedcount

def gitGetUpdate(mode='Background'):
    '''Download the latest updates into the local copy of the GSAS-II 
    repository from the remote master, but don't actually update the 
    GSAS-II files being used. This can be done immediately or in background. 

    In 'Background' mode, a background process is launched. The results 
    from the process are recorded in file in ~/GSASII_bkgUpdate.log 
    (located in %HOME% on Windows). A pointer to the created process is 
    returned.

    In 'immediate' mode, the update is performed immediately. The 
    function does not return until after the update is downloaded.

    :returns: In 'Background' mode, returns a Popen object (see subprocess). 
      In 'immediate' mode nothing is returned.
    '''
    if mode == 'Background':
        return subprocess.Popen([sys.executable, __file__, '--git-fetch'])
    else:
        g2repo = openGitRepo(path2GSAS2)
        g2repo.remote().fetch()
        if GetConfigValue('debug'): print(f'Updates fetched')

def gitHistory(values='tag',g2repo=None,maxdepth=100):
    '''Provides the history of commits to the master, either as tags 
    or hash values

    :param str values: specifies what type of values are returned. 
      If values=='hash', then hash values or for values=='tag', a 
      list of list of tag(s). 
    :param str g2repo: git.Rwpo connecton to GSAS-II installation. If 
       None (default) it will be opened. 
    :returns: a list of str values where each value is a hash for 
      a commit (values=='hash'), 
      for values=='tag', a list of lists, where a list of tags is provided
      for each commit. When tags are provided, for any commit that does 
      not have any associated tag(s), that entry is omitted from the list.
      for values=='both', a list of lists, where a hash is followed by a 
      list of tags (if any) is provided
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    history = list(g2repo.iter_commits('master'))
    if values.lower().startswith('h'):
        return [i.hexsha for i in history]
    elif values.lower().startswith('t'):
        tagmap = {} # generate lookup table for to get tags
        for t in g2repo.tags:
            tagmap.setdefault(t.commit.hexsha, []).append(t.name)
        return [tagmap[i.hexsha] for i in history if i.hexsha in tagmap]
    elif values.lower().startswith('b'):
        # slow with history >thousands
        # tagmap = {} # generate lookup table for to get tags
        # for t in g2repo.tags:
        #     tagmap.setdefault(t.commit.hexsha, []).append(t.name)
        # return [[i.hexsha]+tagmap.get(i.hexsha,[]) for i in history]
        
        # potentially faster code
        r1 = [[i.hexsha]+g2repo.git.tag('--points-at',i).split('\n')
                    for i in history[:maxdepth]]
        return [[i[0]] if i[1]=='' else i for i in r1]
    else:
        raise ValueError(f'gitHistory has invalid value specified: {value}')

def getGitBinaryReleases(cache=False):
    '''Retrieves the binaries and download urls of the latest release

    :param bool cache: when cache is True and the binaries file names 
       are retrieved (does not always succeed when done via GitHub 
       Actions), the results are saved in a file for reuse should the 
       retrieval fail. Default is False so the file is not changed. 

    :returns: a URL dict for GSAS-II binary distributions found in the newest 
      release in a GitHub repository. The repo location is defined in global 
      `G2binURL`.

      The dict keys are references to binary distributions, which are named 
      as f"{platform}_p{pver}_n{npver}" where platform is determined 
      in :func:`GSASIIpath.GetBinaryPrefix` (linux_64, mac_arm, win_64,...)
      and where `pver` is the Python version (such as "3.10") and `npver` is 
      the numpy version (such as "1.26").

      The value associated with each key contains the full URL to 
      download a tar containing that binary distribution. 
    '''
    # Get first page of releases
    releases = []
    tries = 0
    while tries < 5: # this has been known to fail, so retry
        tries += 1
        releases = requests.get(
            url=f"{G2binURL}/releases", 
            headers=BASE_HEADER
        ).json()
        try:        
            # Get assets of latest release
            assets = requests.get(
                url=f"{G2binURL}/releases/{releases[-1]['id']}/assets",
                headers=BASE_HEADER
                ).json()

            versions = []
            URLs = []
            count = 0
            for asset in assets:
                if asset['name'].endswith('.tgz'):
                    versions.append(asset['name'][:-4]) # Remove .tgz tail
                    URLs.append(asset['browser_download_url'])
                    count += 1
            # Cache the binary releases for later use in case GitHub 
            # prevents us from using a query to get them
            if cache and count > 4:
                fp = open(os.path.join(path2GSAS2,'inputs','BinariesCache.txt'),'w')
                for key in dict(zip(versions,URLs)):
                    fp.write(f'{key} : {res[key]}\n')
                fp.close()
            return dict(zip(versions,URLs))
        except:
            print('Attempt to list GSAS-II binary releases failed, sleeping for 10 sec and then retrying')
            import time
            time.sleep(10)  # this does not seem to help when GitHub is not letting the queries through

    print(f'Could not get releases from {G2binURL}. Using cache')
    res = {}
    try:
        fp = open(os.path.join(path2GSAS2,'inputs','BinariesCache.txt'),'r')
        for line in fp.readlines():
            key,val = line.split(':',1)[:2]
            res[key.strip()] = val.strip()
        fp.close()
        return res
    except:
        raise IOError('Cache read of releases failed too.')
    
def getGitBinaryLoc(npver=None,pyver=None,verbose=True):
    '''Identify the best GSAS-II binary download location from the 
    distributions in the latest release section of the github repository 
    on the CPU platform, and Python & numpy versions. The CPU & Python
    versions must match, but the numpy version may only be close.
    
    :param str npver: Version number to use for numpy, if None (default)
      the version is taken from numpy in the current Python interpreter.
    :param str pyver: Version number to use for Python, if None (default)
      the version is taken from the current Python interpreter.
    :param bool verbose: if True (default), status messages are printed
    :returns: a URL for the tar file (success) or None (failure)
    '''    
    bindir = GetBinaryPrefix(pyver)
    if npver:
        inpver = intver(npver)
    else:
        npver = np.__version__
        inpver = intver(np.__version__)
    # get binaries matching the required install, approximate match for numpy
    URLdict = getGitBinaryReleases()
    versions = {}
    for d in URLdict:
        if d.startswith(bindir):
            v = intver(d.rstrip('/').split('_')[3].lstrip('n'))
            versions[v] = d
    intVersionsList = sorted(versions.keys())
    if not intVersionsList:
        print('No binaries located to match',bindir)
        return
    elif inpver < min(intVersionsList):
        vsel = min(intVersionsList)
        if verbose: print(
                f'Warning: The requested numpy, version, {npver},'
                f' is older than\n\tthe oldest dist version, {fmtver(vsel)}')
    elif inpver >= max(intVersionsList):
        vsel = max(intVersionsList)
        if verbose and inpver == max(intVersionsList):
            print(
                f'The requested numpy version, {npver},'
                f' matches the binary dist, version {fmtver(vsel)}')
        elif verbose:
            print(
                f'Note: using a binary dist for numpy version {fmtver(vsel)} '
                f'which is older than the requested numpy, version {npver}')
    else:
        vsel = min(intVersionsList)
        for v in intVersionsList:
            if v <= inpver:
                vsel = v
            else:
                if verbose: print(
                        f'FYI: Selecting dist version {fmtver(vsel)}'
                        f' as the requested numpy, version, {npver},'
                        f'\n\tis older than the next dist version {fmtver(v)}')
                break
    return URLdict[versions[vsel]]

def InstallGitBinary(tarURL, instDir, nameByVersion=False, verbose=True):
    '''Install the GSAS-II binary files into the location
    specified. 
    
    :param str tarURL: a URL for the tar file.
    :param str instDir: location directory to install files. This directory
        may not exist and will be created if needed.
    :param bool nameByVersion: if True, files are put into a subdirectory
        of `instDir`, named to match the tar file (with plaform, Python & 
        numpy versions). 
        Default is False, where the binary files are put directly into 
        `instDir`.
    :param bool verbose: if True (default), status messages are printed.
    :returns: None
    '''
    # packages not commonly used so import them here not on startup
    import requests
    import tempfile
    import tarfile
    # download to scratch
    tar = tempfile.NamedTemporaryFile(suffix='.tgz',delete=False)
    try:
        tar.close()
        if verbose: print(f'Downloading {tarURL}')
        r = requests.get(tarURL, allow_redirects=True)
        open(tar.name, 'wb').write(r.content)
        # open in tar
        tarobj = tarfile.open(name=tar.name)
        if nameByVersion:
            binnam = os.path.splitext(os.path.split(tarURL)[1])[0]
            install2dir = os.path.join(instDir,binnam)
        else:
            install2dir = instDir
        for f in tarobj.getmembers(): # loop over files in repository
            # do a bit of sanity checking for safety. Don't install anything
            #  unless it goes into in the specified directory
            if sys.platform == "win32" and f.name.startswith('._'): continue # clean up Mac cruft for Windows
            if '/' in f.name or '\\' in f.name:
                print(f'skipping file {f.name} -- path alteration not allowed')
                continue
            if f.name != os.path.basename(f.name):
                print(f'skipping file {f.name} -- how did this happen?')
                continue
            newfil = os.path.normpath(os.path.join(install2dir,f.name))
            tarobj.extract(f, path=install2dir, set_attrs=False)
            # set file mode and mod/access times (but not ownership)
            os.chmod(newfil,f.mode)
            os.utime(newfil,(f.mtime,f.mtime))
            if verbose: print(f'Created GSAS-II binary file {newfil}')
    finally:
        del tarobj
        os.unlink(tar.name)

def GetRepoUpdatesInBackground():
    '''Wrapper to make sure that :func:`gitGetUpdate` is called only 
    if git has been used to install GSAS-II.
    
    :returns: returns a Popen object (see subprocess)
    '''
    if HowIsG2Installed().startswith('git'):
        return gitGetUpdate(mode='Background')

def gitStartUpdate(cmdopts):
    '''Update GSAS-II in a separate process, by running this script with the 
    options supplied in the call to this function and then exiting GSAS-II. 
    '''
    cmd = [sys.executable, __file__] + cmdopts
    if GetConfigValue('debug'): print('Starting updates with command\n\t'+
                                      f'{" ".join(cmd)}')
    proc = subprocess.Popen(cmd)
    # on windows the current process needs to end so that the source files can
    # be written over. On unix the current process needs to stay running
    # so the child is not killed.
    if sys.platform != "win32": proc.wait()
    sys.exit()

def dirGitHub(dirlist,orgName=gitTutorialOwn, repoName=gitTutorialRepo):
    '''Obtain a the contents of a GitHub repository directory using 
    the GitHub REST API.

    :param str dirlist: a list of sub-directories `['parent','child',sub']` 
      for `parent/child/sub` or `[]` for a file in the top-level directory.     
    :param str orgName: the name of the GitHub organization
    :param str repoName: the name of the GitHub repository
    :returns: a list of file names or None if the dirlist info does not 
      reference a directory

    examples::

        dirGitHub([], 'GSASII', 'TutorialTest')
        dirGitHub(['TOF Sequential Single Peak Fit', 'data'])

    The first example will get the contents of the top-level 
    directory for the specified repository 

    The second example will provide the contents of the 
    "TOF Sequential Single Peak Fit"/data directory. 
    '''
    dirname = ''
    for item in dirlist:
        dirname += item + '/'
    URL = f"https://api.github.com/repos/{orgName}/{repoName}/contents/{dirname}"
    r = requests.get(URL, allow_redirects=True)
    try:
        return [rec['name'] for rec in r.json()]
    except:
        return None

def rawGitHubURL(dirlist,filename,orgName=gitTutorialOwn, repoName=gitTutorialRepo,
                 branchname="master"):
    '''Create a URL that can be used to view/downlaod the raw version of 
    file in a GitHub repository. 

    :param str dirlist: a list of sub-directories `['parent','child',sub']` 
      for `parent/child/sub` or `[]` for a file in the top-level directory.     
    :param str filename: the name of the file
    :param str orgName: the name of the GitHub organization
    :param str repoName: the name of the GitHub repository
    :param str branchname: the name of the GitHub branch. Defaults 
       to "master".

    :returns: a URL-encoded URL
    '''
    import urllib.parse  # not used very often, import only when needed
    dirname = ''
    for item in dirlist:
        # it's not clear that the URLencode is needed for the directory name
        dirname += urllib.parse.quote(item) + '/'
        #filename = urllib.parse.quote(filename)
    return f"https://raw.githubusercontent.com/{orgName}/{repoName}/{branchname}/{dirname}{filename}"

def downloadDirContents(dirlist,targetDir,orgName=gitTutorialOwn, repoName=gitTutorialRepo):
    '''Download the entire contents of a directory from a repository 
    on GitHub. Used to download data for a tutorial.
    '''
    filList = dirGitHub(dirlist, orgName=orgName, repoName=repoName)
    if filList is None:
        print(f'Directory {"/".join(dirlist)!r} does not have any files')
        return None
    for fil in filList:
        if fil.lower() == 'index.html': continue
        URL = rawGitHubURL(dirlist,fil,orgName=orgName,repoName=repoName)
        r = requests.get(URL, allow_redirects=True)
        outfil = os.path.join(targetDir,fil)
        if r.status_code == 200:
            open(outfil, 'wb').write(r.content)
            print(f'wrote {outfil}')
        elif r.status_code == 404:
            print(f'Warning: {fil} is likely a subdirectory of directory {"/".join(dirlist)!r}')
        else:
            print(f'Unexpected web response for {fil}: {r.status_code}')
    return

#==============================================================================
#==============================================================================
# routines to interface with subversion
g2home = 'https://subversion.xray.aps.anl.gov/pyGSAS' # 'Define the location of the GSAS-II subversion repository'
proxycmds = [] # 'Used to hold proxy information for subversion, set if needed in whichsvn'
svnLocCache = None  # 'Cached location of svn to avoid multiple searches for it'

def MakeByte2str(arg):
    '''Convert output from subprocess pipes (bytes) to str (unicode) in Python 3.
    In Python 2: Leaves output alone (already str). 
    Leaves stuff of other types alone (including unicode in Py2)
    Works recursively for string-like stuff in nested loops and tuples.

    typical use::

        out = MakeByte2str(out)

    or::

        out,err = MakeByte2str(s.communicate())
    
    '''
    if isinstance(arg,str): return arg
    if isinstance(arg,bytes):
        try:
            return arg.decode()
        except:
            if GetConfigValue('debug'): print('Decode error')
            return arg
    if isinstance(arg,list):
        return [MakeByte2str(i) for i in arg]
    if isinstance(arg,tuple):
        return tuple([MakeByte2str(i) for i in arg])
    return arg
                
def getsvnProxy():
    '''Loads a proxy for subversion from the proxyinfo.txt file created 
    by bootstrap.py or File => Edit Proxy...; If not found, then the 
    standard http_proxy and https_proxy environment variables are scanned
    (see https://docs.python.org/3/library/urllib.request.html#urllib.request.getproxies) 
    with case ignored and that is used. 
    '''
    global proxycmds
    proxycmds = []
    proxyinfo = os.path.join(os.path.expanduser('~/.G2local/'),"proxyinfo.txt")
    if not os.path.exists(proxyinfo):
        proxyinfo = os.path.join(path2GSAS2,"proxyinfo.txt")
    if os.path.exists(proxyinfo):
        fp = open(proxyinfo,'r')
        host = fp.readline().strip()
        # allow file to begin with comments
        while host.startswith('#'):
            host = fp.readline().strip()
        port = fp.readline().strip()
        etc = []
        line = fp.readline()
        while line:
            etc.append(line.strip())
            line = fp.readline()
        fp.close()
        setsvnProxy(host,port,etc)
        return host,port,etc
    import urllib.request
    proxdict = urllib.request.getproxies()
    varlist = ("https","http")
    for var in proxdict:
        if var.lower() in varlist:
            proxy = proxdict[var]
            pl = proxy.split(':')
            if len(pl) < 2: continue
            host = pl[1].strip('/')
            port = ''
            if len(pl) == 3:
                port = pl[2].strip('/').strip()
            return host,port,''
    return '','',''

def setsvnProxy(host,port,etc=[]):
    '''Sets the svn commands needed to use a proxy
    '''
    global proxycmds
    proxycmds = []
    host = host.strip()
    port = port.strip()
    if host: 
        proxycmds.append('--config-option')
        proxycmds.append('servers:global:http-proxy-host='+host)
        if port:
            proxycmds.append('--config-option')
            proxycmds.append('servers:global:http-proxy-port='+port)
    for item in etc:
        proxycmds.append(item)
        
def whichsvn():
    '''Returns a path to the subversion exe file, if any is found.
    Searches the current path after adding likely places where GSAS-II
    might install svn. 

    :returns: None if svn is not found or an absolute path to the subversion
      executable file.
    '''
    # use a previosuly cached svn location
    global svnLocCache
    if svnLocCache: return svnLocCache
    # prepare to find svn
    is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    svnprog = 'svn'
    if sys.platform.startswith('win'): svnprog += '.exe'
    host,port,etc = getsvnProxy()
    if GetConfigValue('debug') and host:
        print('DBG_Using proxy host {} port {}'.format(host,port))
    if GetConfigValue('svn_exec'):
        exe_file = GetConfigValue('svn_exec')
        print('Using ',exe_file)
        if is_exe(exe_file):
            try:
                p = subprocess.Popen([exe_file,'help'],stdout=subprocess.PIPE)
                res = p.stdout.read()
                if not res: return
                p.communicate()
                svnLocCache = os.path.abspath(exe_file)
                return svnLocCache
            except:
                pass
    # add likely places to find subversion when installed with GSAS-II
    pathlist = os.environ["PATH"].split(os.pathsep)
    pathlist.insert(0,os.path.split(sys.executable)[0])
    pathlist.insert(1,path2GSAS2)
    for rpt in ('..','bin'),('..','Library','bin'),('svn','bin'),('svn',),('.'):
        pt = os.path.normpath(os.path.join(path2GSAS2,*rpt))
        if os.path.exists(pt):
            pathlist.insert(0,pt)    
    # search path for svn or svn.exe
    for path in pathlist:
        exe_file = os.path.join(path, svnprog)
        if is_exe(exe_file):
            try:
                p = subprocess.Popen([exe_file,'help'],stdout=subprocess.PIPE)
                res = p.stdout.read()
                if not res: return
                p.communicate()
                svnLocCache = os.path.abspath(exe_file)
                return svnLocCache
            except:
                pass        
    svnLocCache = None

svn_version = None
def svnVersion(svn=None):
    '''Get the version number of the current subversion executable.
    The result is cached, as this takes a bit of time to run and 
    is done a fair number of times.

    :returns: a string with a version number such as "1.6.6" or None if
      subversion is not found.

    '''
    global svn_version
    if svn_version is not None:
        return svn_version
    if not svn: svn = whichsvn()
    if not svn: return

    cmd = [svn,'--version','--quiet']
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print ('subversion error!\nout=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        return None
    svn_version = out.strip()
    return svn_version

def svnVersionNumber(svn=None):
    '''Get the version number of the current subversion executable

    :returns: a fractional version number such as 1.6 or None if
      subversion is not found.

    '''
    ver = svnVersion(svn)
    if not ver: return 
    M,m = ver.split('.')[:2]
    return int(M)+int(m)/10.

def svnGetLog(fpath=os.path.split(__file__)[0],version=None):
    '''Get the revision log information for a specific version of the specified package

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located.
    :param int version: the version number to be looked up or None (default)
       for the latest version.

    :returns: a dictionary with keys (one hopes) 'author', 'date', 'msg', and 'revision'

    '''
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    if version is not None:
        vstr = '-r'+str(version)
    else:
        vstr = '-rHEAD'

    cmd = [svn,'log',fpath,'--xml',vstr]
    if proxycmds: cmd += proxycmds
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print ('out=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        return None
    x = ET.fromstring(out)
    d = {}
    for i in x.iter('logentry'):
        d = {'revision':i.attrib.get('revision','?')}
        for j in i:
            d[j.tag] = j.text
        break # only need the first
    return d

svnLastError = ''
def svnGetRev(fpath=os.path.split(__file__)[0],local=True,verbose=True):
    '''Obtain the version number for the either the last update of the local version
    or contacts the subversion server to get the latest update version (# of Head).

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    :param bool local: determines the type of version number, where
       True (default): returns the latest installed update 
       False: returns the version number of Head on the server

    :Returns: the version number as an str or 
       None if there is a subversion error (likely because the path is
       not a repository or svn is not found). The error message is placed in
       global variable svnLastError
    '''

    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    if local:
        cmd = [svn,'info',fpath,'--xml']
    else:
        cmd = [svn,'info',fpath,'--xml','-rHEAD']
    if svnVersionNumber() >= 1.6:
        cmd += ['--non-interactive', '--trust-server-cert']
    if proxycmds: cmd += proxycmds
    # if GetConfigValue('debug'):
    #     s = 'subversion command:\n  '
    #     for i in cmd: s += i + ' '
    #     print(s)
    s = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        if verbose:
            print ('svn failed\n%s'%out)
            print ('err=%s'%err)
            print('\nsvn command:',' '.join(cmd))
        global svnLastError
        svnLastError = err
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
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if not svn: return
    cmd = [svn,'status',fpath,'--xml']
    if proxycmds: cmd += proxycmds
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err: return None
    x = ET.fromstring(out)
    changed = []
    for i in x.iter('entry'):
        if i.find('wc-status').attrib.get('item') == 'modified': 
            changed.append(i.attrib.get('path'))
    return changed

def svnCleanup(fpath=os.path.split(__file__)[0],verbose=True):
    '''This runs svn cleanup on a selected local directory. 

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    '''
    svn = whichsvn()
    if not svn: return
    if verbose: print(u"Performing svn cleanup at "+fpath)
    cmd = [svn,'cleanup',fpath]
    if verbose: showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print(60*"=")
        print("****** An error was noted, see below *********")
        print(60*"=")
        print(err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        #raise Exception('svn cleanup failed')
        return False
    elif verbose:
        print(out)
    return True
        
def svnUpdateDir(fpath=os.path.split(__file__)[0],version=None,verbose=True):
    '''This performs an update of the files in a local directory from a server. 

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    :param version: the number of the version to be loaded. Used only
       cast as a string, but should be an integer or something that corresponds to a
       string representation of an integer value when cast. A value of None (default)
       causes the latest version on the server to be used.
    '''
    svn = whichsvn()
    if not svn: return
    if version:
        verstr = '-r' + str(version)
    else:
        verstr = '-rHEAD'
    if verbose: print(u"Updating files at "+fpath)
    cmd = [svn,'update',fpath,verstr,
           '--non-interactive',
           '--accept','theirs-conflict','--force']
    if svnVersionNumber() >= 1.6:
        cmd += ['--trust-server-cert']
    if proxycmds: cmd += proxycmds
    #if verbose or GetConfigValue('debug'):
    if verbose:
        s = 'subversion command:\n  '
        for i in cmd: s += i + ' '
        print(s)
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print(60*"=")
        print("****** An error was noted, see below *********")
        print(60*"=")
        print(err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        if svnCleanup(fpath):
            s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = MakeByte2str(s.communicate())
            if err:
                print(60*"=")
                print("****** Drat, failed again: *********")
                print(60*"=")
                print(err)
            else:
                return
        if 'Checksum' in err:  # deal with Checksum problem
            err = svnChecksumPatch(svn,fpath,verstr)
            if err:
                print('error from svnChecksumPatch\n\t',err)
            else:
                return
        raise Exception('svn update failed')
    elif verbose:
        print(out)

def showsvncmd(cmd):
    s = '\nsvn command:  '
    for i in cmd: s += i + ' '
    print(s)

def svnChecksumPatch(svn,fpath,verstr):
    '''This performs a fix when svn cannot finish an update because of
    a Checksum mismatch error. This seems to be happening on OS X for 
    unclear reasons. 
    '''
    print('\nAttempting patch for svn Checksum mismatch error\n')
    svnCleanup(fpath)
    cmd = [svn,'update','--set-depth','empty',
               os.path.join(fpath,'bindist')]
    showsvncmd(cmd)        
    if svnVersionNumber() >= 1.6:
        cmd += ['--non-interactive', '--trust-server-cert']
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    cmd = [svn,'switch',g2home+'/trunk/bindist',
               os.path.join(fpath,'bindist'),
               '--non-interactive', '--trust-server-cert', '--accept',
               'theirs-conflict', '--force', '-rHEAD', '--ignore-ancestry']
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    DownloadG2Binaries(g2home,verbose=True)
    cmd = [svn,'update','--set-depth','infinity',
               os.path.join(fpath,'bindist')]
    if svnVersionNumber() >= 1.6:
        cmd += ['--non-interactive', '--trust-server-cert']
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    cmd = [svn,'update',fpath,verstr,
                       '--non-interactive',
                       '--accept','theirs-conflict','--force']
    if svnVersionNumber() >= 1.6:
        cmd += ['--trust-server-cert']
    if proxycmds: cmd += proxycmds
    showsvncmd(cmd)        
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    #if err: print('error=',err)
    return err
        
def svnUpgrade(fpath=os.path.split(__file__)[0]):
    '''This reformats subversion files, which may be needed if an upgrade of subversion is
    done. 

    :param str fpath: path to repository dictionary, defaults to directory where
       the current file is located
    '''
    svn = whichsvn()
    if not svn: return
    cmd = [svn,'upgrade',fpath,'--non-interactive']
    if proxycmds: cmd += proxycmds
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print("svn upgrade did not happen (this is probably OK). Messages:")
        print (err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)

def svnUpdateProcess(version=None,projectfile=None,branch=None):
    '''perform an update of GSAS-II in a separate python process'''
    if not projectfile:
        projectfile = ''
    else:
        projectfile = os.path.realpath(projectfile)
        print ('restart using %s'%projectfile)
    if branch:
        version = branch
    elif not version:
        version = ''
    else:
        version = str(version)
    # start the upgrade in a separate interpreter (avoids loading .pyd files)
    ex = sys.executable
    if sys.platform == "darwin": # mac requires pythonw which is not always reported as sys.executable
        if os.path.exists(ex+'w'): ex += 'w'
    proc = subprocess.Popen([ex,__file__,projectfile,version])
    if sys.platform != "win32":
        proc.wait()
    sys.exit()

def svnSwitchDir(rpath,filename,baseURL,loadpath=None,verbose=True):
    '''This performs a switch command to move files between subversion trees.
    Note that if the files were previously downloaded, 
    the switch command will update the files to the newest version.
    
    :param str rpath: path to locate files, relative to the GSAS-II
      installation path (defaults to path2GSAS2)
    :param str URL: the repository URL
    :param str loadpath: the prefix for the path, if specified. Defaults to path2GSAS2
    :param bool verbose: if True (default) diagnostics are printed
    '''
    svn = whichsvn()
    if not svn: return
    URL = baseURL[:]
    if baseURL[-1] != '/':
        URL = baseURL + '/' + filename
    else:
        URL = baseURL + filename
    if loadpath:
        fpath = os.path.join(loadpath,rpath,filename)
        svntmp = os.path.join(loadpath,'.svn','tmp')
    else:
        fpath = os.path.join(path2GSAS2,rpath,filename)
        svntmp = os.path.join(path2GSAS2,'.svn','tmp')
    # fix up problems with missing empty directories
    if not os.path.exists(fpath):
        print('Repairing missing directory',fpath)
        cmd = [svn,'revert',fpath]
        s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
        out,err = MakeByte2str(s.communicate())
        if out: print(out)
        if err: print(err)
    if not os.path.exists(svntmp):
        print('Repairing missing directory',svntmp)
        cmd = ['mkdir',svntmp]
        s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
        out,err = MakeByte2str(s.communicate())
        if out: print(out)
        if err: print(err)
        
    cmd = [svn,'switch',URL,fpath,
           '--non-interactive','--trust-server-cert',
           '--accept','theirs-conflict','--force','-rHEAD']
    if svnVersionNumber(svn) > 1.6: cmd += ['--ignore-ancestry']
    if proxycmds: cmd += proxycmds
    if verbose:
        print(u"Loading files to "+fpath+u"\n  from "+URL)
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print(60*"=")
        print ("****** An error was noted, see below *********")
        print(60*"=")
        print ('out=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        if svnCleanup(fpath):
            s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = MakeByte2str(s.communicate())
            if err:
                print(60*"=")
                print("****** Drat, failed again: *********")
                print(60*"=")
                print(err)
            else:
                return True
        return False
    if verbose:
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        print('\n=== Output from svn switch'+(43*'='))
        print(out.strip())
        print((70*'=')+'\n')
    return True

def svnInstallDir(URL,loadpath):
    '''Load a subversion tree into a specified directory

    :param str URL: the repository URL
    :param str loadpath: path to locate files

    '''
    svn = whichsvn()
    if not svn: return
    cmd = [svn,'co',URL,loadpath,'--non-interactive']
    if svnVersionNumber() >= 1.6: cmd += ['--trust-server-cert']
    print("Loading files from "+URL)
    if proxycmds: cmd += proxycmds
    s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())   #this fails too easily
    if err:
        print(60*"=")
        print ("****** An error was noted, see below *********")
        print(60*"=")
        print (err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        if svnCleanup(loadpath):
            s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = MakeByte2str(s.communicate())
            if err:
                print(60*"=")
                print("****** Drat, failed again: *********")
                print(60*"=")
                print(err)
                return False
        else:
            return False
    print ("Files installed at: "+loadpath)
    return True

def svnGetFileStatus(fpath=os.path.split(__file__)[0],version=None):
    '''Compare file status to repository (svn status -u)

    :returns: updatecount,modcount,locked where 
       updatecount is the number of files waiting to be updated from 
       repository 
       modcount is the number of files that have been modified locally
       locked  is the number of files tagged as locked
    '''
    import xml.etree.ElementTree as ET
    svn = whichsvn()
    if version is not None:
        vstr = '-r'+str(version)
    else:
        vstr = '-rHEAD'
    cmd = [svn,'st',fpath,'--xml','-u',vstr]
    if svnVersionNumber() >= 1.6:
        cmd += ['--non-interactive', '--trust-server-cert']
    if proxycmds: cmd += proxycmds
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print ('out=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        return None

    locked = 0
    updatecount = 0
    modcount = 0
    x = ET.fromstring(out)
    for i0 in x.iter('entry'):
        filename = i0.attrib.get('path','?')
        wc_rev = ''
        status = ''
        switched = ''
        for i1 in i0.iter('wc-status'):
            wc_rev = i1.attrib.get('revision','')
            status = i1.attrib.get('item','')
            switched = i1.attrib.get('switched','')
            if i1.attrib.get('wc-locked',''): locked += 1
        if status == "unversioned": continue
        if switched == "true": continue
        if status == "modified":
            modcount += 1
        elif status == "normal":
            updatecount += 1
        file_rev = ''
        for i2 in i1.iter('commit'):
            file_rev = i2.attrib.get('revision','')
        local_status = ''
        for i1 in i0.iter('repos-status'):
            local_status = i1.attrib.get('item','')
        #print(filename,wc_rev,file_rev,status,local_status,switched)
    return updatecount,modcount,locked

def svnList(URL,verbose=True):
    '''Get a list of subdirectories from and svn repository
    '''    
    svn = whichsvn()
    if not svn:
        print('**** unable to load files: svn not found ****')
        return ''
    # get binaries matching the required type -- other than for the numpy version
    cmd = [svn, 'list', URL,'--non-interactive', '--trust-server-cert']
    if proxycmds: cmd += proxycmds
    if verbose:
        s = 'Running svn command:\n  '
        for i in cmd: s += i + ' '
        print(s)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res,err = MakeByte2str(p.communicate())
    return res

def DownloadG2Binaries(g2home,verbose=True):
    '''Download GSAS-II binaries from appropriate section of the
    GSAS-II svn repository based on the platform, numpy and Python
    version
    '''    
    bindir = GetBinaryPrefix()
    #npver = 'n{}.{}'.format(*np.__version__.split('.')[0:2])
    inpver = intver(np.__version__)
    svn = whichsvn()
    if not svn:
        print('**** unable to load files: svn not found ****')
        return ''
    # get binaries matching the required type -- other than for the numpy version
    cmd = [svn, 'list', g2home + '/Binaries/','--non-interactive', '--trust-server-cert']
    if proxycmds: cmd += proxycmds
    if verbose:
        s = 'Running svn command:\n  '
        for i in cmd: s += i + ' '
        print(s)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res,err = MakeByte2str(p.communicate())
    versions = {}
    for d in res.split():
        if d.startswith(bindir):
            v = intver(d.rstrip('/').split('_')[3].lstrip('n'))
            versions[v] = d
    intVersionsList = sorted(versions.keys())
    if not intVersionsList:
        print('No binaries located matching',bindir)
        return
    elif inpver < min(intVersionsList):
        vsel = min(intVersionsList)
        print('Warning: The current numpy version, {}, is older than\n\tthe oldest dist version, {}'
              .format(np.__version__,fmtver(vsel)))
    elif inpver >= max(intVersionsList):
        vsel = max(intVersionsList)
        if verbose: print(
                'FYI: The current numpy version, {}, is newer than the newest dist version {}'
                .format(np.__version__,fmtver(vsel)))
    else:
        vsel = min(intVersionsList)
        for v in intVersionsList:
            if v <= inpver:
                vsel = v
            else:
                if verbose: print(
                        'FYI: Selecting dist version {} as the current numpy version, {},\n\tis older than the next dist version {}'
                        .format(fmtver(vsel),np.__version__,fmtver(v)))
                break
    distdir = g2home + '/Binaries/' + versions[vsel]
    # switch reset command: distdir = g2home + '/trunk/bindist'
    svnSwitchDir('bindist','',distdir,verbose=verbose)
    return os.path.join(path2GSAS2,'bindist')

# def svnTestBranch(loc=None):
#     '''Returns the name of the branch directory if the installation has been switched.
#     Returns none, if not a branch
#     the test 2frame branch. False otherwise
#     '''
#     if loc is None: loc = path2GSAS2
#     svn = whichsvn()
#     if not svn:
#         print('**** unable to load files: svn not found ****')
#         return ''
#     cmd = [svn, 'info', loc]
#     if proxycmds: cmd += proxycmds
#     p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#     res,err = MakeByte2str(p.communicate())
#     for l in res.split('\n'):
#         if "Relative URL:" in l: break
#     if "/branch/" in l:
#         return l[l.find("/branch/")+8:].strip()
#     else:
#         return None
    
def svnSwitch2branch(branch=None,loc=None,svnHome=None):
    '''Switch to a subversion branch if specified. Switches to trunk otherwise.
    '''
    if svnHome is None: svnHome = g2home
    svnURL = svnHome + '/trunk'
    if branch:
        if svnHome.endswith('/'):
            svnURL = svnHome[:-1]
        else:
            svnURL = svnHome
        if branch.startswith('/'):
            svnURL += branch
        else:
            svnURL += '/' + branch
    svnSwitchDir('','',svnURL,loadpath=loc)
#==============================================================================
#==============================================================================
    
def runScript(cmds=[], wait=False, G2frame=None):
    '''run a shell script of commands in an external process
    
    :param list cmds: a list of str's, each ietm containing a shell (cmd.exe
      or bash) command
    :param bool wait: if True indicates the commands should be run and then 
      the script should return. If False, then the currently running Python 
      will exit. Default is False
    :param wx.Frame G2frame: provides the location of the current .gpx file
      to be used to restart GSAS-II after running the commands, if wait 
      is False. Default is None which prevents restarting GSAS-II regardless of
      the value of wait.
    '''
    import tempfile
    if not cmds:  #debug
        print('nothing to do in runScript')
        return
    if sys.platform != "win32":
        suffix = '.sh'
    else:
        suffix = '.bat'
        
    fp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    shellname = fp.name
    for line in cmds:
        fp.write(line)
        fp.write('\n')

    if not wait:
        if G2frame:
            projectfile = ''
            if G2frame.GSASprojectfile:
                projectfile = os.path.realpath(G2frame.GSASprojectfile)
            main = os.path.join(path2GSAS2,'GSASII.py')
            ex = sys.executable
            if sys.platform == "darwin": # mac requires pythonw which is not always reported as sys.executable
                if os.path.exists(ex+'w'): ex += 'w'
            print ('restart using ',' '.join([ex,main,projectfile]))
            fp.write(' '.join([ex,main,projectfile]))
            fp.write('\n')
    fp.close()

    # start the upgrade in a separate interpreter (avoids loading .pyd files)
    if sys.platform != "win32":
        proc = subprocess.Popen(['bash',shellname])
    else:
        proc = subprocess.Popen([shellname],shell=True)
    if wait:
        proc.wait()
    else:
        if sys.platform != "win32": proc.wait()
        sys.exit()

def IPyBreak_base(userMsg=None):
    '''A routine that invokes an IPython session at the calling location
    This routine is only used when debug=True is set in config.py
    '''
    savehook = sys.excepthook # save the exception hook
    try: 
        from IPython.terminal.embed import InteractiveShellEmbed
    except ImportError:
        try:
            # try the IPython 0.12 approach
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
        except ImportError:
            print ('IPython InteractiveShellEmbed not found')
            return
    import inspect
    #from IPython import __version__
    #if __version__.startswith('8.12.'): # see https://github.com/ipython/ipython/issues/13966
    from IPython.core import getipython
    if getipython.get_ipython() is None:
        ipshell = InteractiveShellEmbed.instance()
    else:
        ipshell = InteractiveShellEmbed()

    frame = inspect.currentframe().f_back
    msg   = 'Entering IPython console inside {0.f_code.co_filename} at line {0.f_lineno}\n'.format(frame)
    if userMsg: msg += userMsg
    # globals().update(locals()) # This might help with vars inside list comprehensions, etc.
    ipshell(msg,stack_depth=2) # Go up one level, to see the calling routine
    sys.excepthook = savehook # reset IPython's change to the exception hook

def exceptHook(*args):
    '''A routine to be called when an exception occurs. It prints the traceback
    with fancy formatting and then calls an IPython shell with the environment
    of the exception location.
    
    This routine is only used when debug=True is set in config.py    
    '''
    try:
        from IPython.core import ultratb
    except:
        pass

    try: 
        from IPython.terminal.embed import InteractiveShellEmbed
        import IPython.core
        if sys.platform.startswith('win'):
            IPython.core.ultratb.FormattedTB(call_pdb=False,color_scheme='NoColor')(*args)
        else:
            IPython.core.ultratb.FormattedTB(call_pdb=False,color_scheme='LightBG')(*args)
        from IPython.core import getipython
        if getipython.get_ipython() is None:
            ipshell = InteractiveShellEmbed.instance()
        else:
            ipshell = InteractiveShellEmbed()
    except ImportError:
        print ('IPython not installed or is really old')
        return

    import inspect
    frame = inspect.getinnerframes(args[2])[-1][0]
    msg   = 'Entering IPython console at {0.f_code.co_filename} at line {0.f_lineno}\n'.format(frame)
    savehook = sys.excepthook # save the exception hook
    try:
        ipshell(msg,local_ns=frame.f_locals,global_ns=frame.f_globals) # newest (IPython >= 8)
    except DeprecationWarning: # IPython <=7
        try: # IPython >=5
            class c(object): pass
            pseudomod = c() # create something that acts like a module
            pseudomod.__dict__ = frame.f_locals
            InteractiveShellEmbed(banner1=msg)(module=pseudomod,global_ns=frame.f_globals)
        except: # 'IPython <5
            InteractiveShellEmbed(banner1=msg)(local_ns=frame.f_locals,global_ns=frame.f_globals)
    sys.excepthook = savehook # reset IPython's change to the exception hook

def DoNothing():
    '''A routine that does nothing. This is called in place of IPyBreak and pdbBreak
    except when the debug option is set True in config.py
    '''
    pass 

def InvokeDebugOpts():
    'Called in GSASII.py to set up debug options'
    if any('SPYDER' in name for name in os.environ):
        print('Running from Spyder, keeping breakpoint() active & skipping exception trapping')
    elif GetConfigValue('debug'):
        try:
            import pdb
            global pdbBreak
            pdbBreak = pdb.set_trace
            import IPython
            global IPyBreak
            IPyBreak = IPyBreak_base
            sys.excepthook = exceptHook
            os.environ['PYTHONBREAKPOINT'] = 'GSASIIpath.IPyBreak_base'
            print ('Debug on: IPython: Exceptions and G2path.IPyBreak(); pdb: G2path.pdbBreak()')
        except:
            print ('Debug on failed. IPython not installed?')
    else: # not in spyder or debug enabled, hide breakpoints
        os.environ['PYTHONBREAKPOINT'] = '0'

def TestSPG(fpth):
    '''Test if pyspg.[so,.pyd] can be run from a location in the path
    '''
    try:
        if not os.path.exists(fpth): return False
        if not glob.glob(os.path.join(fpth,'pyspg.*')): return False
    except:
        return False
    savpath = sys.path[:]
    sys.path = [fpth]
    # test to see if a shared library can be used
    try:
        import pyspg
    except ModuleNotFoundError as err:
        print(70*'=')
        print(f'Binary module pyspg not found in {fpth!r}\nerror msg: {err}')
        print(70*'=')
        sys.path = savpath
        return False
    except ImportError as err:
        print(70*'=')
        print(f'Module pyspg in {fpth!r} could not be loaded\nerror msg: {err}')
        print(70*'=')
        sys.path = savpath
        return False
    except Exception as err:
        print(70*'=')
        print(f'Error importing module pyspg in {fpth!r}\nerror msg: {err}')
        print(70*'=')
        sys.path = savpath
        return False
    try:
        pyspg.sgforpy('P -1')
    except Exception as err:
        print(70*'=')
        print(f'Module pyspg in {fpth!r} could not be run\nerror msg: {err}')
        print(70*'=')
        sys.path = savpath
        return False
    sys.path = savpath
    return True
    
def SetBinaryPath(printInfo=False, loadBinary=False):
    '''
    Add location of GSAS-II shared libraries (binaries: .so or 
    .pyd files) to path
    
    This routine must be executed after GSASIIpath is imported 
    and before any other GSAS-II imports are done, since 
    they assume binary files are in path

    :param bool printInfo: When True, information is printed to show
      has happened (default is False)
    :param bool loadBinary: no longer in use. This is now done in 
      :func:`GSASIIdataGUI.ShowVersions`.
    '''
    # run this routine only once no matter how many times it is called
    # using the following globals to check for repeated calls
    global BinaryPathLoaded,binaryPath,BinaryPathFailed
    if BinaryPathLoaded or BinaryPathFailed: return
    try:
        inpver = intver(np.__version__)
    except (AttributeError,TypeError): # happens on building docs
        return
    if path2GSAS2 not in sys.path:
        sys.path.insert(0,path2GSAS2)  # make sure current path is used
    binpath = None
    binprfx = GetBinaryPrefix()
    # places to look for the GSAS-II binary directory
    binseapath = [os.path.abspath(sys.path[0])]  # where Python is installed
    binseapath += [os.path.abspath(os.path.dirname(__file__))]  # directory where this file is found
    binseapath += [os.path.dirname(binseapath[-1])]  # parent of above directory
    binseapath += [os.path.expanduser('~/.GSASII')]       # directory in user's home
    def appendIfExists(searchpathlist,loc,subdir):
        newpath = os.path.join(loc,subdir)
        if os.path.exists(newpath):
            if newpath in searchpathlist: return
            searchpathlist.append(newpath)
    searched = []
    for loc in binseapath:
        if loc in searched: continue
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
            v = intver(d.rstrip('/').split('_')[-1].lstrip('n'))
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
            if TestSPG(fpth):
                binpath = fpth   # got one that works, look no farther!
                break
        else:
            continue
        break
    if binpath:  # were GSAS-II binaries found?
        binaryPath = binpath
        BinaryPathLoaded = True
    else:
        print('*** ERROR: Unable to find GSAS-II binaries. Much of GSAS-II cannot function')
        BinaryPathFailed = True
        return None

    # add the data import and export directory to the search path
    if binpath not in sys.path: sys.path.insert(0,binpath)
    if printInfo: print(f'GSAS-II binary directory: {binpath}')
    newpath = os.path.join(path2GSAS2,'imports')
    if newpath not in sys.path: sys.path.append(newpath)
    newpath = os.path.join(path2GSAS2,'exports')
    if newpath not in sys.path: sys.path.append(newpath)
    LoadConfig(printInfo)

def LoadConfig(printInfo=True):
    # setup read of config.py, if present
    global configDict
    try:
        import config
        configDict = config.__dict__
        import inspect
        vals = [True for i in inspect.getmembers(config) if '__' not in i[0]]
        if printInfo:
            print (str(len(vals))+' values read from config file '+os.path.abspath(config.__file__))
    except ImportError:
        configDict = {'Clip_on':True}
    except Exception as err:
        print(60*'*',"\nError reading config.py file")
        if printInfo:
            import traceback
            print(traceback.format_exc())
        print(60*'*')
        configDict = {'Clip_on':True}

# def MacStartGSASII(g2script,project=''):
#     '''Start a new instance of GSAS-II by opening a new terminal window and starting
#     a new GSAS-II process. Used on Mac OS X only.

#     :param str g2script: file name for the GSASII.py script
#     :param str project: GSAS-II project (.gpx) file to be opened, default is blank
#       which opens a new project
#     '''
#     if project and os.path.splitext(project)[1] != '.gpx':
#         print(f'file {project} cannot be used. Not GSAS-II project (.gpx) file')
#         return
#     if project and not os.path.exists(project):
#         print(f'file {project} cannot be found.')
#         return 
#     elif project:
#         project = os.path.abspath(project)
#         if not os.path.exists(project): 
#             print(f'lost project {project} with abspath')
#             raise Exception(f'lost project {project} with abspath')
#     g2script = os.path.abspath(g2script)
#     pythonapp = sys.executable
#     if os.path.exists(pythonapp+'w'): pythonapp += 'w'
#     script = f'''
# set python to "{pythonapp}"
# set appwithpath to "{g2script}"
# set filename to "{project}"
# set filename to the quoted form of the POSIX path of filename

# tell application "Terminal"
#      activate
#      do script python & " " & appwithpath & " " & filename & "; exit"
# end tell
# '''
#     subprocess.Popen(["osascript","-e",script])

def MacRunScript(script):
    '''Start a bash script in a new terminal window.
    Used on Mac OS X only.

    :param str script: file name for a bash script
    '''
    script = os.path.abspath(script)
    osascript = '''
set bash to "/bin/bash"
set filename to "{}"

tell application "Terminal"
     activate
     do script bash & " " & filename & "; exit"
end tell
'''.format(script)
    subprocess.Popen(["osascript","-e",osascript])

#==============================================================================
#==============================================================================
# conda/pip routines
def findConda():
    '''Determines if GSAS-II has been installed as g2conda or gsas2full
    with conda located relative to this file. 
    We could also look for conda relative to the python (sys.executable)
    image, but I don't want to muck around with python that someone else
    installed.
    '''
    parent = os.path.split(path2GSAS2)[0]
    if sys.platform != "win32":
        activate = os.path.join(parent,'bin','activate')
        conda = os.path.join(parent,'bin','conda')
    else:
        activate = os.path.join(parent,'Scripts','activate.bat')
        conda = os.path.join(parent,'condabin','conda.bat')
    if os.path.exists(activate) and os.path.exists(conda):
        return conda,activate
    else:
        return None
        
def condaTest(requireAPI=False):
    '''Returns True if it appears that Python is being run under Anaconda 
    Python with conda present. Tests for conda environment vars and that 
    the conda package is installed in the current environment.

    :returns: True, if running under Conda
    '''    
    if not all([(i in os.environ) for i in ('CONDA_DEFAULT_ENV','CONDA_EXE', 'CONDA_PREFIX', 'CONDA_PYTHON_EXE')]): return False
    if requireAPI:
        # is the conda package available?
        try:
            import conda.cli.python_api
        except:
            print('You do not have the conda package installed in this environment',
                  '\nConsider using the "conda install conda" command')
            return False

    # There is no foolproof way to check if someone activates conda
    # but then calls a different Python using its path...
    # ...If we are in the base environment then the conda Python 
    # should be the same path as the one currently being run:
    if os.environ['CONDA_DEFAULT_ENV'] == 'base':
        try:
            if os.path.samefile(os.environ['CONDA_PYTHON_EXE'],
                                sys.executable): return True
        except:
            return False

    # ...If not in the base environment, what we can do is check if the
    # python we are running in shares the beginning part of its path with
    # the one in the base installation:
    dir1 = os.path.dirname(os.environ['CONDA_PYTHON_EXE'])
    dir2 = os.path.dirname(sys.executable)
    if sys.platform != "win32": # python in .../bin/..
        dir1 = os.path.dirname(dir1)
        dir2 = os.path.dirname(dir2)    
    return commonPath(dir1,dir2)

def condaInstall(packageList):
    '''Installs one or more packages using the anaconda conda package
    manager. Can be used to install multiple packages and optionally 
    use channels.

    :param list packageList: a list of strings with name(s) of packages
      and optionally conda options. 
      Examples:: 

       packageList=['gsl']
       packageList=['-c','conda-forge','wxpython']
       packageList=['numpy','scipy','matplotlib']

    :returns: None if the the command ran normally, or an error message
      if it did not.
    '''
    try:
        import conda.cli.python_api
    except:
        print('You do not have the conda package installed in this environment',
                  '\nConsider using the "conda install conda" command')
        return None
    try:
        (out, err, rc) = conda.cli.python_api.run_command(
            conda.cli.python_api.Commands.INSTALL,packageList,
            search_path=('conda-forge'),
#    use_exception_handler=True#, stdout=sys.stdout,
            stderr=sys.stderr)
        #print('rc=',rc)
        print('Ran conda. output follows...')
        print(70*'='+'\n'+out+'\n'+70*'=')
        #print('err=',err)
        if rc != 0: return str(out)
    except Exception as msg:
        print(f"\nConda error occurred, see below\n{msg}")
        return "error occurred"
    return None
    
def fullsplit(fil,prev=None):
    '''recursive routine to split all levels of directory names
    '''
    if prev is None: # first call: normalize and drop file name
        fil = os.path.normcase(os.path.abspath(os.path.dirname(fil)))
        prev = []
    i,j = os.path.split(fil)
    if j:
        prev.insert(0,j)
        out = fullsplit(i,prev)
    else:
        return [i]+prev
    return out

def commonPath(dir1,dir2):
    '''Check if two directories share a path. Note that paths 
    are considered the same if either directory is a subdirectory 
    of the other, but not if they are in different subdirectories
    /a/b/c shares a path with /a/b/c/d but /a/b/c/d and /a/b/c/e do not.

    :returns: True if the paths are common
    '''

    for i,j in zip(fullsplit(dir1),fullsplit(dir2)):
        if i != j: return False
    return True

def pipInstall(packageList):
    '''Installs one or more packages using the pip package installer.
    Use of this should be avoided if conda can be used (see :func:`condaTest`
    to test for conda). Can be used to install multiple packages together. 
    One can use pip options, but this is probably not needed. 

    :param list packageList: a list of strings with name(s) of packages
      Examples:: 

       packageList=['gsl']
       packageList=['wxpython','matplotlib','scipy']
       packageList=[r'\\Mac\\Home\\Scratch\\wheels\\pygsl-2.3.3-py3-none-any.whl']
       packageList=['z:/Scratch/wheels/pygsl-2.3.3-py3-none-any.whl']

    :returns: None if the the command ran normally, or an error message
      if it did not.
    '''
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install']+packageList)
    except Exception as msg:
        return msg
    return None
    
def condaEnvCreate(envname, packageList, force=False):
    '''Create a Python interpreter in a new conda environment. Use this
    when there is a potential conflict between packages and it would 
    be better to keep the packages separate (which is one of the reasons
    conda supports environments). Note that conda should be run from the 
    case environment; this attempts to deal with issues if it is not. 

    :param str envname: the name of the environment to be created. 
      If the environment exists, it will be overwritten only if force is True.
    :param list packageList: a list of conda install create command 
      options, such as::

            ['python=3.7', 'conda', 'gsl', 'diffpy.pdffit2',
                '-c', 'conda-forge', '-c', 'diffpy']

    :param bool force: if False (default) an error will be generated 
      if an environment exists

    :returns: (status,msg) where status is True if an error occurs and
      msg is a string with error information if status is True or the 
      location of the newly-created Python interpreter.
    '''
    if not all([(i in os.environ) for i in ('CONDA_DEFAULT_ENV',
                            'CONDA_EXE', 'CONDA_PREFIX', 'CONDA_PYTHON_EXE')]):
        return True,'not running under conda?'
    try:
        import conda.cli.python_api
    except:
        return True,'conda package not available (in environment)'
    # workaround for bug that avoids nesting packages if running from an
    # environment (see https://github.com/conda/conda/issues/11493)
    p = os.path.dirname(os.path.dirname(os.environ['CONDA_EXE']))
    if not os.path.exists(os.path.join(p,'envs')):
        msg = ('Error derived installation path not found: '+
                  os.path.join(p,'envs'))
        print(msg)
        return True,msg
    newenv = os.path.join(p,'envs',envname)
    if os.path.exists(newenv) and not force:
        msg = 'path '+newenv+' already exists and force is not set, aborting'
        print(msg)
        return True,msg
    pathList = ['-p',newenv]
    try:
        (out, err, rc) = conda.cli.python_api.run_command(
            conda.cli.python_api.Commands.CREATE,
            packageList + pathList,
    use_exception_handler=True, stdout=sys.stdout, stderr=sys.stderr
            )
        #print('rc=',rc)
        #print('out=',out)
        #print('err=',err)
        if rc != 0: return True,str(out)
        if sys.platform == "win32":
            newpython = os.path.join(newenv,'python.exe')
        else:
            newpython = os.path.join(newenv,'bin','python')
        if os.path.exists(newpython):
            return False,newpython
        return True,'Unexpected, '+newpython+' not found'
    except Exception as msg:
        print("Error occurred, see below\n",msg)
        return True,'Error: '+str(msg)
    
def addCondaPkg():
    '''Install the conda API into the current conda environment using the 
    command line, so that the API can be used in the current Python interpreter

    Attempts to do this without a shell failed on the Mac because it seems that
    the environment was inherited; seems to work w/o shell on Windows. 
    '''
    if not all([(i in os.environ) for i in ('CONDA_DEFAULT_ENV','CONDA_EXE',
                        'CONDA_PREFIX', 'CONDA_PYTHON_EXE')]):
        return None
    condaexe = os.environ['CONDA_EXE']
    currenv = os.environ['CONDA_DEFAULT_ENV']
    if sys.platform == "win32":
        cmd = [os.environ['CONDA_EXE'],'install','conda','-n',currenv,'-y']
        p = subprocess.Popen(cmd,
                         #stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    else:
        script = 'source ' + os.path.join(
            os.path.dirname(os.environ['CONDA_PYTHON_EXE']),
            'activate') + ' base; '
        script += 'conda install conda -n '+currenv+' -y'
        p = subprocess.Popen(script,shell=True,env={},
                         #stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out,err = MakeByte2str(p.communicate())
    #print('Output from adding conda:\n',out)
    if err:
        print('Note error/warning:')
        print(err)
    if currenv == "base":
        print('\nUnexpected action: adding conda to base environment???')
#==============================================================================
#==============================================================================
# routines for reorg of GSAS-II directory layout
def getIconFile(imgfile):
    '''Looks in either the main GSAS-II install location (old) or subdirectory 
    icons (after reorg) for an icon

    :returns: the full path for the icon file
    '''
    if os.path.exists(os.path.join(path2GSAS2,'icons',imgfile)):
        return os.path.join(path2GSAS2,'icons',imgfile)
    if os.path.exists(os.path.join(path2GSAS2,imgfile)): # patch 3/2024 for svn dir organization
        return os.path.join(path2GSAS2,imgfile)
    print(f'getIconFile Warning: file {imgfile} not found')
    return None

#==============================================================================
#==============================================================================
def makeScriptShortcut():
    '''Creates a shortcut to GSAS-II in the current Python installation
    so that "import G2script" (or "import G2script as GSASIIscripting")
    can be used without having to add GSASII to the path.

    The new shortcut is then tested.

    :returns: returns the name of the created file if successful. None
      indicates an error. 
    '''
    import datetime as dt
    for p in sys.path:
        if 'site-packages' in p: break
    else:
        print('No site-packages directory found in Python path')
        return
    newfil = os.path.join(p,'G2script.py')
    fp = open(newfil,'w')
    fp.write(f'#Created in makeScriptShortcut from {__file__}')
    fp.write(dt.datetime.strftime(dt.datetime.now(),
                                      " at %Y-%m-%dT%H:%M\n"))

    fp.write(f"""import sys,os
Path2GSASII='{path2GSAS2}'
if os.path.exists(os.path.join(Path2GSASII,'GSASIIscriptable.py')):
    print('setting up GSASIIscriptable from',Path2GSASII)
    if Path2GSASII not in sys.path:
        sys.path.insert(0,Path2GSASII)
    from GSASIIscriptable import *
else:
    print('GSASIIscriptable not found in ',Path2GSASII)
    print('Rerun "Install GSASIIscriptable shortcut" from inside GSAS-II')
    sys.exit()
""")
    fp.close()
    print('Created file',newfil)
    try:
        import G2script
    except ImportError:
        print('Unexpected error: import of G2script failed!')
        return
    return newfil

# see if a directory for local modifications is defined. If so, stick that in the path
if os.path.exists(os.path.expanduser('~/.G2local/')):
    sys.path.insert(0,os.path.expanduser('~/.G2local/'))
    fl = glob.glob(os.path.expanduser('~/.G2local/GSASII*.py*'))
    files = ""
    prev = None
    for f in sorted(fl): # make a list of files, dropping .pyc files where a .py exists
        f = os.path.split(f)[1]
        if os.path.splitext(f)[0] == prev: continue
        prev = os.path.splitext(f)[0]
        if files: files += ", "
        files += f
    if files:
        print("*"*75)
        print("Warning: the following source files are locally overridden in "+os.path.expanduser('~/.G2local/'))
        print("  "+files)
        print("*"*75)

BinaryPathFailed = False
BinaryPathLoaded = False
binaryPath = ''
IPyBreak = DoNothing
pdbBreak = DoNothing

if __name__ == '__main__':
    '''What follows is called to update (or downdate) GSAS-II in a 
    separate process. 
    '''
    # check what type of update is being called for
    gitUpdate = False
    preupdateType = None
    updateType = None
    regressversion = None
    help = False
    project = None
    version = None
    
    for arg in sys.argv[1:]:
        if '--git-fetch' in arg:   # pulls latest updates from server but does not apply them
            if preupdateType or updateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            updateType = 'fetch'
        elif '--git-reset' in arg:   # restores locally changed GSAS-II files to distributed versions also updates
            gitUpdate = True
            if preupdateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            preupdateType = 'reset'
        elif '--git-stash' in arg:   # saves locally changed GSAS-II files in "stash"
            argsplit = arg.split('=')
            if len(argsplit) == 1:
                message = None
            elif len(argsplit) == 2:
                message=argsplit[1]
                if message.startswith('"') and message.endswith('"'):
                    message = message.strip('"')
                if message.startswith("'") and message.endswith("'"):
                    message = message.strip("'")
                message = message.replace('"',"'") # double quote not allowed
            else:
                print('invalid form for --git-stash')
                help = True
                break
            gitUpdate = True
            if preupdateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            preupdateType = 'stash'
        elif '--git-update' in arg:  # update to latest downloaded version
            gitUpdate = True
            if updateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            updateType = 'update'
        elif '--git-regress' in arg:
            argsplit = arg.split('=')
            if len(argsplit) != 2:
                print('invalid form for --git-regress')
                help = True
                break
            gitversion = argsplit[1]
            # make sure the version or tag supplied is valid and convert to
            # a full sha hash
            g2repo = openGitRepo(path2GSAS2)
            try:
                regressversion = g2repo.commit(gitversion).hexsha
            except git.BadName:
                print(f'invalid version specified ({version}) for GitHub regression')
                help = True
                break
            if updateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            updateType = 'regress'
            gitUpdate = True
        elif '--help' in arg:
            help = True
            break
        elif os.path.exists(arg):   # svn args parsed later; this is just checking
            project = arg
            pass
        else:   # for old-style svn update
            if arg.isdecimal() or not arg: 
                #version = arg
                pass
            else:
                print(f'unknown arg {arg}')
                help = True
        if gitUpdate and version:
            print('Conflicting arguments (git & svn opts combined?)')
            help = True

    if help or len(sys.argv) == 1:
        print('''Options when running GSASIIpath.py standalone

to update/regress repository from svn repository:
   python GSASIIpath.py <project> <version>
            where <project> is an optional path reference to a .gpx file
            and <version> is a specific GSAS-II version to install 
                (default is latest)

to update/regress repository from git repository:
   python GSASIIpath.py option <project>
       where option will be one or more of the following:
            --git-fetch            downloads lastest changes from repo
                                   any other options will be ignored

            --git-stash="message"  saves local changes 

            --git-reset            discards local changes 

            --git-update

            --git-regress=version

       and where <project> is an optional path reference to a .gpx file

       Note: --git-reset and --git-stash cannot be used together. Likewise
             --git-update and --git-regress cannot be used together.
             However either --git-reset or --git-stash can be used 
             with either --git-update or --git-regress.
             --git-fetch cannot be used with any other options.
''')
        sys.exit()


    if updateType == 'fetch':
        # download the latest updates from GitHub to the local repository
        # in background while GSAS-II runs no updates are applied
        logfile = os.path.join(os.path.expanduser('~'),'GSASII_bkgUpdate.log')
        mode = 'a'
        # don't let log file get too large (20K bytes)
        if os.path.exists(logfile) and os.path.getsize(logfile) > 20000:
            mode = 'w'
        # if file open fails, there is probably a concurent update process
        try:
            fp = open(logfile,mode)
        except:
            print('background git update was unable to open log file')
            sys.exit()
        fp.write('Starting background git update')
        fp.write(dt.datetime.strftime(dt.datetime.now(),
                                      " at %Y-%m-%dT%H:%M\n"))
        try:
            import git
        except:
            fp.write('git import failed')
            fp.close()
            sys.exit()
        try:
            g2repo = openGitRepo(path2GSAS2)
            g2repo.remote().fetch()
            fp.write(f'Updates fetched\n')
        except Exception as msg:
            fp.write(f'Update failed with message {msg}\n')

        if g2repo.head.is_detached:
            fp.write(f'Status: reverted to an old install\n')
        else:
            try:
                rc,lc,_ = gitCheckForUpdates(False,g2repo)
                if len(rc) == 0:
                    fp.write('Status: no unapplied commits\n')
                else:
                    fp.write(f'Status: unapplied commits now {len(rc)}\n')
            except Exception as msg:
                fp.write(f'\ngitCheckForUpdates failed with message {msg}\n')
        fp.write('update done at')
        fp.write(dt.datetime.strftime(dt.datetime.now(),
                                      " at %Y-%m-%dT%H:%M\n\n"))
        fp.close()
        sys.exit()
        
    if gitUpdate:
        import time
        time.sleep(1) # delay to give the main process a chance to exit
                      # so we don't change code for a running process
                      # windows does not like that
        try:
            import git
        except:
            print('git import failed')
            sys.exit()
        try:
            g2repo = openGitRepo(path2GSAS2)
        except Exception as msg:
            print(f'Update failed with message {msg}\n')
            sys.exit()
        print('git repo opened')
                  
    if preupdateType == 'reset':
        # --git-reset   (preupdateType = 'reset')
        print('Restoring locally-updated GSAS-II files to original status')
        openGitRepo(path2GSAS2).git.reset('--hard','origin/master')
        try:
            if g2repo.active_branch.name != 'master': 
                g2repo.git.switch('master')
        except TypeError:   # fails on detached head
            pass

    elif preupdateType == 'stash':
        # --git-stash   (preupdateType = 'stash')
        print('Stashing locally-updated GSAS-II files')
        if message:
            g2repo.git.stash(f'-m"{message}"')
        else:
            g2repo.git.stash()
            
    # Update to the latest GSAS-II version. This assumes that a fetch has
    # been done prior, or this will only update to the last time that
    # it was done.
    if updateType == 'update':
        # --git-update  (updateType = 'update')
        if g2repo.is_dirty():
            print('Cannot update a directory with locally-made changes')
            sys.exit()
        print('Updating to latest GSAS-II version')
        if g2repo.head.is_detached:
            g2repo.git.switch('master')
        g2repo.git.merge('--ff-only')
        print('git: updated to latest version')

    # Update or regress to a specific GSAS-II version.
    # this will always cause a "detached head" status
    elif updateType == 'regress':
        # --git-regress (updateType = 'regress')
        if g2repo.is_dirty():
            print('Cannot regress a directory with locally-made changes')
            sys.exit()
        print(f'Regressing to git version {regressversion[:6]}')
        g2repo.git.checkout(regressversion)

    if gitUpdate:
        # now restart GSAS-II with the new version
        # G2scrpt = os.path.join(path2GSAS2,'GSASII.py')
        if project:
            print(f"Restart GSAS-II with project file {project!r}")
            # subprocess.Popen([sys.executable,G2scrpt,project])
        else:
            print("Restart GSAS-II without a project file ")
            # subprocess.Popen([sys.executable,G2scrpt])
        import GSASIIfiles
        GSASIIfiles.openInNewTerm(project)
        print ('exiting update process')
        sys.exit()
        
    else:
        # this is the old svn update process
        LoadConfig()
        import time
        time.sleep(1) # delay to give the main process a chance to exit
        # perform an update and restart GSAS-II
        try:
            project,version = sys.argv[1:3]
        except ValueError:
            project = None
            version = 'trunk'
        loc = os.path.dirname(__file__)
        if version == 'trunk':
            svnSwitch2branch('')
        elif '/' in version:
            svnSwitch2branch(version)
        elif version:
            print("Regress to version "+str(version))
            svnUpdateDir(loc,version=version)
        else:
            print("Update to current version")
            svnUpdateDir(loc)
        ex = sys.executable
        if sys.platform == "darwin": # mac requires pythonw which is not always reported as sys.executable
            if os.path.exists(ex+'w'): ex += 'w'
        if project:
            print("Restart GSAS-II with project file "+str(project))
            subprocess.Popen([ex,os.path.join(loc,'GSASII.py'),project])
        else:
            print("Restart GSAS-II without a project file ")
            subprocess.Popen([ex,os.path.join(loc,'GSASII.py')])
        print ('exiting update process')
        sys.exit()
