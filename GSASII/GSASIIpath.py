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
#import pathlib
import subprocess
import datetime as dt
try:
    import numpy as np
except ImportError:
    print("numpy import failed in GSASIIpath")

# fix up path before using git. Needed when using conda without
#   activate (happens on MacOS w/GSAS-II.app)
pyPath = os.path.dirname(os.path.realpath(sys.executable))
if sys.platform != "win32" and pyPath not in os.environ['PATH'].split(':'):
    os.environ['PATH'] = pyPath + ':' + os.environ['PATH']

# hard-coded github repo locations
G2binURL = "https://api.github.com/repos/AdvancedPhotonSource/GSAS-II-buildtools"
g2URL = "https://github.com/AdvancedPhotonSource/GSAS-II.git"
# tutorial repo owner & Repo name
gitTutorialOwn,gitTutorialRepo = 'AdvancedPhotonSource', 'GSAS-II-Tutorials'

path2GSAS2 = os.path.dirname(os.path.abspath(os.path.expanduser(__file__))) # location of this file; save before any changes in pwd

# convert version numbers as '1.2.3' to integers (1002) and back (to 1.2)
fmtver = lambda v: str(v//1000)+'.'+str(v%1000)
intver = lambda vs: sum([int(i) for i in vs.split('.')[0:2]]*np.array((1000,1)))

def GetConfigValue(key,default=None,getDefault=False):
    '''Return the configuration file value for key or a default value
    if not specified.

    :param str key: a value to be found in the configuration settings
    :param any default: a value to be supplied if a value for key is
      not specified in the config file or the config file is not found.
      Defaults to None.
    :param bool getDefault: If True looks up the default value from the
      config_example.py file (default value is False). Do not specify a
      getDefault=True if a value is provided for default.
    :returns: the value found or the default.
    '''
    if getDefault:
        if default is not None:
            raise ValueError('Cannot use default and getDefault together')
        default = GetConfigDefault(key)
    try:
        return configDict.get(key,default)
    except NameError: # this happens when building docs
        return None

def SetConfigValue(parmdict):
    '''Set configuration variables. Note that parmdict is a dictionary
    from :func:`GSASIIctrlGUI.GetConfigValsDocs` where each element is a
    lists. The first item in list is the default value, the second is
    the value to use for that configuration variable. Most of the
    information gathered in GetConfigValsDocs is no longer used.
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

def GetConfigDefault(key):
    '''Return the default value for a config value

    :param str key: a value to be found in the configuration (config_example.py) file
    :returns: the default value or None
    '''
    try:
        from . import config_example
    except:
        return None
    return config_example.__dict__.get(key)

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
    fsource/SConstruct or GSASII-buildtools/compile/nameTar.py
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
G2_installed_result = None
def HowIsG2Installed():
    '''Determines if GSAS-II was installed with git.
    Result is cached to avoid time needed for multiple calls of this.

    :returns:
      * a string starting with 'git' from git,
        if installed from the GSAS-II GitHub repository (defined in g2URL),
        the string is 'github', if the post-3/2024 directory structure is
        in use '-rev' is appended.
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
        elif g2URL.replace('https://github.com/',
                           'git@github.com:') in g2repo.remote().urls:
            return 'github'+rev
        G2_installed_result = 'git'+rev
        return G2_installed_result
    except:
        pass
    G2_installed_result = 'noVCS'
    return G2_installed_result

def getSavedVersionInfo():
    '''Get version number information from a file written by install
    routines. This is faster than getting the information from git. Also,
    when GSAS-II is installed into Python, the files are no longer in 
    a git repository so querying git is not possible.

    The saved_version.py file is written by install/save_versions.py.
    The git_verinfo.py file is written by install/tag-version.py or
    by install/incr-mini-version.py. If both are present, use the 
    saved_version.py file preferentially. 

    :returns: a reference to the version variables or None if no 
      version info file is found.
    '''
    try:
        from . import saved_version as gv
        return gv
    except:
        pass
    try:
        from . import git_verinfo as gv
        return gv
    except:
        pass
    return None  # this is unexpected as all dists should have git_verinfo.py

def GetVersionNumber():
    '''Obtain a numeric (sequential) version number for GSAS-II from version
    files, or directly from git if no other choice.

    :returns: an int value normally, but unexpected error conditions could result in 
      a value of 'unknown' or '?'.
    '''
    # look for a previously recorded tag -- this is quick
    gv = getSavedVersionInfo()
    if gv is not None:
        for item in gv.git_tags+gv.git_prevtags:
            if item.isnumeric(): return int(item)

    if HowIsG2Installed().startswith('git'):
        # unexpected: should always find a version from getSavedVersionInfo()
        #   from a version file, but if that fails ask Git for the most
        #   recent tag that starts & ends with a digit and does not contain a '.'
        try:
            g2repo = openGitRepo(path2GSAS2)
            tag,vers,gnum = g2repo.git.describe('--tags','--match','[0-9]*[0-9]',
                                                         '--exclude','*.*').split('-')
            if tag.isnumeric(): return int(tag)
        except:
            pass
        return "unknown"

    # Not installed from git and no version file -- very strange!
    return "?"

def GetVersionTag():
    '''Obtain a release (X.Y.Z) version number for GSAS-II from version
    files, or directly from git if no other choice.

    :returns: a string of form <Major>.<Minor>.<mini> normally, but unexpected 
      error conditions could result in a value of '?.?.?' or '?'.
    '''
    # look for a previously recorded tag -- this is quick
    gv = getSavedVersionInfo()
    if gv is not None:
        if '?' not in gv.git_versiontag: return gv.git_versiontag

    if HowIsG2Installed().startswith('git'):
        # unexpected: should always find a version from getSavedVersionInfo()
        #   from a version file, but if that fails ask Git for the most
        #   recent tag that has the X.Y.Z format.
        try:
            g2repo = openGitRepo(path2GSAS2)
            tag,vers,gnum = g2repo.git.describe('--tags','--match','*.*.*').split('-')
            return tag
        except:
            pass
        return "?.?.?"

    # Not installed from git and no version file -- very strange!
    return "?"

def getG2Branch():
    '''Get name of current branch, as named on local computer
    '''
    if HowIsG2Installed().startswith('git'):
        g2repo = openGitRepo(path2GSAS2)
        try:
            return g2repo.active_branch.name
        except TypeError: # likely on a detached head
            return '?'
    
def getG2VersionInfo():
    '''Get the git version information. This can be a bit slow, so reading
    .../GSASII/saved_version.py may be faster (in main but not master branch)
    '''
    gv = getSavedVersionInfo()
    if HowIsG2Installed().startswith('git'):
        g2repo = openGitRepo(path2GSAS2)
        commit = g2repo.head.commit
        ctim = commit.committed_datetime.strftime('%d-%b-%Y %H:%M')
        now = dt.datetime.now().replace(
            tzinfo=commit.committed_datetime.tzinfo)
        delta = now - commit.committed_datetime
        age = delta.total_seconds()/(60*60*24.)
        gversion = f"Tag: #{GetVersionNumber()}, {GetVersionTag()}"
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
            elif str(g2repo.active_branch) != 'main':
                msg += f'\n\tUsing development branch "{g2repo.active_branch}"'
            elif age > 60 and len(rc) > 0:
                msg += f"\n\t**** This version is really old. Please update. >= {len(rc)} updates have been posted ****"
            elif age > 5 and len(rc) > 0:
                msg += f"\n\t**** Please consider updating. >= {len(rc)} updates have been posted"
            elif len(rc) > 0:
                msg += f"\n\tThis GSAS-II version is ~{len(rc)} updates behind current."
        return f"  GSAS-II:    {commit.hexsha[:8]}, {ctim} ({age:.1f} days old). {gversion}{msg}"
    elif gv is not None:
        vt = ''
        cvt = ''
        try:
            if gv.git_versiontag:
                vt = gv.git_versiontag
                cvt,cvn = getGitHubVersion()
        except:
            pass
                
        for item in gv.git_tags+gv.git_prevtags:
            if item.isnumeric():
                tags = item
                if vt:
                    tags += f', {vt}'
                msg = f"GSAS-II version: Git: {gv.git_version[:8]}, #{tags} (reinstall to update)"
                if vt != cvt and cvt is not None:
                    msg += f'\n\tNote that the current GSAS-II version is {cvt}'
                return msg
    # Failed to get version info, fallback on old version number routine; should not happen anymore
    return f"GSAS-II installed without git, last tag: #{GetVersionNumber()}, {GetVersionTag()}"

#==============================================================================
#==============================================================================
# routines to interface with GitHub.
def saveGitHubVersion():
    '''Get the latest GSAS-II version tags from the GitHub site
    and place them into the config.ini file. This is always done in 
    background so that app startup time is minimally delayed. 

    :returns: Returns a Popen object (see subprocess).
    '''
    try:
        import requests
    except:
        print('Unable to use requests module')
        return
    return subprocess.Popen([sys.executable, __file__, '--github-tags'])
    if GetConfigValue('debug'): print('Updates fetched')

def getGitHubVersion():
    '''Get the latest git version info when not accessible from git
    as saved by saveGitHubVersion
    '''
    import configparser
    cfgfile = os.path.expanduser(os.path.normpath('~/.GSASII/config.ini'))
    if not os.path.exists(cfgfile):
        if GetConfigValue('debug'): print(f"{cfgfile} not found")
        return None,None
    try:
        cfg = configparser.ConfigParser()
        # Read the configuration file
        cfg.read(cfgfile)
    except Exception as err:
        if GetConfigValue('debug'): print(f"Error reading {cfgfile}\n",err)
        return None,None
    if 'version info' not in cfg:
        if GetConfigValue('debug'): print(f"no saved version number in {cfgfile}")
        return None,None
    return cfg['version info'].get('lastVersionTag'),cfg['version info'].get('lastVersionNumber')

# routines to interface with git.
BASE_HEADER = {'Accept': 'application/vnd.github+json',
               'X-GitHub-Api-Version': '2022-11-28'}

def openGitRepo(repo_path):
    try:
        import git
    except:
        return None
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
    try:
        import git
    except:
        return None
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
       * value&8==8: repository has has been switched to other than main branch

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
        #count_modified_files = len(g2repo.index.diff(None))
    if g2repo.head.is_detached:
        code += 2                             # detached
    else:
        if g2repo.active_branch.name != 'main':
            code += 8                         # not on main branch
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
    try:
        import git
    except:
        print('Failed to import git in gitCheckForUpdates()')
        return (None,None,None)
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
            print(f'Failed to get updates from {g2repo.remote().url}\nerror: {msg}')
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
    a commit that is containined in the main branch

    returns the count and the commit object for the
    parent commit that connects the current stranded
    branch to the main branch.

    None is returned if no connection is found
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    if not g2repo.head.is_detached:
        return 0,g2repo.commit()
    # is detached head in main branch?
    if g2repo.commit() in g2repo.iter_commits('main'):
        return 0,g2repo.commit()
    # count number of commits since leaving main branch
    masterList = list(g2repo.iter_commits('main'))
    for c,i in enumerate(g2repo.commit().iter_parents()):
        if i in masterList:
            return c+1,i
    else:
        return None,None

def gitCountRegressions(g2repo=None):
    '''Count the number of new check ins on the main branch since
    the head was detached as well as any checkins made on the detached
    head.

    :returns: maincount,detachedcount, where

      * maincount is the number of check ins made on the main branch
        remote repository since the reverted check in was first made.
      * detachedcount is the number of check ins made locally
        starting from the detached head (hopefully 0)

      If the connection between the current head and the main branch
      cannot be established, None is returned for both.
      If the connection from the reverted check in to the newest version
      (I don't see how this could happen) then only maincount will be None.
    '''
    if g2repo is None:
        g2repo = openGitRepo(path2GSAS2)
    # get parent of current head that is in main branch
    detachedcount,parent = countDetachedCommits(g2repo)
    if detachedcount is None: return None,None
    maincount = 0
    for h in g2repo.iter_commits('main'):
        if h == parent:
            return maincount,detachedcount
        maincount += 1
    return None,detachedcount

def gitGetUpdate(mode='Background'):
    '''Download the latest updates into the local copy of the GSAS-II
    repository from the remote main, but don't actually update the
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
        if GetConfigValue('debug'): print('Updates fetched')

def gitHistory(values='tag',g2repo=None,maxdepth=100):
    '''Provides the history of commits to the main, either as tags
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
    history = list(g2repo.iter_commits('main'))
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
        raise ValueError(f'gitHistory has invalid values specified: {values}')

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
    try:
        import requests
    except:
        print('Unable to install binaries in getGitBinaryReleases():\n requests module not available')
        return
    # Get first page of releases. (Could there be more than one?)
    releases = []
    tries = 0
    while tries < 5: # this has been known to fail, so retry
        tries += 1
        releases = requests.get(
            url=f"{G2binURL}/releases",
            headers=BASE_HEADER
        ).json()
        try:
            # loop over assets of latest release (will [-1] always get this?)
            versions = []
            URLs = []
            for asset in releases[-1]['assets']:
                if not asset['name'].endswith('.tgz'): continue
                versions.append(asset['name'][:-4]) # Remove .tgz tail
                URLs.append(asset['browser_download_url'])
            count = len(versions)
            # Cache the binary releases for later use in case GitHub
            # prevents us from using a query to get them
            if cache and count > 4:
                fp = open(os.path.join(path2GSAS2,'inputs','BinariesCache.txt'),'w')
                res = dict(zip(versions,URLs))
                for key in res:
                    fp.write(f'{key} : {res[key]}\n')
                fp.close()
            return dict(zip(versions,URLs))
        except:
            print('Attempt to get GSAS-II binary releases/assets failed, sleeping for 10 sec and then retrying')
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

def getGitBinaryLoc(npver=None,pyver=None,verbose=True,debug=False):
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
    if debug:
        print('URLdict:')
        for k in URLdict: print(k,URLdict[k])
    versions = {}
    for d in URLdict:
        if d.startswith(bindir):
            v = intver(d.rstrip('/').split('_')[3].lstrip('n'))
            versions[v] = d
    if debug: print('versions:',versions)
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
    import tempfile
    import tarfile
    try:
        import requests
    except:
        print('Unable to install binaries in InstallGitBinary():\n requests module not available')
        return
    # download to scratch
    tarobj = None
    tar = tempfile.NamedTemporaryFile(suffix='.tgz',delete=False)
    try:
        tar.close()
        if verbose: print(f'Downloading {tarURL}')
        r = requests.get(tarURL, allow_redirects=True)
        with open(tar.name, 'wb') as fp:
            fp.write(r.content)
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
            if verbose: print(f'Created GSAS-II binary file {os.path.split(newfil)[1]}')
        if verbose: print(f'Binary files created in {os.path.split(newfil)[0]}')

    finally:
        if tarobj: del tarobj
        os.unlink(tar.name)

def GetRepoUpdatesInBackground():
    '''Get the latest GSAS-II version info.
    This serves to make sure that :func:`gitGetUpdate` is called only
    if git has been used to install GSAS-II.

    :returns: returns a Popen object (see subprocess)
    '''
    if HowIsG2Installed().startswith('git'):
        return gitGetUpdate(mode='Background')
    else:
        return saveGitHubVersion()

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
    try:
        import requests
    except:
        print('Unable to search GitHub in dirGitHub():\n requests module not available')
        return
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
                 branchname="main"):
    '''Create a URL that can be used to view/downlaod the raw version of
    file in a GitHub repository.

    :param str dirlist: a list of sub-directories `['parent','child',sub']`
      for `parent/child/sub` or `[]` for a file in the top-level directory.
    :param str filename: the name of the file
    :param str orgName: the name of the GitHub organization
    :param str repoName: the name of the GitHub repository
    :param str branchname: the name of the GitHub branch. Defaults
       to "main".

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
    try:
        import requests
    except:
        print('Unable to download Tutorial data in downloadDirContents():\n requests module not available')
        return
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
            with open(outfil, 'wb') as fp:
                fp.write(r.content)
            print(f'wrote {outfil}')
        elif r.status_code == 404:
            print(f'Warning: {fil} is likely a subdirectory of directory {"/".join(dirlist)!r}')
        else:
            print(f'Unexpected web response for {fil}: {r.status_code}')
    return

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
            main = os.path.join(path2GSAS2,'G2.py')
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
    This routine is only used when debug=True is set in the configuration
    settings
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

    This routine is only used when debug=True is set in the configuration settings
    '''
    try:
        #from IPython.core import ultratb
        import IPython.core.ultratb
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
    except TypeError:  # Ipython 9.x removes color_scheme
        try:
            IPython.core.ultratb.FormattedTB(call_pdb=False)(*args)
            from IPython.core import getipython
            if getipython.get_ipython() is None:
                ipshell = InteractiveShellEmbed.instance()
            else:
                ipshell = InteractiveShellEmbed()
        except Exception as msg:
            print('IPython patch failed, msg=',msg)
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
    except when the debug option is set True in the configuration settings
    '''
    pass

def InvokeDebugOpts():
    'Called to set up debug options'
    if any('SPYDER' in name for name in os.environ):
        print('Running from Spyder, keeping breakpoint() active & skipping exception trapping')
    elif GetConfigValue('debug'):
        try:
            import pdb
            global pdbBreak
            pdbBreak = pdb.set_trace
            import IPython
            IPython
            global IPyBreak
            IPyBreak = IPyBreak_base
            sys.excepthook = exceptHook
            os.environ['PYTHONBREAKPOINT'] = 'GSASIIpath.IPyBreak_base'
            print ('Debug on: IPython: Exceptions and G2path.IPyBreak(); pdb: G2path.pdbBreak()')
        except:
            print ('Debug on failed. IPython not installed?')
    else: # not in spyder or debug enabled, hide breakpoints
        os.environ['PYTHONBREAKPOINT'] = '0'

def TestSPG():
    '''Test if pyspg.[so,.pyd] can be run from a location in the existing path
    Do not modify the path if not.
    '''
    def showVersion():
        try:
            f = os.path.join(os.path.dirname(pyspg.__file__),'GSASIIversion.txt')
            with open(f,'r') as fp:
                version = fp.readline().strip()
                vnum = fp.readline().strip()
            print(f'  Binary ver: {vnum}, {version}')
        except:
            if GetConfigValue('debug'): 
                print('  Binaries:   undated')
    try:
        from . import pyspg
        pyspg
        showVersion()
        return True
    except ImportError:
        pass
    try:
        import pyspg
        pyspg
    except ImportError:
        return False
    try:
        pyspg.sgforpy('P -1')
    except Exception as err:
        print(70*'=')
        print(f'Module pyspg in {pyspg.__file__} could not be run\nerror msg: {err}')
        print(70*'=')
        return False
    showVersion()
    return True

def pathhack_TestSPG(fpth):
    '''Test if pyspg.[so,.pyd] can be run from a specified location. If so 
    modify the path to include it.
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

def SetBinaryPath(showConfigMsg=False):
    '''
    Add location of GSAS-II shared libraries (binaries: .so or
    .pyd files) to path (when needed). When GSAS-II is installed by
    pixi, no change in the path is needed.

    This routine must be executed after GSASIIpath is imported
    and before any other GSAS-II imports are done, since
    they may assume binary files are in path

    :param bool showConfigMsg: When True, config info is shown (default is False)
    '''
    # cache the results of this routine so that repeated calls
    # only search for binary routines once
    global BinaryPathLoaded,binaryPath,BinaryPathFailed
    if BinaryPathLoaded or BinaryPathFailed: return
    try:
        intver(np.__version__)
    except (AttributeError,TypeError): # happens on building docs
        return
    try:
        from GSASII import pypowder
        pypowder
        binaryPath = None   # special value to indicate that binaries have been installed into package
        if showConfigMsg:
            print(f'GSAS-II binaries co-located with GSAS-II: {os.path.dirname(__file__)}')
        BinaryPathLoaded = True
        LoadConfig(showConfigMsg)
        return
    except ImportError:
        pass

    try:
        from . import pathHacking
    except ImportError:
        print('Binary load failed and module pathHacking not present')
        BinaryPathFailed = True
        return

    LoadConfig(showConfigMsg)
    BinaryPathFailed = pathHacking._path_discovery(showConfigMsg)

def WriteConfig(configDict):
    '''Write the configDict information to the GSAS-II ini settings
    into file ~/.GSASII/config.ini. Called from
    :func:`GSASIIctrlGUI.SaveConfigVars`.
    '''
    import configparser

    localdir = os.path.expanduser(os.path.normpath('~/.GSASII'))
    if not os.path.exists(localdir):
        try:
            os.mkdir(localdir)
            print(f'Created directory {localdir}')
        except Exception as msg:
            print(f'Error trying to create directory {localdir}\n{msg}')
            return True
    cfgfile = os.path.join(localdir,'config.ini')
    cfgP = configparser.ConfigParser()
    if os.path.exists(cfgfile): 
        cfgP.read(cfgfile)  # read previous file so other sections are retained
    cfgP['GUI settings'] = configDict

    # Write the configuration file
    with open(cfgfile, 'w') as configfile:
        cfgP.write(configfile)
    print(f"Configuration settings saved as {cfgfile}")

def LoadConfig(printInfo=True):
    '''Read configuration settings from ~/.GSASII/config.ini, if present.
    Place the values into global dict configDict.

    :param bool printInfo: if printInfo is True (default) then a message
      is shown with the number of settings read (upon startup).
    '''
    def XferConfigIni():
        '''copy the contents of the config.py file to file ~/.GSASII/config.ini.
        This "patch code" used for master->main transition and can eventually
        be removed.
        '''
        import types
        configDict = {}
        try:
            import config
            #import config_example as config
            for i in config.__dict__:
                if i.startswith('__') and i.endswith('__'): continue
                if isinstance(config.__dict__[i],types.ModuleType): continue
                configDict.update({i:str(config.__dict__[i])})
        except ImportError as err:
            print("New install: start without a config.py file")
            return
        except Exception as err:
            print("Error reading config.py file\n",err)
            return
        print(f"Contents of {config.__file__} to be written from config.py...")
        WriteConfig(configDict)

    import configparser
    global configDict
    configDict = {}
    cfgfile = os.path.expanduser(os.path.normpath('~/.GSASII/config.ini'))
    if not os.path.exists(cfgfile):
        print(f'N.B. Configuration file {cfgfile} does not exist')
        # patch 2/7/25: transform GSAS-II config.py contents to config.ini
        XferConfigIni()
    try:
        from . import config_example
    except ImportError as err:
        try:
            import GSASII.config_example as config_example
        except ImportError as err:
            print("Error importing config_example.py file\n",err)
            return

    # get the original capitalization (lost by configparser)
    capsDict = {key.lower():key for key in config_example.__dict__ if not key.startswith('__')}

    try:
        cfg = configparser.ConfigParser()
        # Read the configuration file
        cfg.read(cfgfile)
    except Exception as err:
        print("Error reading {cfgfile}\n",err)
        return

    # Access values from the configuration file
    try:
        cfgG = cfg['GUI settings']
    except KeyError:
        cfgG = {}
        
    for key in cfgG:
        key = key.lower()  # not needed... but in case configparser ever changes
        capKey = capsDict.get(key)
        if capKey is None:
            print(f'Item {key} not defined in config_example')
            continue
        try:
            if cfgG[key] == 'None':
                configDict[capKey] = None
            elif key.endswith('_pos') or key.endswith('_size'): # list of integers
                configDict[capKey] = tuple([int(i) for i in
                                    cfgG[key].strip('()').split(',')])
            elif key.endswith('_location') or key.endswith('_directory') or key.endswith('_exec'): # None (above) or str
                configDict[capKey] = cfgG.get(key)
            elif cfgG[key].startswith('[') and cfgG[key].endswith(']'): # list of strings
                s = cfgG[key].strip('[]')
                if s == '':
                    res = []
                else:
                    res = [i.strip("'").replace(r'\\','\\') for i in s.split(', ')]
                configDict[capKey] = res
            elif isinstance(config_example.__dict__[capKey],bool):
                configDict[capKey] = cfgG.getboolean(key)
            elif isinstance(config_example.__dict__[capKey],float):
                configDict[capKey] = cfgG.getfloat(key)
            elif isinstance(config_example.__dict__[capKey],int):
                configDict[capKey] = cfgG.getint(key)
            elif isinstance(config_example.__dict__[capKey],str):
                configDict[capKey] = cfgG.get(key).replace(r'\\','\\')
            else:
                print('*** problem with',type(config_example.__dict__[capKey]))
                continue
        except:
            continue
    if printInfo:
        print (f'{len(configDict)} values read from {cfgfile}')
    # make sure this value is set
    configDict['Clip_on'] = configDict.get('Clip_on',True)

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
            conda.cli.python_api
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
        print(f'Preparing to install package(s): {" ,".join(packageList)}'+
                  '\nThis can take a while')
        # the next line works, but the subsequent cli is considered more stable
        #conda.cli.main('install',  '-y', *packageList)
        # this is considered to be supported in the long term
        (out, err, rc) = conda.cli.python_api.run_command(
            conda.cli.python_api.Commands.INSTALL,packageList,
#            search_path=('conda-forge'),   # broken!
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
    # update conda package specs (pkg=1.2.3) to pip package specs (pkg==1.2.3)
    # at present no other specifiers are used. 
    for i,val in enumerate(packageList):
        if '=' in val and '==' not in val:
            packageList[i] = packageList[i].replace('=','==')
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
    base environment; this attempts to deal with issues if it is not.

    Currently, this is used only to install diffpy.PDFfit2.

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
        p = sys.exec_prefix
    else:
        # workaround for bug that avoids nesting packages if running from an
        # environment (see https://github.com/conda/conda/issues/11493)
        p = os.path.dirname(os.path.dirname(os.environ['CONDA_EXE']))
    try:
        import conda.cli.python_api
    except:
        return True,'conda package not available (in environment)'
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
            use_exception_handler=True) # ,stdout=sys.stdout, stderr=sys.stderr)
        print(out)
        if rc != 0:
            print(err)
            return True,str(out+err)
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
    #condaexe = os.environ['CONDA_EXE']
    currenv = os.environ['CONDA_DEFAULT_ENV']
    if sys.platform == "win32":
        cmd = [os.environ['CONDA_EXE'],'install','conda','-n',currenv,'-y']
        with subprocess.Popen(cmd,
                         #stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         encoding='UTF-8') as p:
            out,err = p.communicate()
    else:
        script = 'source ' + os.path.join(
            os.path.dirname(os.environ['CONDA_PYTHON_EXE']),
            'activate') + ' base; '
        script += 'conda install conda -n '+currenv+' -y'
        with subprocess.Popen(script,shell=True,env={},
                         #stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         encoding='UTF-8') as p:
            out,err = p.communicate()
    
    if out is not None and GetConfigValue('debug'): print('Output from adding conda:\n',out)
    if err and err is not None:
        print('Note error/warning from running conda:\n',err)
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
    if not HowIsG2Installed().startswith('git'):
        print('GSAS-II installed directly, shortcut likely not needed')
        return None
    for p in sys.path:
        if 'site-packages' in p: break
    else:
        print('No site-packages directory found in Python path')
        return
    newfil = os.path.join(p,'G2script.py')
    with open(newfil,'w') as fp:
        fp.write(f'#Created in makeScriptShortcut from {__file__}')
        fp.write(dt.datetime.strftime(dt.datetime.now(),
                                          " at %Y-%m-%dT%H:%M\n"))

        fp.write(f"""
import sys,os
Path2GSASII=r'{path2GSAS2}'
if os.path.exists(os.path.join(Path2GSASII,'GSASIIscriptable.py')):
    print('setting up GSASIIscriptable from',Path2GSASII)
    if os.path.dirname(Path2GSASII) not in sys.path:
        sys.path.insert(0,os.path.dirname(Path2GSASII))
    try:
        from GSASII.GSASIIscriptable import *
    except:
        print('Import of GSASIIscriptable failed.\\nRerun "Install GSASIIscriptable shortcut" from inside GSAS-II?')
        sys.exit()
else:
    print('GSASIIscriptable not found in ',Path2GSASII)
    print('Rerun "Install GSASIIscriptable shortcut" from inside GSAS-II')
    sys.exit()
    """)
        fp.close()
    print('Created file',newfil)
    try:
        import G2script
        G2script
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

def postURL(URL,postdict,getcookie=None,usecookie=None,
                timeout=None,retry=2,mode='get'):
    '''Posts a set of values as from a web form using the "get" or "post"
    protocols.
    If access fails to an https site, the access is retried with http.

    :param str URL: the URL to post; typically something
       like 'http://www.../dir/page?'
    :param dict postdict: contains keywords and values, such
       as {'centrosymmetry': '0', 'crystalsystem': '0', ...}
    :param dict getcookie: dict to save cookies created in call, or None
       (default) if not needed.
    :param dict usecookie: dict containing cookies to be used in call,
       or None (default) if not needed.
    :param int timeout: specifies a timeout period for the get or post (default
      is None, which means the timeout period is set by the server). The value
      when specified is the time in seconds to wait before giving up on the
      request.
    :param int retry: the number of times to retry the request, if it times out.
      This is only used if timeout is specified. The default is 2. Note that
      if retry is left at the default value (2), The timeout is increased by
      25% for the second try.
    :param str mode: either 'get' (default) or 'post'. Determines how
       the request will be submitted.
    :returns: a string with the response from the web server or None
       if access fails.
    '''
    try:
        import requests # delay this until now, since rarely needed
    except:
        # this import seems to fail with the Anaconda pythonw on
        # macs; it should not!
        print('Warning: failed to import requests. Python config error')
        return None

    if mode == 'get':
        reqopt = requests.get
    else:
        reqopt = requests.post

    repeat = True
    count = 0
    while repeat:
        count += 1
        r = None
        repeat = False
        try:
            if timeout is not None:
                r = reqopt(URL,params=postdict,cookies=usecookie,
                                     timeout=timeout)
            else:
                r = reqopt(URL,params=postdict,cookies=usecookie)
            if r.status_code == 200:
                if GetConfigValue('debug'): print('request OK')
                page = r.text
                if getcookie is not None:
                    getcookie.update(r.cookies)
                return page # success
            else:
                print('request to {} failed. Reason={}'.format(URL,r.reason))
        except requests.exceptions.ConnectionError as msg:
            if 'time' in str(msg) and 'out' in str(msg):
                print(f'server timeout accessing {URL}')
                if GetConfigValue('debug'): print('full error=',msg)
                if timeout is not None and count < retry:
                    if retry == 2:
                        timeout *= 1.25
                        print(f'retry with timout={timeout} sec')
                    repeat = True
            else:
                print('connection error - not on internet?')
                if URL.startswith('https:'):
                    print('Retry with http://')
                    repeat = True
                    URL = URL.replace('https:','http:')
        except requests.exceptions.Timeout as msg:
            print(f'timeout accessing {URL}')
            if GetConfigValue('debug'): print('full error=',msg)
            if timeout is not None and count < retry:
                if retry == 2:
                    timeout *= 1.25
                    print(f'retry with timout={timeout} sec')
                repeat = True
            if timeout is not None and count <= retry: repeat = True
        except requests.exceptions.ReadTimeout as msg:
            print(f'timeout reading from {URL}')
            if GetConfigValue('debug'): print('full error=',msg)
            if timeout is not None and count < retry:
                if retry == 2:
                    timeout *= 1.25
                    print(f'retry with timout={timeout} sec')
                repeat = True
        except requests.exceptions.ConnectTimeout:
            print(f'timeout accessing {URL}')
        except Exception as msg:    # other error
            print(f'Error accessing {URL}')
            if GetConfigValue('debug'): print(msg)
        finally:
            if r: r.close()
    else:
        return None

if __name__ == '__main__':
    '''What follows is called to update (or downdate) GSAS-II in a
    separate process.
    '''
    # check what type of update is being called for
    import git
    gitUpdate = False
    preupdateType = None
    updateType = None
    regressversion = None
    help = False
    project = None

    for arg in sys.argv[1:]:
        if '--git-fetch' in arg:   # pulls latest updates from server but does not apply them
            if preupdateType or updateType:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            updateType = 'fetch'
        elif '--github-tags' in arg:   # gets latest tags from github
            if preupdateType or updateType or gitUpdate:
                print(f'previous option conflicts with {arg}')
                help = True
                break
            updateType = 'tags'
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
                print(f'invalid version specified ({regressversion}) for GitHub regression')
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
        elif os.path.exists(arg):   # this is just checking
            project = arg
            pass
        else:
            print(f'unknown arg {arg}')
            help = True

    if help or len(sys.argv) == 1:
        print('''Options when running GSASIIpath.py standalone

to update/regress repository from git repository:
   python GSASIIpath.py option <project>
       where option will be one or more of the following:
            --git-fetch            downloads lastest changes from repo
                                   any other options will be ignored

            --git-stash="message"  saves local changes

            --git-reset            discards local changes

            --git-update

            --git-regress=version

            --github-tags         saves most recent tag info from GitHub

       and where <project> is an optional path reference to a .gpx file

       Note: --git-reset and --git-stash cannot be used together. Likewise
             --git-update and --git-regress cannot be used together.
             However either --git-reset or --git-stash can be used
             with either --git-update or --git-regress.

             --git-fetch cannot be used with any other options.
             --github-tags cannot be used with any other options.
''')
        sys.exit()

    if updateType == 'tags':
        # get the most recent tag numbers from the GitHub site. Do this
        # via GitHub when git access is not available. Use here
        # allows this to be done in the background.
        import requests
        url='https://github.com/AdvancedPhotonSource/GSAS-II/tags'
        releases = requests.get(url=url)
        taglist = [tag.split('"')[0] for tag in releases.text.split('AdvancedPhotonSource/GSAS-II/releases/tag/')[1:]]
        lastver = sorted([t for t in taglist if 'v' in t])[-1]
        lastnum = sorted([t for t in taglist if 'v' not in t],key=int)[-1]
        #print('tags=',lastver,lastnum)
        # add tag info to config file
        import configparser
        cfgfile = os.path.expanduser(os.path.normpath('~/.GSASII/config.ini'))
        cfg = configparser.ConfigParser()
        cfg.read(cfgfile)
        if 'version info' not in cfg:
            cfg.add_section('version info')
        cfg['version info'].update(
            {'lastVersionTag':lastver,'lastVersionNumber':lastnum})
        with open(cfgfile, 'w') as configfile:
            cfg.write(configfile)
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
            fp.write('Updates fetched\n')
        except Exception as msg:
            fp.write(f'Update failed with message {msg}\n')

        if g2repo.head.is_detached:
            fp.write('Status: reverted to an old install\n')
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
        openGitRepo(path2GSAS2).git.reset('--hard','origin/main')
        try:
            if g2repo.active_branch.name != 'main':
                g2repo.git.switch('main')
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
            g2repo.git.switch('main')
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
        # path hack for restart, when needed
        import importlib.util
        try:
            importlib.util.find_spec('GSASII.GSASIIGUI')
        except ModuleNotFoundError:
            print('Adding GSAS-II location to Python system path')
            sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))

        # now restart GSAS-II with the new version
        # G2scrpt = os.path.join(path2GSAS2,'G2.py')
        if project:
            print(f"Restart GSAS-II with project file {project!r}")
        else:
            print("Restart GSAS-II without a project file ")
        from . import GSASIIfiles
        GSASIIfiles.openInNewTerm(project)
        print ('exiting update process')
        sys.exit()
