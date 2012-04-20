# determine a binary path for the pyd files based on the host OS and the python version,  
# path is relative to location of the script that is called as well as this file
# this must be imported before anything that imports any .pyd/.so file for GSASII
import os.path as ospath
import sys
import platform
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
for loc in sys.path[0],ospath.split(__file__)[0]:
    if bindir:
        if ospath.exists(ospath.join(loc,bindir)) and ospath.join(loc,bindir) not in sys.path: 
            sys.path.insert(0,ospath.join(loc,bindir))
        # is there a bin directory? (created by a local compile), if so put
        # that at the top of the path
    if ospath.exists(ospath.join(loc,'bin')):
        bindir = 'bin'
        if ospath.join(loc,'bin') not in sys.path: 
            sys.path.insert(0,ospath.join(loc,bindir))
if bindir == None:
    print "Warning GSAS-II binary libraries not found, some sections of code will not function"
