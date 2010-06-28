# determine a binary path for the pyd files based on the host OS and the python version,  
# path is relative to location of this file
#this must be imported before anything that imports any pyd file for GSASII
import os.path as ospath
import sys
bindir = None
if sys.platform == "win32":
    bindir = 'binwin%d.%d' % sys.version_info[0:2]
elif sys.platform == "darwin":
    bindir = 'binmac%d.%d' % sys.version_info[0:2]
elif sys.platform == "linux2":
    bindir = 'binlinux%d.%d' % sys.version_info[0:2]
if bindir:
    if ospath.exists(ospath.join(sys.path[0],bindir)) and ospath.join(sys.path[0],bindir) not in sys.path: 
        sys.path.insert(0,ospath.join(sys.path[0],bindir))
# is there a bin directory? (created by a local compile), if so put
# that at the top of the path
if ospath.exists(ospath.join(sys.path[0],'bin')):
    bindir = 'bin'
    if ospath.join(sys.path[0],'bin') not in sys.path: 
        sys.path.insert(0,ospath.join(sys.path[0],bindir))
if bindir == None:
    print "Warning GSAS-II binary libraries not found, some sections of code will not function"
