# determine a binary path for the pyd files based on the host OS and the python version,  
# path is relative to location of this file
import os.path as ospath
import sys
if sys.platform == "win32":
    bindir = 'binwin%d.%d' % sys.version_info[0:2]
elif sys.platform == "darwin":
    bindir = 'binmac%d.%d' % sys.version_info[0:2]
else:
    bindir = 'bin'
if ospath.exists(ospath.join(sys.path[0],bindir)) and ospath.join(sys.path[0],bindir) not in sys.path: 
    sys.path.insert(0,ospath.join(sys.path[0],bindir))
