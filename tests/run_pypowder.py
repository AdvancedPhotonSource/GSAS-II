# compute a single peak with pypowder
import os
import sys
import importlib.util
import numpy as np

home = os.path.dirname(__file__)
G2loc = None
try: 
    G2loc = importlib.util.find_spec('GSASII.GSASIIscriptable')
except ModuleNotFoundError:
    print('ModuleNotFound for GSASII.GSASIIscriptable')

if G2loc is None: # fixup path if GSASII not installed into Python
    print('GSAS-II not installed in Python; Hacking sys.path')
    sys.path.append(os.path.dirname(home))

from GSASII import GSASIIpath
GSASIIpath.SetBinaryPath()

if GSASIIpath.binaryPath:
    import  pypowder as pyd
else:
    from . import pypowder as pyd

# these values generate a nearly symmetric peak
xdata = np.arange(18,22,.05)
pos = 20.
sig = 525.9994695543994
gam = 2.1444606947716025
shl = 0.002
ydata = pyd.pypsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)

# these values generate a peak with significant low-angle broadening
#xdata = np.arange(0.1,5,.05)
#pos = 2
#sig = 525.9994695543994
#gam = 2.1444606947716025
#shl = 0.035
#ydata = pyd.pypsvfcjo(len(xdata),xdata-pos,pos,sig,gam,shl)

import matplotlib.pyplot as plt
plt.plot(xdata,ydata)
plt.show()

