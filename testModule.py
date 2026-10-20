# This is an example for how to run a GSAS-II module from the command
# line. This file must be in the same location as the GSASII directory
# or GSASII must be installed as a Python module.

# Here the arguments array is "faked out" so this emulates running
# `python GSASII/GSASIIstrMain.py /Users/.../2Px4G_copy2test.gpx`
# (which does not work due to import problems).

# Note that args init_globals and alter_sys might be useful to
# trigger some special actions.

import runpy
import sys
#sys.argv[1:] = ['/Users/toby/Scratch/GroupPWDR/2Px4G_copy2test.gpx']
sys.argv[1:] = ['/Users/toby/Scratch/PRseqTest/SRpr5.gpx']
runpy.run_module("GSASII.GSASIIstrMain", run_name="__main__")
