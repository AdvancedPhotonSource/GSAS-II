# Starts GSAS-II from clone of GitHub repo when contents is not
# installed into current Python.
#
# This can be used when GSAS-II is installed, but why bother as the command:
#    python -c "from GSASII.GSASIIGUI import main; main()"
# will do just as well.
# Even better:
#    GSASII_NOPATHHACKING="true" python -c "from GSASII.GSASIIGUI import main; main()"

import os
import sys
# patch 4/25: cleanup old GSASII.py script if still around (replaced by this file)
oldg2script = os.path.join(os.path.dirname(__file__),'GSASII.py')
if os.path.exists(oldg2script):
    os.remove(oldg2script)
    print(f'Cleanup: removing old {oldg2script!r} file')
# end patch
import importlib.util
try:
    pkginfo = importlib.util.find_spec('GSASII.GSASIIGUI')
except ModuleNotFoundError:
    pkginfo = None

if __name__ == '__main__':
    if pkginfo is None:  # hack path if GSASII not installed into Python
        print('Adding GSAS-II location to Python system path')
        sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
    from GSASII.GSASIIGUI import main
    main()
