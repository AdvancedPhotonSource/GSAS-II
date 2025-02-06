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
import importlib
try:
    pkginfo = importlib.util.find_spec('GSASII.GSASIIGUI')
except ModuleNotFoundError:
    pkginfo = None

if __name__ == '__main__':
    if pkginfo is None:  # hack path if GSASII not installed into Python
        print('GSAS-II not installed in Python; Hacking sys.path')
        sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
    from GSASII.GSASIIGUI import main
    main()
