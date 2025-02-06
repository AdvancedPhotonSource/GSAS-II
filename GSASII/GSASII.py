# Starts GSAS-II from clone of GitHub repo when contents is not
# installed into current Python.
#
# This can be used when GSAS-II is installed, but why bother as the command:
#    python -c "from GSASII.GSASIIGUI import main; main()"
# will do just as well. 
import os
import sys
import importlib
try:
    pkginfo = importlib.util.find_spec('GSASII.GSASIIGUI')
except ModuleNotFoundError:
    pkginfo = None

if __name__ == '__main__':
    if pkginfo is None:  # fixup path if GSASII not installed into Python
        os.environ["GSASII_YOLO_PATH"] = "True"
        sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
    from GSASII.GSASIIGUI import main
    if pkginfo is None:  # fixup path if GSASII not installed into Python
        print('GSAS-II not installed in Python: Hacking sys.path')
    main()
