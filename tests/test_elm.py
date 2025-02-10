import os
import sys
import importlib.util
if importlib.util.find_spec('GSASII') is None: # hack path if GSASII not installed into Python
    home = os.path.dirname(__file__)
    sys.path.append(os.path.dirname(home))

from GSASII import GSASIIElem

def test_get_xsection():
    xsection = GSASIIElem.GetXsectionCoeff('Fe')
    assert len(xsection) > 0

if __name__ == '__main__':
    # run self-tests
    test_get_xsection()
    print ("OK")
