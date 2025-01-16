########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''At present,
only modules :mod:`GSASIIspc` and :mod:`GSASIIlattice` have self-tests
and these have not been tested or updated in many, many years. 
'''

from . import GSASIIspc
from . import GSASIIlattice
def test_GSASIIspc():
    '''Test registered self-tests in ``GSASIIspc``.
    Takes no input and returns nothing. Throws an Exception if a test fails.
    '''
    #GSASIIspc.selftestquiet = False
    for test in GSASIIspc.selftestlist:
        test()
def test_GSASIIlattice():
    '''Test registered self-tests in ``GSASIIlattice``.
    Takes no input and returns nothing. Throws an Exception if a test fails.
    '''
    #GSASIIlattice.selftestquiet = False
    for test in GSASIIlattice.selftestlist:
        test()

if __name__ == '__main__':
    test_GSASIIspc()
    test_GSASIIlattice()
    print("OK")
