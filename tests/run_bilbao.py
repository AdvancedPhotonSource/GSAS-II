'''
run_bilbao.py
=====================

Run self-tests using the Bilbao Crystallographic Server 
checking that the server is accessible and that the GSAS-II 
routines that access the server work.
No attempt is made for any comprehensive testing of the server. 

Routines tested:

* createStdSetting (uses nph-cif2std)
* subBilbaoCheckLattice (uses nph-pseudolattice)
* GetNonStdSubgroupsmag (uses subgrmag1_general_GSAS.pl)
* GetNonStdSubgroups (uses subgrmag1_general_GSAS.pl)
* GetStdSGset (uses checkgr.pl)

The other major routines are:

* BilbaoSymSearch1 (pseudosym)
* BilbaoSymSearch2 (pseudosym)
* BilbaoReSymSearch (pseudosym)

These are not tested as they use the GUI extensively, but there is 
a test that:

* directly calls pseudosym & checks output with scanBilbaoSymSearch1
'''

import os
import sys
import importlib.util
home = os.path.dirname(__file__)
if importlib.util.find_spec('GSASII') is None: # hack path if GSASII not installed into Python
    sys.path.append(os.path.dirname(home))
from GSASII import GSASIIobj
from GSASII import SUBGROUPS
from GSASII import GSASIIpath
from GSASII import GSASIIspc as G2spc

def test_createStdSetting():
    rd = GSASIIobj.ImportPhase('Null')
    rd.Phase = {}
    rd.Phase['General'] = {}
    rd.SymOps['xyz'] = None
    cifFile = os.path.join(home,'testinp','diamond.cif')
    SUBGROUPS.createStdSetting(cifFile,rd)
    assert len(rd.Phase['Atoms']) == 1
    assert rd.Phase['General']['Cell'][1] == 3.5668
    assert len(rd.Phase['General']['SGData']['SGOps']) == 24

def test_CheckLattice():
    print('test PSEUDOLATTICE')
    cell = [14.259, 22.539, 8.741, 90.0, 114.1, 90.0]
    page = SUBGROUPS.subBilbaoCheckLattice(12,cell,5)
    found = SUBGROUPS.parseBilbaoCheckLattice(page)
    assert found[0][0][0] == found[0][0][1] == 22.5418
    assert found[0][0][5] == 120.

def test_pseudosym():
    postdict = {'formulae': '', 'cifile': '', 'filename': '', 'what': 'minsup',
                    'maxik': '1', 'onlythisik': '1', 'mysuper2': '211',
                    'x1': '1', 'x2': '0', 'x3': '0',
                    'y1': '0', 'y2': '1', 'y3': '0',
                    'z1': '0', 'z2': '0', 'z3': '1', 'x4': '0',
                    'y4': '0', 'z4': '0',
                    'angtol': '5', 'submit': 'Show', 'maxdelta': '2.0',
                    'stru': '# Space Group ITA number\n227\n# Lattice parameters\n5.43123 5.43123 5.43123 90.0 90.0 90.0\n# number of atoms & atoms\n1\nSi   1 - 0.12500 0.12500 0.12500\n'}
    savedcookies = {}
    URL = SUBGROUPS.bilbaoSite+ SUBGROUPS.pseudosym
    print('test PSEUDOSYM')
    page0 = GSASIIpath.postURL(URL,postdict,
                             getcookie=savedcookies,timeout=SUBGROUPS.timeout)
    res = SUBGROUPS.scanBilbaoSymSearch1(page0,postdict)+[savedcookies]
    assert res[2]['1'][0] == 'Pn-3m'
    assert abs(float(res[2]['1'][5].split()[0]) -  5.43123/2) < 2e-5

def test_SUBGROUPSMAG():
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    
    print('test SUBGROUPSMAG')
    results,baseList = SUBGROUPS.GetNonStdSubgroupsmag(
        SGData,('0','0','0',' ',' ',' ',' ',' ',' ',' '))
    assert len(results) == 322
    assert results[0][0] == "Fd'-3'm'"
    assert results[1][0] == "Fd-3m'"
    #print(results)

def test_SUBGROUPS():
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    print('\n\ntest SUBGROUPS')
    results,baseList = SUBGROUPS.GetNonStdSubgroups(
        SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    assert len(results) == 51
    assert results[0][0] == "P21/c"
    assert results[-1][0] == "P1"
    
def test_GetStdSGset():
    print('\n\ntest Bilbao IDENTIFY GROUP')
    sgnum,sgsym,xmat,xoff = SUBGROUPS.GetStdSGset(G2spc.SpcGroup('R 3 C r')[1])
    assert sgnum == 161
    assert sgsym == 'R 3 c'
    assert xoff == [0,0,0]

#GSASIIpath.IPyBreak_base()

    
if __name__ == '__main__':
    # run self-tests
    test_GetStdSGset()
    test_SUBGROUPSMAG()
    test_SUBGROUPS()
    test_pseudosym()
    test_CheckLattice()
    test_createStdSetting()
    print ("OK")
