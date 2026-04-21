'''
run_bilbao.py
=====================

Run self-tests using the Bilbao Crystallographic Server 
checking that the server is accessible and that the GSAS-II 
routines that access the server work.

Note that this needs an API Key, which can either be set as an
environment variable (BCS_API_KEY) or as the value for the 
preverences setting also named BCS_API_KEY.

SUBGROUPS Routines tested:

* createStdSetting (uses nph-cif2std)
* GetNonStdSubgroupsmag (uses subgrmag1_general_GSAS.pl)
* GetNonStdSubgroups (uses subgrmag1_general_GSAS.pl)
* GetStdSGset (uses checkgr.pl)
* BilbaoSymSearch1 (pseudosym)
* BilbaoGetStructure (pseudosym)
* BilbaoReSymSearch (pseudosym)
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
    if page0 is None: assert False,"Web access failed"
    res = SUBGROUPS.scanBilbaoSymSearch1(page0,postdict)+[savedcookies]
    assert res[2]['1'][0] == 'Pn-3m'
    assert abs(float(res[2]['1'][5].split()[0]) -  5.43123/2) < 2e-5

def test_BilbaoSymSearch1():
    print('\n\ntest Bilbao PSEUDO')
    # perform a minsub search in p4mm
    sgnum=99
    pagelist = {}
    phase = {'General':{}}
    phase['General']['Cell'] = [False, 3.999, 3.999, 4.02, 90.0, 90.0, 90.0, 64.2878]
    phase['General']['AtomPtrs'] = [3, 1, 7, 9]
    phase['Atoms'] = [
        ['Ba1', 'Ba+2', '', 0.0, 0.0, 0.0, 1.0, '4mm(z)', 1, 'I', 0.01, 0, 0, 0, 0, 0, 0, 5790755546499520608],
        ['Ti2', 'Ti+4', '', 0.5, 0.5, 0.42, 1.0, '4mm(z)', 1, 'I', 0.01, 0, 0, 0, 0, 0, 0, 8806409101786770603],
        ['O3', 'O-2', '', 0.5, 0.5, 0.03, 1.0, '4mm(z)', 1, 'I', 0.01, 0, 0, 0, 0, 0, 0, 8862174245482840967],
        ['O4', 'O-2', '', 0.5, 0.0, 0.58, 1.0, 'mm2(z)', 2, 'I', 0.01, 0, 0, 0, 0, 0, 0, 1200819297074464544]]

    formDict,csdict,rowdict,stru = SUBGROUPS.BilbaoSymSearch1(sgnum, phase,
                                                            pagelist=pagelist)
    # finds a higher symmetry space group (P4/mmm)
    assert formDict['cs'] == ['9'] and csdict['9'] and sum(csdict.values()) == 1
    assert rowdict['9'][1] == '123'
    # get the P4/mmm structure
    structures,_ = SUBGROUPS.BilbaoGetStructure(
        csdict,rowdict,stru,pagelist=pagelist)
    assert structures['9'][0] == '123' # space group 
    assert structures['9'][2] == '4' # still 4 atoms
    assert all([   # coordinates all 0.0 or 0.5
        all([eval(i) in [0.0,0.5] for i in atom.split()[3:]])
                                  for atom in structures['9'][3:7]])
    # find structures that should be searched again
    ReSearch = SUBGROUPS.find2SearchAgain(pagelist,'')
    key = '9'
    formDict,csdict,rowdict = SUBGROUPS.BilbaoReSymSearch(
                            key+'R',ReSearch[key],pagelist=pagelist)
    # finds a cubic space grooup, Pm-3m
    assert rowdict['9'][1] == '221'
    stru = ReSearch[key]['stru']
    # get the P4/mmm structure
    structures,_ = SUBGROUPS.BilbaoGetStructure(
        csdict,rowdict,stru,pagelist=pagelist)
    assert structures['9'][0] == '221' # space group 
    assert structures['9'][2] == '3' # now 3 atoms
    assert all([
        all([eval(i) in [0.0,0.5] for i in atom.split()[3:]])
                                  for atom in structures['9'][3:7]])
    
def test_SUBGROUPSMAG():
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    print('test SUBGROUPSMAG')
    results,baseList = SUBGROUPS.GetNonStdSubgroupsmag(
        SGData,('0','0','0',' ',' ',' ',' ',' ',' ',' '))
    assert results is not None,'web access failed'
    assert len(results) == 322
    assert results[0][0] == "Fd'-3'm'"
    assert results[1][0] == "Fd-3m'"
    #print(results)

def test_SUBGROUPS():
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    print('\n\ntest Bilbao SUBGROUPS')
    results,baseList = SUBGROUPS.GetNonStdSubgroups(
        SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    assert results is not None,'web access failed'
    assert len(results) == 51
    assert results[0][0] == "P21/c"
    assert results[-1][0] == "P1"


def test_GetStdSGset():
    print('\n\ntest Bilbao IDENTIFY GROUP')
    sgnum,sgsym,xmat,xoff = SUBGROUPS.GetStdSGset(G2spc.SpcGroup('R 3 C r')[1])
    assert sgnum is not None,'web access failed'
    assert sgnum == 161
    assert sgsym == 'R 3 c'
    assert xoff == [0,0,0]

def test_createStdSetting():
    print('\n\ntest Bilbao BCS_CIF2STANDARD')
    rd = GSASIIobj.ImportPhase('Null')
    rd.Phase = {}
    rd.Phase['General'] = {}
    rd.SymOps['xyz'] = None
    cifFile = os.path.join(home,'testinp','diamond.cif')
    SUBGROUPS.createStdSetting(cifFile,rd)
    assert len(rd.Phase['Atoms']) == 1
    assert rd.Phase['General']['Cell'][1] == 3.5668
    assert len(rd.Phase['General']['SGData']['SGOps']) == 24
    
# def test_CheckLattice():
#     print('test PSEUDOLATTICE')
#     cell = [14.259, 22.539, 8.741, 90.0, 114.1, 90.0]
#     page = SUBGROUPS.subBilbaoCheckLattice(12,cell,5)
#     found = SUBGROUPS.parseBilbaoCheckLattice(page)
#     assert found[0][0][0] == found[0][0][1] == 22.5418
#     assert found[0][0][5] == 120.
    
if __name__ == '__main__':
    #GSASIIpath.InvokeDebugOpts()
    # run self-tests
    test_createStdSetting() # nph-cif2std
    test_GetStdSGset()      # checkgr_gsas.pl
    # test_pseudosym()        # nph-pseudosym, no longer available
    test_BilbaoSymSearch1()
    test_SUBGROUPSMAG()     # subgrmag1_general_GSAS.pl 
    test_SUBGROUPS()        # subgrmag1_general_GSAS.pl 
    print ("OK")
