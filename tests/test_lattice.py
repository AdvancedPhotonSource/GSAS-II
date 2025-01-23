'''Unit tests for code in GSASIIlattice.py.
'''

import os
import sys
import numpy as np

import importlib  # fixup path if GSASII not installed into Python
if importlib.util.find_spec('GSASII') is None:
    print('Beware: Path hacking in progress')
    os.environ["GSASII_YOLO_PATH"] = "True"
    home = os.path.dirname(__file__)
    sys.path.append(os.path.dirname(home))
import GSASII.GSASIIspc as G2spc
from GSASII.GSASIIlattice import *

import testinp.sgtbxlattinp as sgtbxlattinp

# self-test materials follow.
selftestlist = []
'''Defines a list of self-tests'''
selftestquiet = True
def _ReportTest():
    'Report name and doc string of current routine when ``selftestquiet`` is False'
    if not selftestquiet:
        import inspect
        caller = inspect.stack()[1][3]
        doc = eval(caller).__doc__
        if doc is not None:
            print(f'testing {os.path.split(__file__)[1]} with {caller}\n\t({doc})')
        else:
            print(f'testing {os.path.split(__file__)[1]} with {caller}')

NeedTestData = True
def TestData():
    array = np.array
    global NeedTestData
    NeedTestData = False
    global CellTestData
    # output from uctbx computed on platform darwin on 2010-05-28
    CellTestData = [
# cell, g, G, cell*, V, V*
  [(4, 4, 4, 90, 90, 90),
   array([[  1.60000000e+01,   9.79717439e-16,   9.79717439e-16],
       [  9.79717439e-16,   1.60000000e+01,   9.79717439e-16],
       [  9.79717439e-16,   9.79717439e-16,   1.60000000e+01]]), array([[  6.25000000e-02,   3.82702125e-18,   3.82702125e-18],
       [  3.82702125e-18,   6.25000000e-02,   3.82702125e-18],
       [  3.82702125e-18,   3.82702125e-18,   6.25000000e-02]]), (0.25, 0.25, 0.25, 90.0, 90.0, 90.0), 64.0, 0.015625],
# cell, g, G, cell*, V, V*
  [(4.0999999999999996, 5.2000000000000002, 6.2999999999999998, 100, 80, 130),
   array([[ 16.81      , -13.70423184,   4.48533243],
       [-13.70423184,  27.04      ,  -5.6887143 ],
       [  4.48533243,  -5.6887143 ,  39.69      ]]), array([[ 0.10206349,  0.05083339, -0.00424823],
       [ 0.05083339,  0.06344997,  0.00334956],
       [-0.00424823,  0.00334956,  0.02615544]]), (0.31947376387537696, 0.25189277536327803, 0.16172643497798223, 85.283666420376008, 94.716333579624006, 50.825714168082683), 100.98576357983838, 0.0099023858863968445],
# cell, g, G, cell*, V, V*
  [(3.5, 3.5, 6, 90, 90, 120),
   array([[  1.22500000e+01,  -6.12500000e+00,   1.28587914e-15],
       [ -6.12500000e+00,   1.22500000e+01,   1.28587914e-15],
       [  1.28587914e-15,   1.28587914e-15,   3.60000000e+01]]), array([[  1.08843537e-01,   5.44217687e-02,   3.36690552e-18],
       [  5.44217687e-02,   1.08843537e-01,   3.36690552e-18],
       [  3.36690552e-18,   3.36690552e-18,   2.77777778e-02]]), (0.32991443953692895, 0.32991443953692895, 0.16666666666666669, 90.0, 90.0, 60.000000000000021), 63.652867178156257, 0.015710211406520427],
  ]
    global CoordTestData
    CoordTestData = [
# cell, ((frac, ortho),...)
  ((4,4,4,90,90,90,), [
 ((0.10000000000000001, 0.0, 0.0),(0.40000000000000002, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(2.4492935982947065e-17, 0.40000000000000002, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(2.4492935982947065e-17, -2.4492935982947065e-17, 0.40000000000000002)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(0.40000000000000013, 0.79999999999999993, 1.2)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(0.80000000000000016, 1.2, 0.40000000000000002)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(1.2, 0.80000000000000004, 0.40000000000000002)),
 ((0.5, 0.5, 0.5),(2.0, 1.9999999999999998, 2.0)),
]),
# cell, ((frac, ortho),...)
  ((4.1,5.2,6.3,100,80,130,), [
 ((0.10000000000000001, 0.0, 0.0),(0.40999999999999998, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(-0.33424955703700043, 0.39834311042186865, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(0.10939835193016617, -0.051013289294572106, 0.6183281045774256)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(0.069695941716497567, 0.64364635296002093, 1.8549843137322766)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(-0.073350319180835066, 1.1440160419710339, 0.6183281045774256)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(0.67089923785616512, 0.74567293154916525, 0.6183281045774256)),
 ((0.5, 0.5, 0.5),(0.92574397446582857, 1.7366491056364828, 3.0916405228871278)),
]),
# cell, ((frac, ortho),...)
  ((3.5,3.5,6,90,90,120,), [
 ((0.10000000000000001, 0.0, 0.0),(0.35000000000000003, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(-0.17499999999999993, 0.3031088913245536, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(3.6739403974420595e-17, -3.6739403974420595e-17, 0.60000000000000009)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(2.7675166561703527e-16, 0.60621778264910708, 1.7999999999999998)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(0.17500000000000041, 0.90932667397366063, 0.60000000000000009)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(0.70000000000000018, 0.6062177826491072, 0.60000000000000009)),
 ((0.5, 0.5, 0.5),(0.87500000000000067, 1.5155444566227676, 3.0)),
]),
]
    global LaueTestData             #generated by GSAS
    LaueTestData = {
    'R 3 m':[(4.,4.,6.,90.,90.,120.),((1,0,1,6),(1,0,-2,6),(0,0,3,2),(1,1,0,6),(2,0,-1,6),(2,0,2,6),
        (1,1,3,12),(1,0,4,6),(2,1,1,12),(2,1,-2,12),(3,0,0,6),(1,0,-5,6),(2,0,-4,6),(3,0,-3,6),(3,0,3,6),
        (0,0,6,2),(2,2,0,6),(2,1,4,12),(2,0,5,6),(3,1,-1,12),(3,1,2,12),(1,1,6,12),(2,2,3,12),(2,1,-5,12))],
    'R 3':[(4.,4.,6.,90.,90.,120.),((1,0,1,6),(1,0,-2,6),(0,0,3,2),(1,1,0,6),(2,0,-1,6),(2,0,2,6),(1,1,3,6),
        (1,1,-3,6),(1,0,4,6),(3,-1,1,6),(2,1,1,6),(3,-1,-2,6),(2,1,-2,6),(3,0,0,6),(1,0,-5,6),(2,0,-4,6),
        (2,2,0,6),(3,0,3,6),(3,0,-3,6),(0,0,6,2),(3,-1,4,6),(2,0,5,6),(2,1,4,6),(4,-1,-1,6),(3,1,-1,6),
        (3,1,2,6),(4,-1,2,6),(2,2,-3,6),(1,1,-6,6),(1,1,6,6),(2,2,3,6),(2,1,-5,6),(3,-1,-5,6))],
    'P 3':[(4.,4.,6.,90.,90.,120.),((0,0,1,2),(1,0,0,6),(1,0,1,6),(0,0,2,2),(1,0,-1,6),(1,0,2,6),(1,0,-2,6),
        (1,1,0,6),(0,0,3,2),(1,1,1,6),(1,1,-1,6),(1,0,3,6),(1,0,-3,6),(2,0,0,6),(2,0,-1,6),(1,1,-2,6),
        (1,1,2,6),(2,0,1,6),(2,0,-2,6),(2,0,2,6),(0,0,4,2),(1,1,-3,6),(1,1,3,6),(1,0,-4,6),(1,0,4,6),
        (2,0,-3,6),(2,1,0,6),(2,0,3,6),(3,-1,0,6),(2,1,1,6),(3,-1,-1,6),(2,1,-1,6),(3,-1,1,6),(1,1,4,6),
        (3,-1,2,6),(3,-1,-2,6),(1,1,-4,6),(0,0,5,2),(2,1,2,6),(2,1,-2,6),(3,0,0,6),(3,0,1,6),(2,0,4,6),
        (2,0,-4,6),(3,0,-1,6),(1,0,-5,6),(1,0,5,6),(3,-1,-3,6),(2,1,-3,6),(2,1,3,6),(3,-1,3,6),(3,0,-2,6),
        (3,0,2,6),(1,1,5,6),(1,1,-5,6),(2,2,0,6),(3,0,3,6),(3,0,-3,6),(0,0,6,2),(2,0,-5,6),(2,1,-4,6),
        (2,2,-1,6),(3,-1,-4,6),(2,2,1,6),(3,-1,4,6),(2,1,4,6),(2,0,5,6),(1,0,-6,6),(1,0,6,6),(4,-1,0,6),
        (3,1,0,6),(3,1,-1,6),(3,1,1,6),(4,-1,-1,6),(2,2,2,6),(4,-1,1,6),(2,2,-2,6),(3,1,2,6),(3,1,-2,6),
        (3,0,4,6),(3,0,-4,6),(4,-1,-2,6),(4,-1,2,6),(2,2,-3,6),(1,1,6,6),(1,1,-6,6),(2,2,3,6),(3,-1,5,6),
        (2,1,5,6),(2,1,-5,6),(3,-1,-5,6))],
    'P 3 m 1':[(4.,4.,6.,90.,90.,120.),((0,0,1,2),(1,0,0,6),(1,0,-1,6),(1,0,1,6),(0,0,2,2),(1,0,-2,6),
        (1,0,2,6),(1,1,0,6),(0,0,3,2),(1,1,1,12),(1,0,-3,6),(1,0,3,6),(2,0,0,6),(1,1,2,12),(2,0,1,6),
        (2,0,-1,6),(0,0,4,2),(2,0,-2,6),(2,0,2,6),(1,1,3,12),(1,0,-4,6),(1,0,4,6),(2,0,3,6),(2,1,0,12),
        (2,0,-3,6),(2,1,1,12),(2,1,-1,12),(1,1,4,12),(2,1,2,12),(0,0,5,2),(2,1,-2,12),(3,0,0,6),(1,0,-5,6),
        (3,0,1,6),(3,0,-1,6),(1,0,5,6),(2,0,4,6),(2,0,-4,6),(2,1,3,12),(2,1,-3,12),(3,0,-2,6),(3,0,2,6),
        (1,1,5,12),(3,0,-3,6),(0,0,6,2),(2,2,0,6),(3,0,3,6),(2,1,4,12),(2,2,1,12),(2,0,5,6),(2,1,-4,12),
        (2,0,-5,6),(1,0,-6,6),(1,0,6,6),(3,1,0,12),(3,1,-1,12),(3,1,1,12),(2,2,2,12),(3,1,2,12),
        (3,0,4,6),(3,1,-2,12),(3,0,-4,6),(1,1,6,12),(2,2,3,12))],
    'P 3 1 m':[(4.,4.,6.,90.,90.,120.),((0,0,1,2),(1,0,0,6),(0,0,2,2),(1,0,1,12),(1,0,2,12),(1,1,0,6),
        (0,0,3,2),(1,1,-1,6),(1,1,1,6),(1,0,3,12),(2,0,0,6),(2,0,1,12),(1,1,2,6),(1,1,-2,6),(2,0,2,12),
        (0,0,4,2),(1,1,-3,6),(1,1,3,6),(1,0,4,12),(2,1,0,12),(2,0,3,12),(2,1,1,12),(2,1,-1,12),(1,1,-4,6),
        (1,1,4,6),(0,0,5,2),(2,1,-2,12),(2,1,2,12),(3,0,0,6),(1,0,5,12),(2,0,4,12),(3,0,1,12),(2,1,-3,12),
        (2,1,3,12),(3,0,2,12),(1,1,5,6),(1,1,-5,6),(3,0,3,12),(0,0,6,2),(2,2,0,6),(2,1,-4,12),(2,0,5,12),
        (2,2,-1,6),(2,2,1,6),(2,1,4,12),(3,1,0,12),(1,0,6,12),(2,2,2,6),(3,1,-1,12),(2,2,-2,6),(3,1,1,12),
        (3,1,-2,12),(3,0,4,12),(3,1,2,12),(1,1,-6,6),(2,2,3,6),(2,2,-3,6),(1,1,6,6))],
    }

    global FLnhTestData
    FLnhTestData = [{
    'C(4,0,0)': (0.965, 0.42760447),
    'C(2,0,0)': (1.0122, -0.80233610),
    'C(2,0,2)': (0.0061, 8.37491546E-03),
    'C(6,0,4)': (-0.0898, 4.37985696E-02),
    'C(6,0,6)': (-0.1369, -9.04081762E-02),
    'C(6,0,0)': (0.5935, -0.18234928),
    'C(4,0,4)': (0.1872, 0.16358127),
    'C(6,0,2)': (0.6193, 0.27573633),
    'C(4,0,2)': (-0.1897, 0.12530720)},[1,0,0]]

def test_gmat():
    if NeedTestData: TestData()
    msg = 'test cell2Gmat, fillgmat, Gmat2cell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = cell2Gmat(cell)
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg
        tcell = Gmat2cell(g)
        assert np.allclose(cell,tcell),msg
        tcell = Gmat2cell(G)
        assert np.allclose(tcell,trcell),msg
selftestlist.append(test_gmat)

def test_Avec():
    'test cell2A and A2Gmat'
    _ReportTest()
    if NeedTestData: TestData()
    msg = 'test cell2A and A2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = A2Gmat(cell2A(cell))
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg
selftestlist.append(test_Avec)

def test_2cell():
    'test Gmat2A, A2cell, A2Gmat, Gmat2cell'
    _ReportTest()
    if NeedTestData: TestData()
    msg = 'test Gmat2A, A2cell, A2Gmat, Gmat2cell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = cell2Gmat(cell)
        tcell = A2cell(Gmat2A(G))
        assert np.allclose(cell,tcell),msg
selftestlist.append(test_2cell)

def test_invcell():
    'test invcell2Gmat'
    _ReportTest()
    if NeedTestData: TestData()
    msg = 'test invcell2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = invcell2Gmat(trcell)
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg
    msg = 'test A2invcell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        rcell = A2invcell(cell2A(cell))
        assert np.allclose(rcell,trcell),msg
selftestlist.append(test_invcell)

def test_V():
    'test calc_rVsq, calc_rV, calc_V'
    _ReportTest()
    if NeedTestData: TestData()
    msg = 'test calc_rVsq, calc_rV, calc_V'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        assert np.allclose(calc_rV(cell2A(cell)),trV), msg
        assert np.allclose(calc_V(cell2A(cell)),tV), msg
selftestlist.append(test_V)

def test_2AB():
    'test cell2AB'
    _ReportTest()
    if NeedTestData: TestData()
    msg = 'test cell2AB'
    for (cell,coordlist) in CoordTestData:
        A,B = cell2AB(cell)
        for (frac,ortho) in coordlist:
            to = np.inner(A,frac)
            tf = np.inner(B,to)
            assert np.allclose(ortho,to), msg
            assert np.allclose(frac,tf), msg
            to = np.sum(A*frac,axis=1)
            tf = np.sum(B*to,axis=1)
            assert np.allclose(ortho,to), msg
            assert np.allclose(frac,tf), msg
selftestlist.append(test_2AB)

def test_Brav():
    'test GetBraviasNum() and GenHBravais()'
    _ReportTest()
    testdir = os.path.join(os.path.split(os.path.abspath( __file__ ))[0],'testinp')
    if os.path.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    derror = 1e-4
    def indexmatch(hklin, hkllist, system):
        for hklref in hkllist:
            hklref = list(hklref)
            # these permutations are far from complete, but are sufficient to
            # allow the test to complete
            if system == 'cubic':
                permlist = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1),]
            elif system == 'monoclinic':
                permlist = [(1,2,3),(-1,2,-3)]
            else:
                permlist = [(1,2,3)]

            for perm in permlist:
                hkl = [abs(i) * hklin[abs(i)-1] / i for i in perm]
                if hkl == hklref: return True
                if [-i for i in hkl] == hklref: return True
        else:
            return False

    for key in sgtbxlattinp.sgtbx7:
        spdict = G2spc.SpcGroup(key)
        cell = sgtbxlattinp.sgtbx7[key][0]
        system = spdict[1]['SGSys']
        center = spdict[1]['SGLatt']

        bravcode = GetBraviasNum(center, system)

        g2list = GenHBravais(sgtbxlattinp.dmin, bravcode, cell2A(cell))

        assert len(sgtbxlattinp.sgtbx7[key][1]) == len(g2list), 'Reflection lists differ for %s' % key
        for h,k,l,d,num in g2list:
            for hkllist,dref in sgtbxlattinp.sgtbx7[key][1]:
                if abs(d-dref) < derror:
                    if indexmatch((h,k,l,), hkllist, system):
                        break
            else:
                assert 0,'No match for %s at %s (%s)' % ((h,k,l),d,key)
selftestlist.append(test_Brav)

def test_Laue():
    'test GenHLaue'
    _ReportTest()
    derror = 1e-4
    dmin = sgtbxlattinp.dmin

    def indexmatch(hklin, hklref, system, axis):
        # these permutations are far from complete, but are sufficient to
        # allow the test to complete
        if system == 'cubic':
            permlist = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1),]
        elif system == 'monoclinic' and axis=='b':
            permlist = [(1,2,3),(-1,2,-3)]
        elif system == 'monoclinic' and axis=='a':
            permlist = [(1,2,3),(1,-2,-3)]
        elif system == 'monoclinic' and axis=='c':
            permlist = [(1,2,3),(-1,-2,3)]
        elif system == 'trigonal':
            permlist = [(1,2,3),(2,1,3),(-1,-2,3),(-2,-1,3)]
        elif system == 'rhombohedral':
            permlist = [(1,2,3),(2,3,1),(3,1,2)]
        else:
            permlist = [(1,2,3)]

        hklref = list(hklref)
        for perm in permlist:
            hkl = [abs(i) * hklin[abs(i)-1] / i for i in perm]
            if hkl == hklref: return True
            if [-i for i in hkl] == hklref: return True
        return False

    for key in sgtbxlattinp.sgtbx8:
        spdict = G2spc.SpcGroup(key)[1]
        cell = sgtbxlattinp.sgtbx8[key][0]
        Axis = spdict['SGUniq']
        system = spdict['SGSys']

        g2list = GenHLaue(dmin,spdict,cell2A(cell))
        #if len(g2list) != len(sgtbxlattinp.sgtbx8[key][1]):
        #    print 'failed',key,':' ,len(g2list),'vs',len(sgtbxlattinp.sgtbx8[key][1])
        #    print 'GSAS-II:'
        #    for h,k,l,d in g2list: print '  ',(h,k,l),d
        #    print 'SGTBX:'
        #    for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]: print '  ',hkllist,dref
        assert len(g2list) == len(sgtbxlattinp.sgtbx8[key][1]), (
            'Reflection lists differ for %s' % key
            )
        #match = True
        for h,k,l,d in g2list:
            for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]:
                if abs(d-dref) < derror:
                    if indexmatch((h,k,l,), hkllist, system, Axis): break
            else:
                assert 0,'No match for %s at %s (%s)' % ((h,k,l),d,key)
                #match = False
        #if not match:
            #for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]: print '  ',hkllist,dref
            #print center, Laue, Axis, system

    for spc in LaueTestData:
        data = LaueTestData[spc]
        cell = data[0]
        hklm = np.array(data[1])
        H = hklm[-1][:3]
        hklO = hklm.T[:3].T
        A = cell2A(cell)
        dmin = 1./np.sqrt(calc_rDsq(H,A))
        SGData = G2spc.SpcGroup(spc)[1]
        hkls = np.array(GenHLaue(dmin,SGData,A))
        hklN = hkls.T[:3].T
        #print spc,hklO.shape,hklN.shape
        err = True
        for H in hklO:
            if H not in hklN:
                print ('%d %s'%(H,' missing from hkl from GSASII'))
                err = False
        assert(err)
selftestlist.append(test_Laue)

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    for test in selftestlist:
        test()
    print ("OK")
