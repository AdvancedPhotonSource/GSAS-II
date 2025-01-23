'''Unit tests for code in GSASIIspc.py. Also exercises the pyspg
Fortran routine
'''
import os
import sys
import numpy as np

# TODO: not sure if this is how Tom wants to handle imports
import testinp.spctestinp as spctestinp
import testinp.sgtbxtestinp as sgtbxtestinp

import importlib  # fixup path if GSASII not installed into Python
if importlib.util.find_spec('GSASII') is None:
    print('Beware: Path hacking in progress')
    os.environ["GSASII_YOLO_PATH"] = "True"
    home = os.path.dirname(__file__)
    sys.path.append(os.path.dirname(home))
from GSASII.GSASIIspc import MoveToUnitCell, SpcGroup, SytSym

# self-test materials follow. Requires files in directory testinp
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

def test_MoveToUnitCell():
    '''self-test #0: exercise MoveToUnitCell'''
    _ReportTest()
    msg = "MoveToUnitCell failed"
    assert (MoveToUnitCell([1,2,3])[0] == [0,0,0]).all, msg
    assert (MoveToUnitCell([2,-1,-2])[0] == [0,0,0]).all, msg
    assert abs(MoveToUnitCell(np.array([-.1]))[0]-0.9)[0] < 1e-6, msg
    assert abs(MoveToUnitCell(np.array([.1]))[0]-0.1)[0] < 1e-6, msg
selftestlist.append(test_MoveToUnitCell)

def test_SpcGroup1():
    '''self-test #1: SpcGroup against previous results'''
    _ReportTest()
    testdir = os.path.join(os.path.split(os.path.abspath( __file__ ))[0],'testinp')
    if os.path.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    def CompareSpcGroup(spc, referr, refdict, reflist):
        'Compare output from GSASIIspc.SpcGroup with results from a previous run'
        msg0 = "CompareSpcGroup failed on space group %s" % spc
        result = SpcGroup(spc)
        if result[0] == referr and referr > 0: return True
        for key in refdict:
            if key == 'SGCen' or key == 'SGOps':
                assert len(refdict[key])==len(result[1][key]),f'{msg0}, {key} length'
            else:
                assert refdict[key]==result[1][key],f'{msg0}, key={key}'
        key = 'SGCen'
        indices = list(range(len(result[1][key])))
        for item in refdict[key]:
            for i,j in enumerate(indices):
                if np.allclose(result[1][key][j],item):
                    indices.pop(i)
                    break
            else:
                assert False,f'{msg0} no {key} matches center {item}'
        key = 'SGOps'
        indices = list(range(len(result[1][key])))
        for k,item in enumerate(refdict[key]):
            for i,j in enumerate(indices):
                if np.allclose(result[1][key][j][0],item[0]) and np.allclose(result[1][key][j][1],item[1]):
                    indices.pop(i)
                    break
            else:
                assert False,f'{msg0} no {key} matches sym op #{k}'
    for spc in spctestinp.SGdat:
        CompareSpcGroup(spc, 0, spctestinp.SGdat[spc], spctestinp.SGlist[spc] )
selftestlist.append(test_SpcGroup1)

def test_SpcGroup2():
    '''self-test #2: SpcGroup against cctbx (sgtbx) computations'''
    _ReportTest()
    testdir = os.path.join(os.path.split(os.path.abspath( __file__ ))[0],'testinp')
    if os.path.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    def CompareWcctbx(spcname, cctbx_in, debug=0):
        'Compare output from GSASIIspc.SpcGroup with results from cctbx.sgtbx'
        cctbx = cctbx_in[:] # make copy so we don't delete from the original
        spc = (SpcGroup(spcname))[1]
        if debug: print (spc['SpGrp'])
        if debug: print (spc['SGCen'])
        latticetype = spcname.strip().upper()[0]
        # lattice type of R implies Hexagonal centering", fix the rhombohedral settings
        if latticetype == "R" and len(spc['SGCen']) == 1: latticetype = 'P'
        assert latticetype == spc['SGLatt'], "Failed: %s does not match Lattice: %s" % (spcname, spc['SGLatt'])
        onebar = [1]
        if spc['SGInv']: onebar.append(-1)
        for (op,off) in spc['SGOps']:
            for inv in onebar:
                for cen in spc['SGCen']:
                    noff = off + cen
                    noff = MoveToUnitCell(noff)[0]
                    mult = tuple((op*inv).ravel().tolist())
                    if debug: print ("\n%s: %s + %s" % (spcname,mult,noff))
                    for refop in cctbx:
                        if debug: print (refop)
                        # check the transform
                        if refop[:9] != mult: continue
                        if debug: print ("mult match")
                        # check the translation
                        reftrans = list(refop[-3:])
                        reftrans = MoveToUnitCell(reftrans)[0]
                        if all(abs(noff - reftrans) < 1.e-5):
                            cctbx.remove(refop)
                            break
                    else:
                        assert False, "failed on %s:\n\t %s + %s" % (spcname,mult,noff)
    for key in sgtbxtestinp.sgtbx:
        CompareWcctbx(key, sgtbxtestinp.sgtbx[key])
selftestlist.append(test_SpcGroup2)

def test_SytSym():
    '''self-test #3: exercise SytSym against selected data from IT Volume A'''
    # N.B. SytSym tests GetOprPtrName, GenAtom, GetKNsym
    _ReportTest()
    def ExerciseSiteSym (spc, crdlist):
        'compare site symmetries and multiplicities for a specified space group'
        msg = f"failed on site sym test for {spc}"
        (E,S) = SpcGroup(spc)
        assert not E, msg
        for t in crdlist:
            symb, m, n, od = SytSym(t[0],S)
            if symb.strip() != t[2].strip():
                #GSASIIpath.IPyBreak_base()
                print(f'for {spc} @ {t[0]}, site sym mismatch {symb} != {t[2]} (warning)')
            assert m == t[1],f'{msg}: multiplicity @ {t[0]}'
            #assert symb.strip() == t[2].strip()

    ExerciseSiteSym('p 1',[
            ((0.13,0.22,0.31),1,'1'),
            ((0,0,0),1,'1'),
            ])
    ExerciseSiteSym('p -1',[
            ((0.13,0.22,0.31),2,'1'),
            ((0,0.5,0),1,'-1'),
            ])
    ExerciseSiteSym('C 2/c',[
            ((0.13,0.22,0.31),8,'1'),
            ((0.0,.31,0.25),4,'2(y)'),
            ((0.25,.25,0.5),4,'-1'),
            ((0,0.5,0),4,'-1'),
            ])
    ExerciseSiteSym('p 2 2 2',[
            ((0.13,0.22,0.31),4,'1'),
            ((0,0.5,.31),2,'2(z)'),
            ((0.5,.31,0.5),2,'2(y)'),
            ((.11,0,0),2,'2(x)'),
            ((0,0.5,0),1,'222'),
            ])
    ExerciseSiteSym('p 4/n',[
            ((0.13,0.22,0.31),8,'1'),
            ((0.25,0.75,.31),4,'2(z)'),
            ((0.5,0.5,0.5),4,'-1'),
            ((0,0.5,0),4,'-1'),
            ((0.25,0.25,.31),2,'4(z)'),
            ((0.25,.75,0.5),2,'-4(z)'),
            ((0.25,.75,0.0),2,'-4(z)'),
            ])
    ExerciseSiteSym('p 31 2 1',[
            ((0.13,0.22,0.31),6,'1'),
            ((0.13,0.0,0.833333333),3,'2(100)'),
            ((0.13,0.13,0.),3,'2(110)'),
            ])
    ExerciseSiteSym('R 3 c',[
            ((0.13,0.22,0.31),18,'1'),
            ((0.0,0.0,0.31),6,'3'),
            ])
    ExerciseSiteSym('R 3 c R',[
            ((0.13,0.22,0.31),6,'1'),
            ((0.31,0.31,0.31),2,'3(111)'),
            ])
    ExerciseSiteSym('P 63 m c',[
            ((0.13,0.22,0.31),12,'1'),
            ((0.11,0.22,0.31),6,'m(100)'),
            ((0.333333,0.6666667,0.31),2,'3m(100)'),
            ((0,0,0.31),2,'3m(100)'),
            ])
    ExerciseSiteSym('I a -3',[
            ((0.13,0.22,0.31),48,'1'),
            ((0.11,0,0.25),24,'2(x)'),
            ((0.11,0.11,0.11),16,'3(111)'),
            ((0,0,0),8,'-3(111)'),
            ])
selftestlist.append(test_SytSym)

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    for test in selftestlist:
        test()
    print ("OK")
