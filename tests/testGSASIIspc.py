# Unit tests for space group codes
import sys
import os
sys.path.append('..')

import numpy as np
import GSASIIspc

def CompareSpcGroup(spc, referr, refdict, reflist): 
    'Compare output from GSASIIspc.SpcGroup with results from a previous run'
    # if an error is reported, the dictionary can be ignored
    msg = "failed on space group %s" % spc
    result = GSASIIspc.SpcGroup(spc)
    if result[0] == referr and referr > 0: return True
    keys = result[1].keys()
    #print result[1]['SpGrp']
    assert len(keys) == len(refdict.keys()), msg
    for key in keys:
        #print key, type(refdict[key])
        if key == 'SGOps':
            assert len(refdict[key]) == len(result[1][key]), msg
            for i in range(len(refdict[key])):
                assert np.allclose(result[1][key][i][0],refdict[key][i][0]), msg
                assert np.allclose(result[1][key][i][1],refdict[key][i][1]), msg
        else:
            assert result[1][key] == refdict[key], msg
    assert reflist == GSASIIspc.SGPrint(result[1]), 'SGPrint ' +msg

def CompareWcctbx(spcname, cctbx_in, debug=0):
    'Compare output from GSASIIspc.SpcGroup with results from cctbx.sgtbx'
    cctbx = cctbx_in[:] # make copy so we don't delete from the original
    spc = (GSASIIspc.SpcGroup(spcname))[1]
    if debug: print spc['SpGrp']
    if debug: print spc['SGCen']
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
                GSASIIspc.MoveToUnitCell(noff)
                mult = tuple((op*inv).ravel().tolist())
                if debug: print "\n%s: %s + %s" % (spcname,mult,noff)
                for refop in cctbx:
                    if debug: print refop
                    # check the transform
                    if refop[:9] != mult: continue
                    if debug: print "mult match"
                    # check the translation
                    reftrans = list(refop[-3:])
                    GSASIIspc.MoveToUnitCell(reftrans)
                    if all(abs(noff - reftrans) < 1.e-5):
                        cctbx.remove(refop)
                        break
                else:
                    assert False, "failed on %s:\n\t %s + %s" % (spcname,mult,noff)

def TestSiteSym (spc, crdlist):
    'compare site symmetries and multiplicities for a specified space group'
    msg = "failed on site sym test for %s" % spc
    (E,S) = GSASIIspc.SpcGroup(spc)
    assert not E, msg

    for t in crdlist:
        symb, m = GSASIIspc.SytSym(t[0],S)
        if symb.strip() != t[2].strip() or m != t[1]:
            print spc,t[0],m,symb
        assert m == t[1]
        #assert symb.strip() == t[2].strip()

#import sgtbxout
#key = 'r -3 c'
#CompareWcctbx(key, sgtbxout.sgtbx[key],1)
#exit()

# test #0: exercise MoveToUnitCell
import numpy as np
msg = "MoveToUnitCell failed"
v = [0,1,2,-1,-2]; GSASIIspc.MoveToUnitCell(v); assert v==[0,0,0,0,0], msg
v = np.array([-.1]); GSASIIspc.MoveToUnitCell(v); assert abs(v-0.9) < 1e-6, msg
v = np.array([.1]); GSASIIspc.MoveToUnitCell(v); assert abs(v-0.1) < 1e-6, msg

# test #1: SpcGroup and SGPrint against previous results
import SGdat
for spc in SGdat.SGdat:
    #print spc
    CompareSpcGroup(spc, 0, SGdat.SGdat[spc], SGdat.SGlist[spc] )

# test #2: SpcGroup against cctbx results
import sgtbxdat
for key in sgtbxdat.sgtbx:
    #print key
    CompareWcctbx(key, sgtbxdat.sgtbx[key])

# test #3: exercise SytSym (includes GetOprPtrName, GenAtom, GetKNsym)
# for selected space groups against info in IT Volume A
TestSiteSym('p 1',[
        ((0.13,0.22,0.31),1,'1'),
        ((0,0,0),1,'1'),
])
TestSiteSym('p -1',[
        ((0.13,0.22,0.31),2,'1'),
        ((0,0.5,0),1,'-1'),
])
TestSiteSym('C 2/c',[
        ((0.13,0.22,0.31),8,'1'),
        ((0.0,.31,0.25),4,'2(010)'),
        ((0.25,.25,0.5),4,'-1'),
        ((0,0.5,0),4,'-1'),
])
TestSiteSym('p 2 2 2',[
        ((0.13,0.22,0.31),4,'1'),
        ((0,0.5,.31),2,'2(001)'),
        ((0.5,.31,0.5),2,'2(010)'),
        ((.11,0,0),2,'2(100)'),
        ((0,0.5,0),1,'222'),
])
TestSiteSym('p 4/n',[
        ((0.13,0.22,0.31),8,'1'),
        ((0.25,0.75,.31),4,'2(001)'),
        ((0.5,0.5,0.5),4,'-1'),
        ((0,0.5,0),4,'-1'),
        ((0.25,0.25,.31),2,'4(001)'),
        ((0.25,.75,0.5),2,'-4(001)'),
        ((0.25,.75,0.0),2,'-4(001)'),
])
TestSiteSym('p 31 2 1',[
        ((0.13,0.22,0.31),6,'1'),
        ((0.13,0.0,0.833333333),3,'2(100)'),
        ((0.13,0.13,0.),3,'2(110)'),
])
TestSiteSym('R 3 c',[
        ((0.13,0.22,0.31),18,'1'),
        ((0.0,0.0,0.31),6,'3'),
])
TestSiteSym('R 3 c R',[
        ((0.13,0.22,0.31),6,'1'),
        ((0.31,0.31,0.31),2,'3(111)'),
])
TestSiteSym('P 63 m c',[
        ((0.13,0.22,0.31),12,'1'),
        ((0.11,0.22,0.31),6,'m(100)'),
        ((0.333333,0.6666667,0.31),2,'3m(100)'),
        ((0,0,0.31),2,'3m(100)'),
])
TestSiteSym('I a -3',[
        ((0.13,0.22,0.31),48,'1'),
        ((0.11,0,0.25),24,'2(100)'),
        ((0.11,0.11,0.11),16,'3(111)'),
        ((0,0,0),8,'-3(111)'),
])

print "OK"
