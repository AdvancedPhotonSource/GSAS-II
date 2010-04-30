'''
 this generates space group tables for all 230 space groups in 
 standard settings using the GSASIIspc.SpcGroup routine. A number 
 of redundant tests are computed, so more than 230 test cases are generated.

 The output from this is placed in spctestinp.py which contains two dictionaries,
 SGdat and SGlist that can be used for testing.

'''
import sys
import datetime
sys.path.append('..')
import GSASIIspc

sgdat = {}
sglist = {}
def GenSGdat(spc):
    (E,D) = GSASIIspc.SpcGroup(spc)
    if E: 
        print >> sys.stderr, "error on: ",spc
    else:
        sgdat[spc] = D
        sglist[spc] = GSASIIspc.SGPrint(D)

duplist = []
######################################################################
#cubic (36)
######################################################################
GenSGdat('p 2 3')
GenSGdat('f 2 3')
GenSGdat('i 2 3')
GenSGdat('p 21 3')
GenSGdat('i 21 3')
GenSGdat('p m 3')
GenSGdat('p n -3')
GenSGdat('f m 3')
GenSGdat('f d 3')
GenSGdat('i m 3')
GenSGdat('p a 3')
GenSGdat('i a 3')
GenSGdat('p 4 3 2')
GenSGdat('p 42 3 2')
GenSGdat('f 4 3 2')
GenSGdat('f 41 3 2')
GenSGdat('i 4 3 2')
GenSGdat('p 43 3 2')
GenSGdat('p 41 3 2')
GenSGdat('i 41 3 2')
GenSGdat('p -4 3 m')
GenSGdat('f -4 3 m')
GenSGdat('i -4 3 m')
GenSGdat('p -4 3 n')
GenSGdat('f -4 3 c')
GenSGdat('i -4 3 d')
GenSGdat('p m -3 m')
GenSGdat('p n -3 n')
GenSGdat('p m -3 n')
GenSGdat('p n -3 m')
GenSGdat('f m -3 m')
GenSGdat('f m -3 c')
GenSGdat('f d -3 m')
GenSGdat('f d -3 c')
GenSGdat('i m -3 m')
GenSGdat('i a -3 d')
# duplicates (IT A naming)
GenSGdat('p m -3') # dup: as before
duplist.append(('p m 3', 'p m -3'))
GenSGdat('f m -3') # dup: as before
duplist.append(('f m 3','f m -3'))
GenSGdat('f d -3') # dup: as before
duplist.append(('f d 3','f d -3'))
GenSGdat('i m -3') # dup: as before
duplist.append(('i m 3','i m -3'))
GenSGdat('p a -3') # dup: as before
duplist.append(('p a 3','p a -3'))
GenSGdat('i a -3') # dup: as before
duplist.append(('i a 3','i a -3'))
######################################################################
# ortho (59)
######################################################################
GenSGdat('p 2 2 2')
GenSGdat('p 2 2 21')
GenSGdat('p 21 21 2')
GenSGdat('p 21 21 21')
GenSGdat('c 2 2 21')
GenSGdat('c 2 2 2')
GenSGdat('f 2 2 2')
GenSGdat('i 2 2 2')
GenSGdat('i 21 21 21')
GenSGdat('p m m 2')
GenSGdat('p m c 21')
GenSGdat('p c c 2')
GenSGdat('p m a 2')
GenSGdat('p c a 21')
GenSGdat('p n c 2')
GenSGdat('p m n 21')
GenSGdat('p b a 2')
GenSGdat('p n a 21')
GenSGdat('p n n 2')
GenSGdat('c m m 2')
GenSGdat('c m c 21')
GenSGdat('c c c 2')
GenSGdat('a m m 2')
GenSGdat('a b m 2')
GenSGdat('a m a 2')
GenSGdat('a b a 2')
GenSGdat('f m m 2')
GenSGdat('f d d 2')
GenSGdat('i m m 2')
GenSGdat('i b a 2')
GenSGdat('i m a 2')
GenSGdat('p m m m')
GenSGdat('p n n n')
GenSGdat('p c c m')
GenSGdat('p b a n')
GenSGdat('p m m a')
GenSGdat('p n n a')
GenSGdat('p m n a')
GenSGdat('p c c a')
GenSGdat('p b a m')
GenSGdat('p c c n')
GenSGdat('p b c m')
GenSGdat('p n n m')
GenSGdat('p m m n')
GenSGdat('p b c n')
GenSGdat('p b c a')
GenSGdat('p n m a')
GenSGdat('c m c m')
GenSGdat('c m c a')
GenSGdat('c m m m')
GenSGdat('c c c m')
GenSGdat('c m m a')
GenSGdat('c c c a')
GenSGdat('f m m m')
GenSGdat('f d d d')
GenSGdat('i m m m')
GenSGdat('i b a m')
GenSGdat('i b c a')
GenSGdat('i m m a')
######################################################################
# tetragonal (68)
######################################################################
GenSGdat('p 4')
GenSGdat('p 41')
GenSGdat('p 42')
GenSGdat('p 43')
GenSGdat('i 4')
GenSGdat('i 41')
GenSGdat('p -4')
GenSGdat('i -4')
GenSGdat('p 4/m')
GenSGdat('p 42/m')
GenSGdat('p 4/n')
GenSGdat('p 42/n')
GenSGdat('i 4/m')
GenSGdat('i 41/a')
GenSGdat('p 4 2 2')
GenSGdat('p 4 21 2')
GenSGdat('p 41 2 2')
GenSGdat('p 41 21 2')
GenSGdat('p 42 2 2')
GenSGdat('p 42 21 2')
GenSGdat('p 43 2 2')
GenSGdat('p 43 21 2')
GenSGdat('i 4 2 2')
GenSGdat('i 41 2 2')
GenSGdat('p 4 m m')
GenSGdat('p 4 b m')
GenSGdat('p 42 c m')
GenSGdat('p 42 n m')
GenSGdat('p 4 c c')
GenSGdat('p 4 n c')
GenSGdat('p 42 m c')
GenSGdat('p 42 b c')
GenSGdat('i 4 m m')
GenSGdat('i 4 c m')
GenSGdat('i 41 m d')
GenSGdat('i 41 c d')
GenSGdat('p -4 2 m')
GenSGdat('p -4 2 c')
GenSGdat('p -4 21 m')
GenSGdat('p -4 21 c')
GenSGdat('p -4 m 2')
GenSGdat('p -4 c 2')
GenSGdat('p -4 b 2')
GenSGdat('p -4 n 2')
GenSGdat('i -4 m 2')
GenSGdat('i -4 c 2')
GenSGdat('i -4 2 m')
GenSGdat('i -4 2 d')
GenSGdat('p 4/m m m')
GenSGdat('p 4/m c c')
GenSGdat('p 4/n b m')
GenSGdat('p 4/n n c')
GenSGdat('p 4/m b m')
GenSGdat('p 4/m n c')
GenSGdat('p 4/n m m')
GenSGdat('p 4/n c c')
GenSGdat('p 42/m m c')
GenSGdat('p 42/m c m')
GenSGdat('p 42/n b c')
GenSGdat('p 42/n n m')
GenSGdat('p 42/m b c')
GenSGdat('p 42/m n m')
GenSGdat('p 42/n m c')
GenSGdat('p 42/n c m')
GenSGdat('i 4/m m m')
GenSGdat('i 4/m c m ')
GenSGdat('i 41/a m d')
GenSGdat('i 41/a c d')
# duplicate -- note gives wrong Laue class
#GenSGdat('i 41 1 1') # dup: as before
#duplist.append(('i 41','i 41 1 1'))
#GenSGdat('p 4/n 1 ') # does not work
######################################################################
# triclinic (2)
######################################################################
GenSGdat('p 1')
GenSGdat('p -1')
######################################################################
# monoclinic (13)
######################################################################
GenSGdat('p 2')
GenSGdat('p 21')
GenSGdat('c 2')
GenSGdat('p m')
GenSGdat('p c')
GenSGdat('c m')
GenSGdat('c c')
GenSGdat('p 2/m')
GenSGdat('p 21/m')
GenSGdat('c 2/m')
GenSGdat('p 2/c')
GenSGdat('p 21/c')
GenSGdat('c 2/c')
# duplicates
GenSGdat('c 1 2 1') # dup: as before
duplist.append(('c 2','c 1 2 1'))
GenSGdat('c 1 2/c 1') # dup: as before
duplist.append(('c 2/c','c 1 2/c 1'))
GenSGdat('p 1 2/m 1') # dup: as before
duplist.append(('p 2/m','p 1 2/m 1'))
######################################################################
# trigonal (25)
######################################################################
GenSGdat('p 3')
GenSGdat('p 31')
GenSGdat('p 32')
GenSGdat('r 3')
GenSGdat('p -3')
GenSGdat('r -3')
GenSGdat('p 3 1 2')
GenSGdat('p 3 2 1')
GenSGdat('p 31 1 2')
GenSGdat('p 31 2 1')
GenSGdat('p 32 1 2')
GenSGdat('p 32 2 1')
GenSGdat('r 3 2')
GenSGdat('p 3 m 1')
GenSGdat('p 3 1 m')
GenSGdat('p 3 c 1')
GenSGdat('p 3 1 c')
GenSGdat('r 3 m')
GenSGdat('r 3 c')
GenSGdat('p -3 1 m')
GenSGdat('p -3 1 c')
GenSGdat('p -3 m 1')
GenSGdat('p -3 c 1')
GenSGdat('r -3 m') 
GenSGdat('r -3 c')
# duplicate
GenSGdat('r 3 r') # dup: rhomb setting
GenSGdat('r -3 r') # dup: rhomb setting
GenSGdat('r 3 2 r') # dup: rhomb setting
GenSGdat('r -3 c r') # dup: rhomb setting
GenSGdat('r 3 m r') # dup: rhomb setting
GenSGdat('r 3 c r') # dup: rhomb setting
GenSGdat('r -3 m r') # dup: rhomb setting
GenSGdat('p 32 1 1') # dup: as before
duplist.append(('p 32','p 32 1 1'))
GenSGdat('r 3 2 h') # dup: hex setting
duplist.append(('r 3 2','r 3 2 h'))
GenSGdat('r 3 m h') # dup: hex setting
duplist.append(('r 3 m','r 3 m h'))
######################################################################
# hexagonal (27)
######################################################################
GenSGdat('p 6')
GenSGdat('p 61')
GenSGdat('p 65')
GenSGdat('p 62')
GenSGdat('p 64')
GenSGdat('p 63')
GenSGdat('p -6')
GenSGdat('p 6/m')
GenSGdat('p 63/m')
GenSGdat('p 6 2 2')
GenSGdat('p 61 2 2')
GenSGdat('p 65 2 2')
GenSGdat('p 62 2 2')
GenSGdat('p 64 2 2')
GenSGdat('p 63 2 2')
GenSGdat('p 6 m m')
GenSGdat('p 6 c c')
GenSGdat('p 63 c m')
GenSGdat('p 63 m c')
GenSGdat('p -6 m 2')
GenSGdat('p -6 c 2')
GenSGdat('p -6 2 m')
GenSGdat('p -6 2 c')
GenSGdat('p 6/m m m')
GenSGdat('p 6/m c c')
GenSGdat('p 63/m c m')
GenSGdat('p 63/m m c')
# duplicate
GenSGdat('p 65 1 1') # dup: as before
duplist.append(('p 65','p 65 1 1'))
GenSGdat('p 6/m 1 1') # dup: as before
duplist.append(('p 6/m','p 6/m 1 1'))
######################################################################
# non-standard space group settings
######################################################################G
GenSGdat('p 1 1 2/m') # dup: non-standard
GenSGdat('p 2/m 1 1') # dup: non-standard
GenSGdat('F -1')      # dup: non-standard
GenSGdat('a 2 2 2')   # dup: non-standard

# do a bit of internal consistency checking
import numpy as np
array = np.array
float32=np.float32
# check for internal agreement with duplicates
for key1,key2 in duplist:
    msg = "checking %s against %s" % (key1, key2)
    keys = sgdat[key1].keys()
    assert len(keys) == len(sgdat[key2].keys()), msg
    for key in keys:
        if key == 'SGOps':
            assert len(sgdat[key2][key]) == len(sgdat[key1][key]), msg
            for i in range(len(sgdat[key2][key])):
                assert np.allclose(sgdat[key1][key][i][0],sgdat[key2][key][i][0]), msg
                assert np.allclose(sgdat[key1][key][i][1],sgdat[key2][key][i][1]), msg
        elif key == 'SGCen':
            assert len(sgdat[key2][key]) == len(sgdat[key1][key]), msg
            for i in range(len(sgdat[key2][key])):
                assert np.allclose(sgdat[key1][key][i][0],sgdat[key2][key][i][0]), msg
                assert np.allclose(sgdat[key1][key][i][1],sgdat[key2][key][i][1]), msg
        elif key == 'SpGrp': # expect this to differ
            pass
        else:
            assert sgdat[key1][key] == sgdat[key2][key], msg+': key = '+key


fp = open('spctestinp.py','w')
fp.write("# output from GSASIIspc computed on platform %s on %s\n" %
         (sys.platform, datetime.date.today()) )
fp.write("import numpy as np\n")
fp.write("array = np.array\n")
fp.write("float32=np.float32\n")
fp.write('# testing %s space groups (25 dups/non-standard)\n' % len(sgdat))
fp.write('SGdat = {\n')
for spc in sgdat:
    fp.write('"%s": %s ,\n' % (spc,sgdat[spc],))
fp.write("}\n")

fp.write('SGlist = {\n')
for spc in sgdat:
    fp.write('"%s": %s ,\n' % (spc,sglist[spc],))
fp.write("}\n")
fp.close()
