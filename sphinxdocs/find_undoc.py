'''Find Python files in GSASII that are not referenced in the sphinx documentation
'''
import glob
import os.path
import subprocess as sp

# get list of documented files (misses out on data files -- with no functions or classes)
loc = os.path.split(os.path.realpath(__file__))[0]
documented = [os.path.splitext(os.path.split(fil)[1])[0] for 
              fil in glob.iglob(os.path.join(loc,'build/html/_modules/','*.html'))]
                      
# loop over python files in subversion
proc = sp.Popen(["svn","list",os.path.join(loc,'..')],stdout=sp.PIPE)
undoc = []
for fil in proc.stdout.readlines():
    fil = fil.strip()
    print fil+'...',
    if os.path.splitext(fil.strip())[1] != ".py": continue
    if os.path.splitext(os.path.split(fil)[1])[0] in documented:
        print 'doc'
        continue
    else:
        print '\n'+fil+' undocumented'
        undoc.append(fil)
for fil in undoc:
    print ".. automodule:: "+os.path.splitext(os.path.split(fil)[1])[0]
    print "    :members: "
    print ""
