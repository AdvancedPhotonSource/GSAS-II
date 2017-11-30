'''Find Python files in GSASII that are not referenced in the sphinx documentation
'''
import glob
import os.path
import subprocess as sp

# get list of documented files (misses out on data files -- with no functions or classes)
import glob
loc = os.path.split(os.path.realpath(__file__))[0] # where is the documentation?

# prepare a list of the files that show up in the documentation already
# look for a reference in the .rst files
documented = []
for fil in glob.glob('source/*.rst'):
    for line in open(fil,'r').readlines():
        if line.strip().startswith('.. automodule::'):
            documented.append(line.strip().split('::')[1].strip())

# old method looks for an html version of the source
# documented = [os.path.splitext(os.path.split(fil)[1])[0] for 
#               fil in glob.iglob(os.path.join(loc,'build/html/_modules/','*.html'))]
                      
# loop over python files referenced in subversion
proc = sp.Popen(["svn","list",
                 os.path.join(loc,'..'),
                 os.path.join(loc,'..','exports'),
                 os.path.join(loc,'..','imports'),
                 ],stdout=sp.PIPE)
undoc = []
for fil in proc.stdout.readlines():
    fil = fil.strip()
    #print fil+'...',
    if os.path.splitext(fil.strip())[1] != ".py": continue
    if os.path.splitext(os.path.split(fil)[1])[0] in documented:
        #print fil+' doc'
        continue
    else:
        print(fil+' undocumented')
        undoc.append(fil)
# generate code to place in a .rst file
if undoc:
    print("\n# place this code somewhere in the .rst files\n#")
for fil in undoc:
    print(".. automodule:: "+os.path.splitext(os.path.split(fil)[1])[0])
    print("    :members: ")
    print("")
