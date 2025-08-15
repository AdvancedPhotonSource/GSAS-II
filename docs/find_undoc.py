"""Find Python files in GSASII that are not referenced in the sphinx
documentation. Git-reorg version.
"""

import glob
import os.path

import git

# get list of documented files (misses out on data files -- with no functions or classes)
loc = os.path.dirname(os.path.realpath(__file__))  # where is the documentation?
G2loc = os.path.join(loc, "..")

# prepare a list of the files that show up in the documentation already
# look for a reference in the .rst files
documented = []
for fil in glob.glob(os.path.join(loc, "source/*.rst")):
    with open(fil) as file:
        for line in file.readlines():
            if line.strip().startswith(".. automodule::"):
                documented.append(line.strip().split("::")[1].strip())

undoc = []
G2repo = git.Repo(G2loc)
for entry in G2repo.commit().tree.traverse():
    fil = entry.path
    if not fil.endswith(".py"):
        continue
    if os.path.dirname(fil) not in (
        "GSASII",
        "GSASII/imports",
        "GSASII/exports",
        "GSASII/install",
    ):
        continue
    if os.path.splitext(os.path.split(fil)[1])[0] in documented:
        # print(fil+' doc')
        continue
    else:
        print(fil + " undocumented")
        undoc.append(fil)
# generate code to place in a .rst file
if undoc:
    print("\n# place this code somewhere in the .rst files\n#")
for fil in undoc:
    print(".. automodule:: " + os.path.splitext(os.path.split(fil)[1])[0])
    print("    :members: ")
    print("")
