#!/usr/bin/env /Users/toby/build/cctbx_build/bin/python
import datetime
import sys

from cctbx import sgtbx

# generate input for test 7 in GSASIIlattice
dmin = 3
with open("sgtbxlattinp.py", "w") as fp:
    fp.write(
        f"# output from sgtbx computed on platform {sys.platform} on {datetime.date.today()}\n"
    )
    fp.write(f"dmin = {dmin}\n")
    fp.write("sgtbx7 = {\n")
    for sg in (
        "P-1",  # space groups for all 14 Bravais lattices
        "P2/m",
        "C2/m",
        "Pmmm",
        "Cmmm",
        "Fmmm",
        "Immm",
        "P4/mmm",
        "I4/mmm",
        "R-3m",
        "P6/mmm",
        "Pm-3m",
        "Im-3m",
        "Fm-3m",
    ):
        sgi = sgtbx.space_group_info(symbol=sg)
        cs = sgi.any_compatible_crystal_symmetry(volume=1000)
        ms = cs.build_miller_set(anomalous_flag=False, d_min=dmin)
        # ms.show_summary()
        spcg = (f"{ms.space_group_info()}").split(":")[0]
        print(spcg)
        fp.write(f"'{spcg}': [\n")
        fp.write(f"{ms.unit_cell()} ,\n")
        hkllist = {}
        for hkl, d in ms.d_spacings():
            if hkllist.has_key(d):
                hkllist[d].append(hkl)
                print(hkllist[d])
            else:
                hkllist[d] = [
                    hkl,
                ]
        fp.write("  (\n")
        for d in sorted(hkllist.keys(), reverse=True):
            fp.write(f"    ({hkllist[d]}, {d}),\n")
        fp.write("  )],\n\n")
    fp.write("}\n\n")

    # generate input for test 8 in GSASIIlattice
    spg = (
        "P -1",
        "P 2/m",
        "C 2/m",
        "B 1 1 2/m",
        "B 2/m 1 1",
        "I 2/m",
        "P m m m",
        "A m m m",
        "B m m m",
        "C m m m",
        "F m m m",
        "I m m m",
        "P -4",
        "P 4/m m m",
        "I 4/m m m",
        "P 3",
        "P 3 2 1",
        "P 3 1 m",
        "P -3 1 m",
        "R 3",
        "R 3 m",
        "R -3 m",
        "R 3 R",
        "R -3 m R",
        "P -6",
        "P 6/m m m",
        "P m -3 m",
        "I m -3 m",
        "F m -3 m",
    )  # test w/space groups with only centering extinctions
    # N.B. (cctbx does not accept F -1, etc)

    fp.write("sgtbx8 = {\n")
    for sg in spg:
        sgi = sgtbx.space_group_info(symbol=sg)
        cs = sgi.any_compatible_crystal_symmetry(volume=1000)
        ms = cs.build_miller_set(anomalous_flag=False, d_min=dmin)
        # ms.show_summary()
        spcglist = (f"{ms.space_group_info()}").split(":")
        spcg = spcglist[0]
        if len(spcglist) > 1:
            if spcglist[1] == "R":
                spcg += " R"
        print(spcg)
        # fp.write("'%s': [\n" % spcg)
        fp.write(f"'{sg}': [\n")
        fp.write(f"{ms.unit_cell()} ,\n")
        fp.write("  (\n")
        for hkl, d in ms.d_spacings():
            fp.write(f"    ({hkl}, {d}),\n")
        fp.write("  )],\n\n")
    fp.write("}\n")
