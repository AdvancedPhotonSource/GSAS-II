#!/Users/toby/build/cctbx_build/bin/python
import sys
import cctbx.uctbx as uc
import numpy as np
import datetime

fp = sys.stdout
fp.write("# output from uctbx (cctbx) computed on platform %s on %s\n" %
         (sys.platform, datetime.date.today()) )
fp.write("#import numpy as np\n")
fp.write("array = np.array\n\n")
fp.write("CellTestData = [\n")

for cell in [
    (4,4,4,90,90,90),
    (4.1,5.2,6.3,100,80,130),
    (3.5,3.5,6,90,90,120),
    ]:
    fp.write("# cell, g, G, cell*, V, V*\n")
    result = []
    result.append(cell)
    mm = uc.unit_cell(cell).metrical_matrix()
    result.append(
        np.array([mm[i] for i in (0,3,4,3,1,5,4,5,2)]).reshape(3,3)
        )
    rmm = uc.unit_cell(cell).reciprocal_metrical_matrix()
    result.append(
        np.array([rmm[i] for i in (0,3,4,3,1,5,4,5,2)]).reshape(3,3)
        )
    result.append(
        uc.unit_cell(cell).reciprocal().parameters()
        )
    result.append(
        uc.unit_cell(cell).volume()
        )
    result.append(
        uc.unit_cell(cell).reciprocal().volume()
        )
    fp.write("  %s,\n" % result)
fp.write("]\n")

fp.write("CoordTestData = [\n")
for cell in [
    (4,4,4,90,90,90),
    (4.1,5.2,6.3,100,80,130),
    (3.5,3.5,6,90,90,120),
    ]:
    fp.write("# cell, ((frac, ortho),...)\n")
    fp.write("  ((%s,%s,%s,%s,%s,%s,), [\n" % cell)
    for frac in [
        (0.1,0.,0.),
        (0.,0.1,0.),
        (0.,0.,0.1),
        (0.1,0.2,0.3),
        (0.2,0.3,0.1),
        (0.3,0.2,0.1),
        (0.5,0.5,0.5),
        ]:
        fp.write(" (%s,%s),\n" % (frac, uc.unit_cell(cell).orthogonalize(frac)))
    fp.write("]),\n")
fp.write("]\n")
    
sys.exit()
