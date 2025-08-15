"""
test_tofref.py
===============
Perform a TOF GSAS-II refinement using GSASIIscriptable and tutorial
data. This tests the TOF diffraction code and that Marquardt damping works.
"""

import os
import sys
import tempfile

import numpy.testing as npt

home = os.path.dirname(__file__)
work = tempfile.gettempdir()

import importlib.util

G2loc = None
try:
    G2loc = importlib.util.find_spec("GSASII.GSASIIscriptable")
except ModuleNotFoundError:
    print("ModuleNotFound for GSASII.GSASIIscriptable")

if G2loc is None:  # fixup path if GSASII not installed into Python
    print("GSAS-II not installed in Python; Hacking sys.path")
    sys.path.append(os.path.dirname(home))

import GSASII.GSASIIscriptable as G2sc


def test_refine():
    """Tests a TOF refinement with POWGEN data"""

    def testR(msg, w1):
        print(f"*** {msg}: Rwp(h1)={h1.residuals['wR']:.5f}")
        npt.assert_allclose([h1.residuals["wR"]], [w1], rtol=0.0001)

    print("test_refine(): test a small TOF refinement")

    def dataloc(fil):
        return os.path.join(home, "testinp", fil)

    def workloc(fil):
        return os.path.join(work, fil)

    gpx = G2sc.G2Project(newgpx=workloc("test_scripting.gpx"))
    # setup step 1: add two phases from a original GSAS .EXP file on the web
    URL = "https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW%20Joint%20Refinement/data/NAC.cif"
    phase0 = gpx.add_phase(URL, phasename="NAC", URL=True)
    URL = "https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW%20Joint%20Refinement/data/CaF2.cif"
    phase1 = gpx.add_phase(URL, phasename="CaF2", URL=True)
    URL = "https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW%20Joint%20Refinement/data/PG3_22048.gsa"
    URLprm = "https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW%20Joint%20Refinement/data/POWGEN_2665.instprm"
    # setup step 2: a TOF histogram to the project
    h1 = gpx.add_powder_histogram(
        URL, URLprm, fmthint="GSAS powder", URL=True, phases="all"
    )
    h1.set_refinements({"Limits": [11000.0, 100000]})
    gpx.set_Controls("cycles", 2)
    h1.set_refinements({"Background": {"no. coeffs": 6, "refine": True}})
    gpx.refine()
    testR("After first refinement", 27.34119)

    phase0.set_HAP_refinements({"Scale": True}, [h1])
    gpx.refine()
    testR("2nd refinement w/Phase fraction", 23.272757)

    phase0.set_refinements({"Cell": True})
    phase1.set_refinements({"Cell": True})
    gpx.refine()
    testR("3rd refinement w/cell", 20.966857)

    phase0.set_HAP_refinements({"Mustrain": {"refine": True}}, [h1])
    gpx.refine()
    testR("4th refinement w/Mstrain", 17.202073)
    print("OK")


if __name__ == "__main__":
    import time

    start = time.time()
    test_refine()
    print("elapsed=", time.time() - start)
