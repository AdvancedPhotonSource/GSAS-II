"""
test_image.py
==============
Performs tests that image integration and pixel masking are working.
Also tests the imports that require GSAS-II-compiled binaries.
"""

import os
import sys
import tempfile

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


def test_image():
    """Read in a MAR image (tests pack_f binary module) and then
    integrate it using GeneratePixelMask (tests fmask binary module)
    """
    print("Testing read, mask & integration of MAR image")
    try:
        import requests
    except ModuleNotFoundError:
        print("Module requests not installed, test_image cannot be run")
        return
    workloc = lambda fil: os.path.join(work, fil)
    gpx = G2sc.G2Project(newgpx=workloc("test_image.gpx"))
    img = "https://advancedphotonsource.github.io/GSAS-II-tutorials/2DStrain/data/nx09_strain_011.mar2300"
    img = gpx.add_image(
        img, fmthint="MAR", cacheImage=True, URL=True, download_loc=None
    )
    #    img = gpx.add_image(fmthint='MAR', cacheImage=True,
    #                      imagefile='/tmp/download.mar2300')
    testdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testinp")
    img[0].loadControls(os.path.join(testdir, "nx09_strain_011.imctrl"))
    img[0].GeneratePixelMask()
    masked = sum(img[0].data["Masks"]["SpotMask"]["spotMask"].flatten())
    masked_ratio = masked / img[0].data["Masks"]["SpotMask"]["spotMask"].size
    msg = "Confirm number of masked pixels is about the same"
    assert abs(masked_ratio - 0.02372608695652174) < 0.0001, (
        msg
    )  # + f" {masked_ratio}, 0.02372608695652174"
    # breakpoint()
    plist = img[0].Integrate()
    i = plist[0].getdata("yobs")
    x = plist[0].getdata("x")
    tt1, i1 = (x[233] + x[280]) / 2, sum(i[233:281])
    # print(tt1,i1)
    msg = "Test 1st peak"
    assert abs(tt1 - 2.2695) < 0.0002, msg + f" pos ({tt1})"
    assert abs(i1 - 4656.698709171848) < 0.1, msg + f" intensity ({i1})"
    tt6, i6 = (x[780] + x[850]) / 2, sum(i[780:851])
    # print(tt6,i6)
    msg = "Test 6th peak"
    assert abs(tt6 - 3.9450) < 0.0002, msg + f" pos ({tt6})"
    assert abs(i6 - 26423.739680809293) < 0.1, msg + f" intensity ({i6})"
    print("OK")


def test_CBF():
    "Read in a CBF image (tests unpack_cbf module)"
    print("Testing read of CBF image")
    try:
        import requests
    except ModuleNotFoundError:
        print("Module requests not installed, test_CBF cannot be run")
        return
    workloc = lambda fil: os.path.join(work, fil)
    gpx = G2sc.G2Project(newgpx=workloc("test_CBF.gpx"))
    img = "https://advancedphotonsource.github.io/GSAS-II-tutorials/selftestdata/130mm_0001.cbf"
    img = gpx.add_image(
        img, fmthint="CBF", cacheImage=True, URL=True, download_loc=None
    )
    #    img = gpx.add_image(fmthint='CBF', cacheImage=True,
    #                      imagefile='/tmp/130mm_0001.cbf')
    assert img[0].image.sum() == 184080517, "Sum for entire image"
    assert img[0].image[1000, :].sum() == 96380, "Sum for one column"
    assert img[0].image[:, 1000].sum() == 146435, "Sum for one row"
    print("OK")


if __name__ == "__main__":
    # import GSASII.GSASIIpath as GSASIIpath
    # GSASIIpath.InvokeDebugOpts()
    # import time
    # start = time.time()
    test_image()
    test_CBF()
    # print('elapsed=',time.time()-start)
