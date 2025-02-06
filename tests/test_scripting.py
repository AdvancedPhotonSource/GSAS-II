# perform a GSAS-II refinement using GSASIIscriptable and tutorial
# data

import os
import sys
import tempfile
import numpy.testing as npt
home = os.path.dirname(__file__)
work = tempfile.gettempdir()

skip = False
import importlib.util  # fixup path if GSASII not installed into Python
G2loc = None
try: 
    G2loc = importlib.util.find_spec('GSASII.GSASIIscriptable')
except ModuleNotFoundError:
    print('ModuleNotFound for GSASII.GSASIIscriptable')

if G2loc is None:
    print('GSAS-II not installed in Python: Hacking sys.path')
    os.environ["GSASII_YOLO_PATH"] = "True"
    sys.path.append(os.path.dirname(home))
else:
    # TODO: figure out why this fails with a "path hacking" error from
    # inside the testing suite; Would prefer to skip only from inside pytest
    if "pytest" in sys.modules:
        print("unable to run inside testing suite?!")
        skip = True

import GSASII
import GSASII.GSASIIscriptable as G2sc

def test_refine():
    if skip:
        print('unable to test from pytest')
        return
    # this does not catch the PathHackingException when run from pytest
    try:
        refinetest()
    except GSASII.PathHackingException:
        print('PathHackingException')

def refinetest():
    def testR(msg,w1,w2):
        print(f"*** {msg}: Rwp(h1)={h1.residuals['wR']:.5f}, Rwp(h2)={h2.residuals['wR']:.5f}")
        npt.assert_allclose([h1.residuals['wR'],h2.residuals['wR']],
                            [w1,w2], rtol=0.0001)

    print('test_refine(): test a small refinement')
    dataloc = lambda fil: os.path.join(home,'testinp',fil)
    workloc = lambda fil: os.path.join(work,fil)
    gpx = G2sc.G2Project(newgpx=workloc('test_scripting.gpx'))
    # setup step 1: add two histograms to the project
    h1 = gpx.add_powder_histogram(dataloc("PBSO4.XRA"),dataloc("INST_XRY.PRM"),
                                         fmthint='GSAS powder')
    h2 = gpx.add_powder_histogram(dataloc("PBSO4.CWN"),dataloc("inst_d1a.prm"),
                                         fmthint='GSAS powder')
    # setup step 2: add a phase and link it to the previous histograms
    phase0 = gpx.add_phase(dataloc("PbSO4-Wyckoff.cif"),
                               phasename="PbSO4", histograms=[h1,h2])
    gpx.set_Controls('cycles', 0)
    gpx.refine()
    testR('Before fitting',96.681098,99.748994)
    # 
    h1.set_refinements({'Limits': [16.,158.4]})
    h2.set_refinements({'Limits': [19.,153.]})
    gpx.set_Controls('cycles', 8)
    h1.set_refinements({"Background": { "no. coeffs": 6, "refine": True }})
    h2.set_refinements({"Background": { "no. coeffs": 3, "refine": True }})
    gpx.refine()
    testR('Fit scale & bkg',40.64193551740201,18.6189841945785)
    #
    phase0.set_refinements({'Cell':True})
    phase0.set_HAP_refinements({'HStrain':True},[h2])
    gpx.refine()
    testR('Fit cells',30.072804646662338,15.014744642359773)
    #
    phase0.set_HAP_refinements({'Mustrain':{'refine':True}},[h1])
    #phase0.set_HAP_refinements({'Size':{'refine':True}},[h1])
    h1.set_refinements({"Sample Parameters": {"Shift": True}})
    h2.set_refinements({"Sample Parameters":["DisplaceX","DisplaceY"]})
    phase0.set_refinements({"Atoms":{"all":"XU"}})
    gpx.refine()
    testR('add Mustrain, Shift, Displace[XY], atomic X & Uiso',
              12.66845815113383,6.695761603085025)
    #
    h1.set_refinements({'Instrument Parameters': ['U', 'V', 'W']})
    h2.set_refinements({'Instrument Parameters': ['U', 'V', 'W']})
    gpx.refine()
    testR('add UVW',10.518485346229324,4.495180030877684)
    print('OK')
    
if __name__ == '__main__':
    test_refine()
