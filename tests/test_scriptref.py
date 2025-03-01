# perform a GSAS-II refinement using GSASIIscriptable and tutorial
# data

import os
import sys
import tempfile
import numpy.testing as npt
home = os.path.dirname(__file__)
work = tempfile.gettempdir()

import importlib.util
G2loc = None
try: 
    G2loc = importlib.util.find_spec('GSASII.GSASIIscriptable')
except ModuleNotFoundError:
    print('ModuleNotFound for GSASII.GSASIIscriptable')

if G2loc is None: # fixup path if GSASII not installed into Python
    print('GSAS-II not installed in Python; Hacking sys.path')
    sys.path.append(os.path.dirname(home))

import GSASII
import GSASII.GSASIIscriptable as G2sc

def test_refine():
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
    #h1.set_refinements({'Limits': [16.,158.4]})
    #h2.set_refinements({'Limits': [19.,153.]})
    h1.set_refinements({'Limits': [16.,110]})   # decreasing range speeds this up by ~40%
    h2.set_refinements({'Limits': [19.,120.]})
    #gpx.set_Controls('cycles', 8)
    gpx.set_Controls('cycles', 3)  # also gains ~x1.5 in speed
    h1.set_refinements({"Background": { "no. coeffs": 6, "refine": True }})
    h2.set_refinements({"Background": { "no. coeffs": 3, "refine": True }})
    gpx.refine()
    testR('Fit scale & bkg',45.811562,17.864834)
    #
    phase0.set_refinements({'Cell':True})
    phase0.set_HAP_refinements({'HStrain':True},[h2])
    gpx.refine()
    testR('Fit cells',32.475886, 15.02412)
    #
    phase0.set_HAP_refinements({'Mustrain':{'refine':True}},[h1])
    #phase0.set_HAP_refinements({'Size':{'refine':True}},[h1])
    h1.set_refinements({"Sample Parameters": {"Shift": True}})
    h2.set_refinements({"Sample Parameters":["DisplaceX","DisplaceY"]})
    phase0.set_refinements({"Atoms":{"all":"XU"}})
    gpx.refine()
    testR('add Mustrain, Shift, Displace[XY], atomic X & Uiso',
              13.407161,  6.360408)
    #
    h1.set_refinements({'Instrument Parameters': ['U', 'V', 'W']})
    h2.set_refinements({'Instrument Parameters': ['U', 'V', 'W']})
    gpx.refine()
    testR('add UVW',10.785432,  4.130126)
    # change to Spherical Harmonics, order=2 for the 1st histogram & refine
    phase0.HAPvalue('PO',2,[h1])
    phase0.set_HAP_refinements({"Pref.Ori.":True})
    gpx.refine()
    POdict = phase0.HAPvalue('PO',targethistlist=[h1])[5]
    print('Spherical harmonics values:',POdict)
    npt.assert_allclose((POdict['C(2,0)'],POdict['C(2,2)']),
                            [0.1171084051086,0.11462063648716], rtol=0.001)
    testR('add PO',10.496639, 4.128754)
    #
    print('OK')
    
if __name__ == '__main__':
    import time
    start = time.time()
    test_refine()
    print('elapsed=',time.time()-start)
    
