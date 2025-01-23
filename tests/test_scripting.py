# perform a GSAS-II refinement using GSASIIscriptable and tutorial
# data

import os
import sys
import tempfile
home = os.path.dirname(__file__)
work = tempfile.gettempdir()

import importlib  # fixup path if GSASII not installed into Python
if importlib.util.find_spec('GSASII') is None:
    print('Path hacking!')
    os.environ["GSASII_YOLO_PATH"] = "True"
    sys.path.append(os.path.dirname(home))

import GSASII.GSASIIscriptable as G2sc

dataloc = lambda fil: os.path.join(home,'testinp',fil)
workloc = lambda fil: os.path.join(work,fil)
gpx = G2sc.G2Project(newgpx=workloc('test_scripting.gpx'))
# setup step 1: add two histograms to the project
hist1 = gpx.add_powder_histogram(dataloc("PBSO4.XRA"),dataloc("INST_XRY.PRM"),
                                     fmthint='GSAS powder')
hist2 = gpx.add_powder_histogram(dataloc("PBSO4.CWN"),dataloc("inst_d1a.prm"),
                                     fmthint='GSAS powder')
# setup step 2: add a phase and link it to the previous histograms
phase0 = gpx.add_phase(dataloc("PbSO4-Wyckoff.cif"),
                      phasename="PbSO4", histograms=[hist1,hist2])
gpx.set_Controls('cycles', 0)
gpx.refine()
print(hist1.residuals['wR'])
print(hist2.residuals['wR'])
gpx.set_Controls('cycles', 8)
refdict0 = {"set": {"Background": { "no. coeffs": 3, "refine": True }}}
gpx.do_refinements([refdict0])
print(hist1.residuals['wR'])
print(hist2.residuals['wR'])
