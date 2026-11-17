'''Demo for the "Simulating Powder Diffraction with GSAS-II"
tutorial. Simulation of a diffraction data requires input of
a set of instrument parameter (GSAS-II provides default values
if no file is available) and one or more crystal structures, which 
are usually read in, but can be entered manually. Here two phases, 
CuCr2O4 and CuO are read from CIF files and a high-resolution 
synchrotron pattern is simulated. Simulation is also commonly performed 
using scripting, (see https://gsas-ii-scripting.readthedocs.io/).

Tutorial location: https://advancedphotonsource.github.io/GSAS-II-tutorials/Simulation/SimTutorial.htm

The input for this demo is downloaded from the tutorials web site
into a selected directory along with the created project file 
(CuCr2O4-sim.gpx).

The tutorial is followed through Step 6 (the simulation). 
Not included in the demo: rearranging the plot and the optional Step 7.
'''
import os
import wx

'''
To run this directly, use:
Self-test:  python -m GSASII.GSASIIGUIdriver demos/Simulation_tutorial.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/Simulation_tutorial.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/Simulation_tutorial.py
Add --keep-open to explore the finished project afterwards.
'''

# the tutorial this demo follows (its directory in GSASII/tutorialIndex.py);
# Help > Run a GUI demo... offers to open its web page
TUTORIAL = 'Simulation'

URL = ('https://advancedphotonsource.github.io/GSAS-II-tutorials/'
       'Simulation/data/')
FILES = ['CuCr2O4.cif', 'CuO.cif', 'BT1.prm']
PHASES = ['CuCr2O4', 'CuO']
# formula weights (per unit cell) given by the tutorial
MASS = {'CuCr2O4': 1852.3, 'CuO': 318.2}
# CuO phase fraction for 90 wt% CuCr2O4 / 10 wt% CuO
CUO_FRACTION = 0.6468

def wt_fractions(d, hist):
    'Weight fractions of the phases in a histogram (as on the Data tab)'
    w = {p: d.phase_data(p)['Histograms'][hist]['Scale'][0]
         * d.phase_data(p)['General']['Mass'] for p in PHASES}
    total = sum(w.values())
    return {p: v / total for p, v in w.items()}

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'CuCr2O4-sim.gpx')

    d.step('Get the tutorial data files', wait=False)
    files = d.get_files(FILES, URL, [])
    d.check(all(os.path.getsize(p) > 500 for p in files.values()),
            f'got {", ".join(FILES)} in {d.tmpdir}')

    # ---- Step 1: the phases ---------------------------------------------
    for phase in PHASES:
        d.step(f'Read phase {phase}: Import > Phase > from CIF file (menu), '
               f'{phase}.cif; keep the name {phase}')
        name = d.import_phase(files[f'{phase}.cif'], 'from CIF file')
        d.check(name == phase, f'phase {name} read')

    # ---- Step 2: the simulated histograms --------------------------------
    d.step('Simulate an X-ray dataset: Import > Powder Data > Simulate a '
           'dataset (menu); Cancel the file dialog, choose the APS 11-BM '
           'defaults; name "11BM sim", start at 5°; add to both phases')
    xray = d.simulate_powder(default='Defaults for APS 30keV 11BM',
                             name='11BM sim', start=5)
    d.check(xray.startswith('PWDR 11BM sim'), f'histogram {xray}')

    d.step('Simulate a CW neutron dataset: Import > Powder Data > Simulate '
           'a dataset (menu); instrument file BT1.prm; 1-120° in 0.025° '
           'steps; name "BT1 sim"; add to both phases')
    neut = d.simulate_powder(instprm=files['BT1.prm'], name='BT1 sim',
                             start=1, end=120, step=0.025)
    d.check(neut.startswith('PWDR BT1 sim'), f'histogram {neut}')
    x = d.tree_data(neut)[1][0]
    d.check(abs(x[0] - 1) < 0.03 and abs(x[-1] - 120) < 0.03
            and abs(x[1] - x[0] - 0.025) < 1e-4,
            f'{len(x)} points, {x[0]:.3f}-{x[-1]:.3f}°')
    for phase in PHASES:
        d.check(sorted(d.phase_data(phase)['Histograms'])
                == sorted([xray, neut]), f'{phase} is in both histograms')

    # ---- Step 3: phase fractions ----------------------------------------
    d.step('General tabs: look at the formula weights (CuCr2O4 1852.3, CuO '
           '318.2 per cell)')
    for phase in PHASES:
        d.phase_tab('General', phase)
        mass = d.phase_data(phase)['General']['Mass']
        d.check(abs(mass - MASS[phase]) < 0.5,
                f'{phase} formula weight {mass:.1f}')
        d.pause()

    d.step(f'{d.hap_location()} of CuO, {xray}: set the phase fraction to '
           f'{CUO_FRACTION} (10 wt%)')
    d.show_hap('CuO')
    d.select_hap_histogram(xray)
    d.type(d.find_after(' Phase fraction'), str(CUO_FRACTION), enter=True)
    d.settle(0.5)
    d.wait_until(lambda: d.find(wx.StaticText, label=lambda s: s.startswith(
        ' Wt. fraction: 0.100'), wait=False), 'weight fraction 0.100 shown')
    d.check(True, 'CuO weight fraction 10%')

    d.step('Copy it to the neutron histogram: Edit Phase > Copy data '
           '(menu), Set All, OK')
    d.copy_hap(phase='CuO')
    d.wait_until(lambda: d.phase_data('CuO')['Histograms'][neut]['Scale'][0]
                 == CUO_FRACTION, 'copied CuO phase fraction')
    for hist in (xray, neut):
        wt = wt_fractions(d, hist)
        d.check(abs(wt['CuCr2O4'] - 0.9) < 0.001,
                f'{hist}: ' + ', '.join(f'{p} {100 * v:.1f} wt%'
                                        for p, v in wt.items()))

    # ---- Step 4: no optimization ----------------------------------------
    d.step('Controls: set Max cycles to 0 (compute only)')
    d.select_tree('Controls')
    d.choose(d.find_after('Max cycles: ', wx.ComboBox), '0')
    d.settle()
    d.check(d.tree_data('Controls')['max cyc'] == 0, 'Max cycles is 0')

    # ---- Step 5: neutron scale factor ------------------------------------
    d.step(f'Sample Parameters of {neut}: histogram scale factor 1000')
    d.select_tree(neut, 'Sample Parameters')
    d.type(d.find_after(lambda s: s.strip().startswith(
        'Histogram scale factor')), '1000', enter=True)
    d.settle()
    scale = d.tree_data(neut, 'Sample Parameters')['Scale'][0]
    d.check(scale == 1000, f'{neut} scale factor {scale}')

    # ---- Step 6: simulate ------------------------------------------------
    d.step('Simulate the patterns: Calculate > Refine (menu); save as '
           'CuCr2O4-sim.gpx')
    res = d.refine(gpx=gpx)
    d.check(os.path.exists(gpx), f'project saved as {gpx}')
    for hist in (xray, neut):
        y = d.tree_data(hist)[1][1]
        d.check(max(y) > 100, f'{hist}: simulated pattern, max '
                f'{max(y):.0f} counts')

    d.step(f'Look at the simulated patterns: {neut}')
    d.select_tree(neut)
    d.pause()
    d.step(f'... and {xray}')
    d.select_tree(xray)
    d.check(res['Rwp'] >= 0, 'simulation done')
