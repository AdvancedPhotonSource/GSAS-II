'''GSAS-II GUI scenario: the "Fitting laboratory X-ray powder data for
fluoroapatite" tutorial,
https://advancedphotonsource.github.io/GSAS-II-tutorials/LabData/Laboratory%20X.htm

Self-test:  python -m GSASII.GSASIIGUIdriver demos/LabData_tutorial.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/LabData_tutorial.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/LabData_tutorial.py
Add --keep-open to explore the finished project afterwards.

The three input files are copied from tests/testinp when this file is in a
GSAS-II source tree, otherwise downloaded from the tutorials web site
(G2_GUI_DOWNLOAD=1 forces the download), into a temporary directory
(d.tmpdir). The project is saved there too (fap.gpx).

The GSAS-II tasks (import_powder, import_phase, set_background, refine,
...) are methods of the driver; see GSASII/GSASIIGUIdriver.py.

The tutorial is followed through its four refinements; each one's Rwp is
logged next to the tutorial's. Not scripted: the structure drawing at the
end (Step 9).
'''
import os
import wx
import wx.grid

# the tutorial this demo follows (its directory in GSASII/tutorialIndex.py);
# Help > Run a GUI demo... offers to open its web page
TUTORIAL = 'LabData'

URL = ('https://advancedphotonsource.github.io/GSAS-II-tutorials/'
       'LabData/data/')
FILES = ['FAP.XRA', 'INST_XRY.PRM', 'FAP.EXP']
# a GSAS-II source tree has INST_XRY.PRM in tests/testinp
TESTINP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '..', 'tests', 'testinp')
PHASE = 'fap'
# Rwp (%) reported by the tutorial after each refinement
TUTORIAL_RWP = [45.2, 14.5, 11.2, 10.4]

def refine(d, results, gpx=None):
    """Refine (see GUIDriver.refine); log and show the Rwp next to the
    tutorial's"""
    n = len(results) + 1
    res = d.refine(gpx=gpx)
    rw = res['Rwp']
    results.append(rw)
    d.detail(f'Rwp: {rw:.2f}%   (tutorial: {TUTORIAL_RWP[n - 1]}%)')
    d.check(rw > 0, f'refinement {n} done: Rwp {rw:.2f}%')
    d.pause()
    return rw

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'fap.gpx')

    d.step('Get the tutorial data files', wait=False)
    files = d.get_files(FILES, URL, [TESTINP])
    d.check(all(os.path.getsize(p) > 500 for p in files.values()),
            f'got {", ".join(FILES)} in {d.tmpdir}')

    # ---- Step 1: read the data ------------------------------------------
    d.step('Read the X-ray data. Import > Powder Data > from GSAS powder '
           '(menu)\ndata file: FAP.XRA; instrument file: INST_XRY.PRM')
    new = d.import_powder(files['FAP.XRA'], files['INST_XRY.PRM'],
                          'from GSAS powder data file')
    d.check(len(new) == 1, f'one histogram read from FAP.XRA: {new}')
    hist = new[0]
    inst = d.tree_data(hist, 'Instrument Parameters')[0]['Type'][0]
    d.check(inst == 'PXC', f'CW X-ray histogram ({inst})')

    # ---- Step 2: the phase ----------------------------------------------
    d.step('Read the phase: Import > Phase > from GSAS .EXP file (menu), '
           'FAP.EXP; keep the name fap; add the histogram')
    name = d.import_phase(files['FAP.EXP'], 'from GSAS .EXP file')
    d.check(name == PHASE, f'phase {name} read')
    data = d.phase_data()
    sg = data['General']['SGData']['SpGrp']
    d.check(sg.replace(' ', '') == 'P63/m' and len(data['Atoms']) == 7,
            f'space group {sg}, {len(data["Atoms"])} atoms')
    d.check(list(data['Histograms']) == [hist],
            'the histogram is used by the phase')

    # ---- Step 3: look at the phase ---------------------------------------
    d.step('Look at the phase: General tab (2 formula units per cell)')
    d.phase_tab('General')
    d.pause()

    # ---- Step 4: background ----------------------------------------------
    d.step('Background: choose "log interpolate" with 9 coefficients')
    d.set_background(hist, function='log interpolate', terms=9)
    d.check(True, 'background is log interpolate, 9 terms')

    # ---- Step 5: refinement 1 --------------------------------------------
    results = []
    d.step('Refinement 1: Calculate > Refine (menu) (refines scale and '
           'background); save as fap.gpx')
    refine(d, results, gpx=gpx)
    d.check(os.path.exists(gpx), f'project saved as {gpx}')

    # ---- Step 6: refinement 2 --------------------------------------------
    d.step('General tab: check "Refine unit cell"')
    d.phase_tab('General')
    d.set_check('Refine unit cell')
    d.check(d.phase_data()['General']['Cell'][0],
            'unit cell refinement flag set')

    d.step('Sample Parameters: check that the type is Bragg-Brentano; add '
           'Sample displacement')
    d.select_tree(hist, 'Sample Parameters')
    s = d.tree_data(hist, 'Sample Parameters')
    d.check(s['Type'] == 'Bragg-Brentano', 'diffractometer type is '
            'Bragg-Brentano')
    d.set_check(' Sample displacement')
    d.check(d.tree_data(hist, 'Sample Parameters')['Shift'][1],
            'sample displacement flagged')

    d.step('Refinement 2: Calculate > Refine (menu)')
    refine(d, results)

    # ---- Step 7: refinement 3 --------------------------------------------
    d.step(f'{d.hap_location()}: add size and microstrain')
    d.show_hap()
    d.select_hap_histogram(hist)
    d.set_check('size(')
    d.set_check(' microstrain')
    hap = d.phase_data()['Histograms'][hist]
    d.check(hap['Size'][2][0] and hap['Mustrain'][2][0],
            'size and microstrain flags set')

    d.step('Refinement 3: Calculate > Refine (menu)')
    refine(d, results)

    # ---- Step 8: refinement 4 --------------------------------------------
    d.step('Refine atoms on Atoms tab: double-click the "refine" column '
           'heading; choose X and U for all atoms')
    d.grid_event(d.atoms_grid(), -1, 'refine',
                 wx.grid.wxEVT_GRID_LABEL_LEFT_DCLICK, modal=True)
    dlg = d.dialog('Refinement controls')
    d.highlight(d.find(wx.CheckListBox, parent=dlg))
    d.gui(dlg.SetSelections, [1, 2])     # X - coordinates, U - thermal
    d.pause(d.delay * 0.5)
    d.click(dlg, wx.ID_OK)
    d.wait_no_modal()
    d.wait_until(lambda: all(a[2] == 'XU' for a in d.phase_data()['Atoms']),
                 'X and U refinement flags on all atoms')
    d.check(True, 'all atoms refine X and U')

    d.step('Refinement 4: Calculate > Refine (menu)')
    rw = refine(d, results)

    # ---- summary --------------------------------------------------------
    print('Rwp (%)   run    tutorial\n' + '\n'.join(
        f'refine {n}  {r:5.2f}  {t:5.2f}'
        for n, (r, t) in enumerate(zip(results, TUTORIAL_RWP), 1)),
        flush=True)
    d.step(f'Final Rwp: {rw:.2f}% (tutorial: 10.38%)')
    d.check(all(b <= a + 0.05 for a, b in zip(results, results[1:])),
            'Rwp never increases: ' + ', '.join(f'{v:.2f}' for v in results))
    d.check(rw < 11.5, f'final Rwp close to the tutorial ({rw:.2f} < 11.5)')
    d.select_tree(hist)
