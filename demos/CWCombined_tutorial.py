'''GSAS-II GUI scenario: the "Combined X-ray/CW-neutron refinement" tutorial
(PbSO4), https://advancedphotonsource.github.io/GSAS-II-tutorials/CWCombined/Combined%20refinement.htm

Self-test:  python -m GSASII.GSASIIGUIdriver demos/CWCombined_tutorial.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/CWCombined_tutorial.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/CWCombined_tutorial.py
Add --keep-open to explore the finished project afterwards.

The four input files are copied from tests/testinp when this file is in a
GSAS-II source tree, otherwise downloaded from the tutorials web site
(G2_GUI_DOWNLOAD=1 forces the download), into a temporary directory
(d.tmpdir). The driver's file-dialog stand-in gives GSAS-II their local
paths (the GUI's file dialogs cannot open URLs). The project is saved
there too (NXPbSO4.gpx).

This scenario is also a self-test: tests/test_gui_scenarios.py (and the
gui-selftest GitHub workflow) run it and require it to pass.

The GSAS-II tasks (import_powder, add_phase, set_atom, refine, ...) are
methods of the driver; see GSASII/GSASIIGUIdriver.py.

The tutorial is followed through its six refinements; each one's Rwp
values are logged next to the tutorial's. Not scripted: the optional plot
key presses, moving the data limits (no values are given) and the Fourier
maps at the end.
'''
import os
import wx
import wx.grid

# the tutorial this demo follows (its directory in GSASII/tutorialIndex.py);
# Help > Run a GUI demo... offers to open its web page
TUTORIAL = 'CWCombined'

URL = ('https://advancedphotonsource.github.io/GSAS-II-tutorials/'
       'CWCombined/data/')
FILES = ['PBSO4.CWN', 'inst_d1a.prm', 'PBSO4.XRA', 'INST_XRY.PRM']
# the same files are in tests/testinp of a GSAS-II source tree; use them
# when present, so the self-test does not need the network
TESTINP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '..', 'tests', 'testinp')
PHASE = 'PbSO4'
ATOMS = [  # type, x, y, z
    ('Pb', '0.1882', '0.25', '0.167'),
    ('S', '0.063', '0.25', '0.686'),
    ('O', '-0.095', '0.25', '0.600'),
    ('O', '0.181', '0.25', '0.543'),
    ('O', '0.085', '0.026', '0.806'),
]
# Rwp (%) reported by the tutorial after each refinement:
# (neutron, x-ray, total); None where the tutorial gives no value
TUTORIAL = [(18, 41, None), (20, 38, None), (15, 31, None), (15, 26, None),
            (7.05, 12.09, None), (4.53, 10.99, 6.71)]

def refine(d, hists, results, gpx=None):
    """Refine (see GUIDriver.refine); log and show the Rwp values
    (neutron, x-ray, total) next to the tutorial's"""
    n = len(results) + 1
    res = d.refine(gpx=gpx)
    rw = [res['wR'][h] for h in hists] + [res['Rwp']]
    results.append(rw)
    fmt = lambda v: '  -  ' if v is None else f'{v:5.2f}'
    tut = TUTORIAL[n - 1]
    d.detail(f'Rwp: neutron {rw[0]:.2f}%, x-ray {rw[1]:.2f}%, '
             f'total {rw[2]:.2f}%   (tutorial: {fmt(tut[0])}, {fmt(tut[1])}'
             f', {fmt(tut[2])})')
    d.check(all(v > 0 for v in rw), f'refinement {n} done: Rwp neutron '
            f'{rw[0]:.2f}%, x-ray {rw[1]:.2f}%, total {rw[2]:.2f}%')
    d.pause()
    return rw

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'NXPbSO4.gpx')

    d.step('Get the tutorial data files',wait=False)
    files = d.get_files(FILES, URL, [TESTINP])
    d.check(all(os.path.getsize(p) > 500 for p in files.values()),
            f'got {", ".join(FILES)} in {d.tmpdir}')

    # ---- Step 1: read the data ------------------------------------------
    # (wait=False: in step mode, go straight on from the download)
    d.step('Read the neutron data: Import > Powder Data > from GSAS powder '
           'data file, PBSO4.CWN with instrument file inst_d1a.prm',
           wait=True)
    reader = 'from GSAS powder data file'
    new = d.import_powder(files['PBSO4.CWN'], files['inst_d1a.prm'], reader)
    d.check(len(new) == 1, f'one histogram read from PBSO4.CWN: {new}')
    hN = new[0]
    d.step('Read the X-ray data: PBSO4.XRA with instrument file '
           'INST_XRY.PRM (separately, so each gets its own instrument file)')
    new = d.import_powder(files['PBSO4.XRA'], files['INST_XRY.PRM'], reader)
    d.check(len(new) == 1, f'one histogram read from PBSO4.XRA: {new}')
    hX = new[0]
    hists = [hN, hX]
    instN = d.tree_data(hN, 'Instrument Parameters')[0]['Type'][0]
    instX = d.tree_data(hX, 'Instrument Parameters')[0]['Type'][0]
    d.check(instN == 'PNC' and instX == 'PXC',
            f'neutron and X-ray histograms ({instN}, {instX})')

    # ---- Step 2: the phase ----------------------------------------------
    d.step(f'Add a phase named {PHASE}: Data > Add new phase')
    d.add_phase(PHASE)

    d.step('General tab: space group P n m a')
    d.set_space_group('P n m a')
    d.check(True, 'space group is P n m a')

    d.step('General tab: unit cell a = 8.480, b = 5.398, c = 6.958 Å')
    d.set_cell(a='8.480', b='5.398', c='6.958')
    d.check(True, 'cell is 8.480, 5.398, 6.958')

    d.step('Atoms tab: Edit Atoms > Append atom, five times')
    d.append_atoms(len(ATOMS))

    for row, (elem, x, y, z) in enumerate(ATOMS):
        d.step(f'Atom {row+1}: {elem} at {x}, {y}, {z} (double-click the '
               'Type cell to pick the element, then type x, y, z)')
        d.set_atom(row, elem, (x, y, z))
    atoms = d.phase_data()['Atoms']
    d.check([a[1] for a in atoms] == [a[0] for a in ATOMS],
            'atoms are ' + ', '.join(f'{a[0]} {a[1]}' for a in atoms))

    # ---- Step 3: link the histograms ------------------------------------
    d.step(f'{d.hap_location()}: Edit Phase > Add powder histograms, '
           'Set All, OK')
    d.add_histograms()
    d.check(sorted(d.phase_data()['Histograms']) == sorted(hists),
            'both histograms are used by the phase')

    results = []
    # ---- Step 4: refinement 1 -------------------------------------------
    d.step('Refinement 1: Calculate > Refine, saving as NXPbSO4.gpx '
           '(refines scale and background)')
    refine(d, hists, results, gpx=gpx)
    d.check(os.path.exists(gpx), f'project saved as {gpx}')

    # ---- Step 5: refinement 2 -------------------------------------------
    d.step('Refinement 2: General tab, tick "Refine unit cell", refine')
    d.phase_tab('General')
    d.set_check('Refine unit cell')
    d.check(d.phase_data()['General']['Cell'][0],
            'unit cell refinement flag set')
    refine(d, hists, results)

    # ---- Step 6: refinement 3 -------------------------------------------
    d.step(f'Refinement 3: {d.hap_location()}, neutron histogram: tick D11, '
           'D22, D33 (hydrostatic strain), refine')
    d.show_hap()
    d.select_hap_histogram(hN)
    for lbl in ('D11', 'D22', 'D33'):
        d.set_check(lbl)
    d.check(d.phase_data()['Histograms'][hN]['HStrain'][1][:3] == [True] * 3,
            'D11, D22, D33 flags set for the neutron histogram')
    refine(d, hists, results)

    # ---- Step 7: refinement 4 -------------------------------------------
    d.step(f'Refinement 4: {d.hap_location()}, X-ray histogram: tick size '
           'and microstrain, refine')
    d.show_hap()
    d.select_hap_histogram(hX)
    d.set_check('size(')
    d.set_check(' microstrain')
    hap = d.phase_data()['Histograms'][hX]
    d.check(hap['Size'][2][0] and hap['Mustrain'][2][0],
            'size and microstrain flags set for the X-ray histogram')
    refine(d, hists, results)

    # ---- Step 8: refinement 5 -------------------------------------------
    d.step('Refinement 5: X-ray Sample Parameters: tick Sample displacement')
    d.select_tree(hX, 'Sample Parameters')
    d.check(d.tree_data(hX, 'Sample Parameters')['Type'] == 'Bragg-Brentano',
            'X-ray diffractometer type is Bragg-Brentano')
    d.set_check(' Sample displacement')

    d.step('Neutron Sample Parameters: goniometer radius 650 mm, tick '
           'Sample X and Y displacement')
    d.select_tree(hN, 'Sample Parameters')
    d.check(d.tree_data(hN, 'Sample Parameters')['Type'] == 'Debye-Scherrer',
            'neutron diffractometer type is Debye-Scherrer')
    d.type(d.find_after(' Goniometer radius'), '650', enter=True)
    d.set_check(' Sample X displ')
    d.set_check(' Sample Y displ')
    sN = d.tree_data(hN, 'Sample Parameters')
    d.check(sN['Gonio. radius'] == 650 and sN['DisplaceX'][1]
            and sN['DisplaceY'][1],
            'neutron goniometer radius 650 mm, X and Y displacement flagged')
    d.check(d.tree_data(hX, 'Sample Parameters')['Shift'][1],
            'X-ray sample displacement flagged')

    d.step('Atoms tab: double-click the "refine" column heading, choose '
           'X and U for all atoms')
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
    refine(d, hists, results)

    # ---- Step 9: refinement 6 -------------------------------------------
    for h, kind in ((hN, 'neutron'), (hX, 'X-ray')):
        d.step(f'Refinement 6: {kind} Instrument Parameters: refine U, V, W')
        d.select_tree(h, 'Instrument Parameters')
        for p in 'UVW':
            d.set_check(d.find_after(f' {p} (', wx.CheckBox))
        inst = d.tree_data(h, 'Instrument Parameters')[0]
        d.check(all(inst[p][2] for p in 'UVW'),
                f'{kind} U, V, W refinement flags set')
    d.step('Refinement 6: refine')
    rw = refine(d, hists, results)

    # ---- summary --------------------------------------------------------
    fmt = lambda v: '   - ' if v is None else f'{v:5.2f}'
    lines = ['Rwp (%)   neutron        x-ray          total',
             '          run  tutorial  run  tutorial  run  tutorial']
    for n, (r, t) in enumerate(zip(results, TUTORIAL), 1):
        lines.append(f'refine {n}  ' + '  '.join(
            f'{r[i]:5.2f} {fmt(t[i])}' for i in range(3)))
    print('\n'.join(lines), flush=True)
    d.step('Final Rwp: neutron {:.2f}%, X-ray {:.2f}%, total {:.2f}% '
           '(tutorial: 4.53%, 10.99%, 6.71%)'.format(*rw))
    tot = [r[2] for r in results]
    d.check(all(b <= a + 0.05 for a, b in zip(tot, tot[1:])),
            'total Rwp never increases: ' + ', '.join(f'{v:.2f}' for v in tot))
    d.check(rw[0] < 6 and rw[1] < 13 and rw[2] < 8.5,
            f'final Rwp close to the tutorial (neutron {rw[0]:.2f} < 6, '
            f'X-ray {rw[1]:.2f} < 13, total {rw[2]:.2f} < 8.5)')
    d.select_tree(hX)
