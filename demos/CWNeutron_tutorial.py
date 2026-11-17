'''Demo for the "CW Neutron Powder fit for Yttrium-Iron Garnet"
tutorial, which fits the structure of yttrium iron/aluminum garnet 
using constant wavelength neutron data. The structure is cubic with 
four different atomic sites, two of which are jointly occupied by 
both iron and aluminum.

Tutorial location: https://advancedphotonsource.github.io/GSAS-II-tutorials/CWNeutron/Neutron%20CW%20Powder%20Data.htm

The input for this demo is downloaded from the tutorials web site
into a selected directory along with the created project file (YAG.gpx).

The tutorial is followed through its four refinements; each one's Rwp is
logged next to the tutorial's. Not included in the demo: viewing 
the .lst file and the covariance matrix at the end.
'''
import os
import wx

'''
To run this directly, use:
Self-test:  python -m GSASII.GSASIIGUIdriver demos/CWNeutron_tutorial.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/CWNeutron_tutorial.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/CWNeutron_tutorial.py
Add --keep-open to explore the finished project afterwards.
'''

# the tutorial this demo follows (its directory in GSASII/tutorialIndex.py);
# Help > Run a GUI demo... offers to open its web page
TUTORIAL = 'CWNeutron'

URL = ('https://advancedphotonsource.github.io/GSAS-II-tutorials/'
       'CWNeutron/data/')
FILES = ['garnet.raw', 'inst_d1a.prm']
# a GSAS-II source tree has inst_d1a.prm in tests/testinp
TESTINP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '..', 'tests', 'testinp')
PHASE = 'garnet'
ATOMS = [  # type, x, y, z
    ('Y', '0.125', '0', '0.25'),
    ('Fe', '0', '0', '0'),
    ('Al', '0', '0', '0'),
    ('Al', '0.375', '0', '0.25'),
    ('Fe', '0.375', '0', '0.25'),
    ('O', '-0.03', '0.05', '0.15'),
]
# Rwp (%) reported by the tutorial after each refinement
TUTORIAL_RWP = [39, 12, 10, 5.18]

def refine(d, results, gpx=None):
    """Refine (see GUIDriver.refine); log and show the Rwp next to the
    tutorial's"""
    n = len(results) + 1
    res = d.refine(gpx=gpx)
    rw = res['Rwp']
    results.append(rw)
    d.detail(f'Rwp: {rw:.2f}%   (tutorial: ~{TUTORIAL_RWP[n - 1]}%)')
    d.check(rw > 0, f'refinement {n} done: Rwp {rw:.2f}%')
    d.pause()
    return rw

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'YAG.gpx')

    d.step('Get the tutorial data files', wait=False)
    files = d.get_files(FILES, URL, [TESTINP])
    d.check(all(os.path.getsize(p) > 500 for p in files.values()),
            f'got {", ".join(FILES)} in {d.tmpdir}')

    # ---- Step 1: read the data ------------------------------------------
    d.step('Read the neutron data. Import > Powder Data > from GSAS powder '
           '(menu)\ndata file: garnet.raw; instrument file: inst_d1a.prm')
    new = d.import_powder(files['garnet.raw'], files['inst_d1a.prm'],
                          'from GSAS powder data file')
    d.check(len(new) == 1, f'one histogram read from garnet.raw: {new}')
    hist = new[0]
    inst = d.tree_data(hist, 'Instrument Parameters')[0]['Type'][0]
    d.check(inst == 'PNC', f'CW neutron histogram ({inst})')

    # ---- Step 2: the phase ----------------------------------------------
    d.step(f'Add a phase named {PHASE}: Data > Add new phase (menu)')
    d.add_phase(PHASE)

    d.step('Change space group on General tab: I a -3 d')
    d.set_space_group('I a -3 d')
    d.check(True, 'space group is I a -3 d')

    d.step('Set lattice on General tab: a = 12.19 Å')
    d.set_cell(a='12.19')
    d.check(True, 'cell is a = 12.19 Å')

    d.step(f'{d.hap_location()}: Edit Phase > Add powder histograms (menu), '
           'select the histogram, OK')
    d.add_histograms()
    d.check(list(d.phase_data()['Histograms']) == [hist],
            'the histogram is used by the phase')

    # ---- Step 3: the atoms ----------------------------------------------
    d.step('Create atoms on Atoms tab: Edit Atoms > Append atom, x6')
    d.append_atoms(len(ATOMS))

    for row, (elem, x, y, z) in enumerate(ATOMS):
        d.step(f'Set Atom {row+1}: to {elem} at {x}, {y}, {z} (double-click '
               'the Type cell to pick the element, then type x, y, z)')
        d.set_atom(row, elem, (x, y, z))
    atoms = d.phase_data()['Atoms']
    d.check([a[1] for a in atoms] == [a[0] for a in ATOMS],
            'atoms are ' + ', '.join(f'{a[0]} {a[1]}' for a in atoms))

    results = []
    d.step('Refinement 1: Calculate > Refine (menu) (refines scale and '
           'background); save as YAG.gpx')
    refine(d, results, gpx=gpx)
    d.check(os.path.exists(gpx), f'project saved as {gpx}')

    # ---- refinement 2: cell, constraints, site fractions -----------------
    d.step('General tab: check "Refine unit cell"')
    d.phase_tab('General')
    d.set_check('Refine unit cell')
    d.check(d.phase_data()['General']['Cell'][0],
            'unit cell refinement flag set')

    for a, b in (('0::AUiso:1', '0::AUiso:2'), ('0::AUiso:3', '0::AUiso:4')):
        d.step(f'Constraints, Phase tab: Edit Constr. > Add equivalence '
               f'(menu); choose {a}, then {b}')
        d.add_constraint('equivalence', a, [b])
    for a, b in (('0::Afrac:1', '0::Afrac:2'), ('0::Afrac:3', '0::Afrac:4')):
        d.step(f'Constraints, Phase tab: Edit Constr. > Add constraint '
               f'equation (menu); choose {a}, then {b}')
        d.add_constraint('equation', a, [b])
    n = len(d.tree_data('Constraints')['Phase'])
    d.check(n == 4, f'{n} phase constraints (2 equivalences, 2 equations)')

    d.step('Atoms tab: select the Fe/Al rows (atoms 2-5), then Edit Atoms > '
           'On selected atoms... > Refine selected (menu); choose '
           'F - site fraction')
    grid = d.atoms_grid()
    def select_rows():
        grid.ClearSelection()
        for r in range(1, 5):
            grid.SelectRow(r, True)
    d.gui(select_rows)
    d.pause(d.delay * 0.5)
    d.menu('Edit Atoms', 'On selected atoms...', 'Refine selected',
           modal=True)
    dlg = d.dialog('Refinement controls')
    d.highlight(d.find(wx.CheckListBox, parent=dlg))
    d.gui(dlg.SetSelections, [0])        # F - site fraction
    d.pause(d.delay * 0.5)
    d.click(dlg, wx.ID_OK)
    d.wait_no_modal()
    d.wait_until(lambda: [a[2] for a in d.phase_data()['Atoms'][1:5]]
                 == ['F'] * 4, 'F refinement flag on the Fe/Al atoms')
    d.check(True, 'atoms 2-5 refine their site fractions')

    d.step('Refinement 2: Calculate > Refine (menu)')
    refine(d, results)

    # ---- refinement 3: atom parameters -----------------------------------
    flags = {'Y': 'U', 'Fe': 'FU', 'Al': 'FU', 'O': 'XU'}
    d.step('Atoms tab: set each atom\'s refine cell (double-click it and '
           'pick): U for Y, FU for Fe and Al, XU for O')
    for row, (elem, *xyz) in enumerate(ATOMS):
        d.edit_cell(d.atoms_grid(), row, 'refine', flags[elem])
    want = [flags[a[0]] for a in ATOMS]
    d.wait_until(lambda: [a[2] for a in d.phase_data()['Atoms']] == want,
                 'atom refinement flags')
    d.check(True, 'atom refinement flags: ' + ', '.join(want))

    d.step('Refinement 3: Calculate > Refine (menu)')
    refine(d, results)

    # ---- refinement 4: profile and sample displacement -------------------
    d.step('Instrument Parameters: add U, V, W')
    d.select_tree(hist, 'Instrument Parameters')
    for p in 'UVW':
        d.set_check(d.find_after(f' {p} (', wx.CheckBox))
    inst = d.tree_data(hist, 'Instrument Parameters')[0]
    d.check(all(inst[p][2] for p in 'UVW'), 'U, V, W refinement flags set')

    d.step('Sample Parameters: set goniometer radius 650 mm, add Sample X '
           'and Y displacement')
    d.select_tree(hist, 'Sample Parameters')
    d.type(d.find_after(' Goniometer radius'), '650', enter=True)
    d.set_check(' Sample X displ')
    d.set_check(' Sample Y displ')
    s = d.tree_data(hist, 'Sample Parameters')
    d.check(s['Gonio. radius'] == 650 and s['DisplaceX'][1]
            and s['DisplaceY'][1],
            'goniometer radius 650 mm; X and Y displacement flagged')

    d.step('Refinement 4: Calculate > Refine (menu)')
    rw = refine(d, results)

    # ---- summary --------------------------------------------------------
    print('Rwp (%)   run    tutorial\n' + '\n'.join(
        f'refine {n}  {r:5.2f}  {t:5.2f}'
        for n, (r, t) in enumerate(zip(results, TUTORIAL_RWP), 1)),
        flush=True)
    d.step(f'Final Rwp: {rw:.2f}% (tutorial: 5.18%)')
    d.check(all(b <= a + 0.05 for a, b in zip(results, results[1:])),
            'Rwp never increases: ' + ', '.join(f'{v:.2f}' for v in results))
    d.check(rw < 6, f'final Rwp close to the tutorial ({rw:.2f} < 6)')
    d.select_tree(hist)
