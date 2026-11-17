'''Demo for the "Combined X-ray/TOF-neutron Rietveld refinement" tutorial,
which fits the structure of Na2Ca3Al2F14 (NAC) to high-resolution 
powder diffraction data collected at both the SNS (POWGEN, TOF neutron)
and the APS (synchrotron, 11-BM) in combination, providing a better
model than would be obtained from either set of data alone. 

Tutorial location: https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW%20Joint%20Refinement/TOF%20combined%20XN%20Rietveld%20refinement%20in%20GSAS.htm

The input for this demo is downloaded from the tutorials web site
into a selected directory along with the created project file 
(NAC.gpx).

The tutorial is followed through Step 5; each refinement's Rwp is logged
next to the tutorial's. The data limits that the tutorial picks on the plot
are typed in the Limits window instead. Not included in the demo: 
Step 6 (further work).
'''
import os
import wx
import wx.grid

'''
To run this directly, use:
Self-test:  python -m GSASII.GSASIIGUIdriver demos/TOFCombined_tutorial.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/TOFCombined_tutorial.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/TOFCombined_tutorial.py
Add --keep-open to explore the finished project afterwards.
'''

# the tutorial this demo follows (its directory in GSASII/tutorialIndex.py);
# Help > Run a GUI demo... offers to open its web page
TUTORIAL = 'TOF-CW Joint Refinement'

URL = ('https://advancedphotonsource.github.io/GSAS-II-tutorials/'
       'TOF-CW%20Joint%20Refinement/data/')
FILES = ['NAC-2015A.EXP', '11BM_NAC.fxye', '11bm_gsas.prm',
         'PG3_22048.gsa', 'POWGEN_1066.instprm',
         'PG3_22049.gsa', 'POWGEN_2665.instprm']
PHASES = ['NAC', 'CaF2']
# Rwp (%) reported by the tutorial after the refinements (None: not given)
TUTORIAL_RWP = [51, 19, None, 8.16, None, 7.822]
# lower TOF limit (μs) for PG3_22048: d ~ 0.52 Å with difC 22600
TOF_MIN = 11750

def refine(d, results, gpx=None):
    """Refine (see GUIDriver.refine); log and show the Rwp next to the
    tutorial's"""
    n = len(results) + 1
    res = d.refine(gpx=gpx)
    rw = res['Rwp']
    results.append(rw)
    tut = TUTORIAL_RWP[n - 1]
    d.detail(f'Rwp: {rw:.2f}%' + (f'   (tutorial: {tut}%)' if tut else ''))
    d.check(rw > 0, f'refinement {n} done: Rwp {rw:.2f}%')
    d.pause()
    return rw

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'NAC.gpx')

    d.step('Get the tutorial data files', wait=False)
    files = d.get_files(FILES, URL, [])
    d.check(all(os.path.getsize(p) > 300 for p in files.values()),
            f'got {len(FILES)} files in {d.tmpdir}')

    # ---- Step 1: the phases ---------------------------------------------
    for phase in PHASES:
        d.step(f'Read phase {phase}: Import > Phase > from GSAS .EXP file '
               f'(menu), NAC-2015A.EXP; choose {phase}, keep the name')
        name = d.import_phase(files['NAC-2015A.EXP'], 'from GSAS .EXP file',
                              choose=phase)
        d.check(name == phase, f'phase {name} read')
    d.check(d.phase_names() == PHASES, f'phases {d.phase_names()}')

    # ---- Step 2: the data -----------------------------------------------
    d.step('Read the X-ray data. Import > Powder Data > from GSAS powder '
           '(menu)\ndata file: 11BM_NAC.fxye; instrument file: '
           '11bm_gsas.prm; add to both phases')
    xray = d.import_powder(files['11BM_NAC.fxye'], files['11bm_gsas.prm'],
                           'from GSAS powder data file')[0]
    d.step('Limits for the X-ray data: Tmin 2.5, Tmax 32 (the tutorial '
           'picks Tmax on the plot)')
    d.set_limits(xray, tmin=2.5, tmax=32)
    d.check(True, f'{xray} limits are 2.5-32°')

    tof = []
    for data, inst in (('PG3_22048.gsa', 'POWGEN_1066.instprm'),
                       ('PG3_22049.gsa', 'POWGEN_2665.instprm')):
        d.step(f'Read TOF neutron data. Import > Powder Data > from GSAS '
               f'powder (menu)\ndata file: {data}; instrument file: {inst}; '
               f'add to both phases')
        new = d.import_powder(files[data], files[inst],
                              'from GSAS powder data file')
        d.check(len(new) == 1, f'one histogram read from {data}: {new}')
        tof.append(new[0])
        if len(tof) == 1:
            d.step(f'Limits for {data}: lower limit at d ~ 0.52 Å '
                   f'(TOF {TOF_MIN} μs)')
            d.set_limits(tof[0], tmin=TOF_MIN)
    hists = [xray] + tof
    for phase in PHASES:
        d.check(sorted(d.phase_data(phase)['Histograms']) == sorted(hists),
                f'{phase} is in all three histograms')

    # ---- Step 3: first refinements ---------------------------------------
    results = []
    d.step('Refinement 1: Calculate > Refine (menu) (refines scales and '
           'background); save as NAC.gpx')
    refine(d, results, gpx=gpx)
    d.check(os.path.exists(gpx), f'project saved as {gpx}')

    for h in tof:
        d.step(f'Background for {h}: 6 coefficients')
        d.set_background(h, terms=6)

    d.step(f'{d.hap_location()} of NAC, X-ray histogram: set microstrain to '
           '100 and refine the phase fraction')
    d.show_hap('NAC')
    d.select_hap_histogram(xray)
    d.type(d.find_after(' microstrain: '), '100', enter=True)
    d.settle()
    d.set_check(' Phase fraction')
    hap = d.phase_data('NAC')['Histograms'][xray]
    d.check(hap['Mustrain'][1][0] == 100 and hap['Scale'][1],
            'NAC microstrain 100, phase fraction flagged')

    d.step('Copy them to the TOF histograms: Edit Phase > Copy data (menu), '
           'Set All, OK')
    d.copy_hap(phase='NAC')
    d.wait_until(lambda: all(d.phase_data('NAC')['Histograms'][h]['Scale'][1]
                             for h in hists), 'copied NAC settings')
    d.check(all(d.phase_data('NAC')['Histograms'][h]['Mustrain'][1][0] == 100
                for h in hists), 'NAC microstrain is 100 in all histograms')

    d.step(f'{d.hap_location()} of CaF2, X-ray histogram: refine the phase '
           'fraction; Edit Phase > Copy flags (menu), Set All, OK')
    d.show_hap('CaF2')
    d.select_hap_histogram(xray)
    d.set_check(' Phase fraction')
    d.copy_hap(flags=True, phase='CaF2')
    d.wait_until(lambda: all(d.phase_data('CaF2')['Histograms'][h]['Scale'][1]
                             for h in hists), 'copied CaF2 flags')
    d.check(True, 'phase fractions of both phases flagged in all histograms')

    for h in range(3):
        d.step(f'Constraints, Histogram/Phase tab: Edit Constr. > Add '
               f'constraint equation (menu); choose 0:{h}:Scale, then '
               f'1:{h}:Scale (phase fractions sum to 1)')
        d.add_constraint('equation', f'0:{h}:Scale', [f'1:{h}:Scale'],
                         tab='Histogram/Phase')
    n = len(d.tree_data('Constraints')['HAP'])
    d.check(n == 3, f'{n} histogram/phase constraints')

    d.step('Refinement 2: Calculate > Refine (menu)')
    refine(d, results)

    # ---- Step 4: nonstructural parameters --------------------------------
    d.step(f'Instrument Parameters of {xray}: add U, V, W')
    d.select_tree(xray, 'Instrument Parameters')
    for p in 'UVW':
        d.set_check(d.find_after(f' {p} (', wx.CheckBox))
    inst = d.tree_data(xray, 'Instrument Parameters')[0]
    d.check(all(inst[p][2] for p in 'UVW'), 'U, V, W refinement flags set')

    d.step(f'Sample Parameters of {xray}: goniometer radius 1000 mm, add '
           'Sample X displacement')
    d.select_tree(xray, 'Sample Parameters')
    d.type(d.find_after(' Goniometer radius'), '1000', enter=True)
    d.set_check(' Sample X displ')
    s = d.tree_data(xray, 'Sample Parameters')
    d.check(s['Gonio. radius'] == 1000 and s['DisplaceX'][1],
            'goniometer radius 1000 mm; X displacement flagged')

    for phase in PHASES:
        d.step(f'General tab of {phase}: check "Refine unit cell"')
        d.phase_tab('General', phase)
        d.set_check('Refine unit cell')
        d.check(d.phase_data(phase)['General']['Cell'][0],
                f'{phase} unit cell refinement flag set')

    for phase in PHASES:
        d.step(f'{d.hap_location()} of {phase}, histogram {tof[0]}: check '
               'D11 (thermal expansion)')
        d.show_hap(phase)
        d.select_hap_histogram(tof[0])
        d.set_check('D11')
        d.check(d.phase_data(phase)['Histograms'][tof[0]]['HStrain'][1][0],
                f'{phase} D11 flagged for {tof[0]}')

    d.step('Refinements 3 and 4: Calculate > Refine (menu), twice')
    refine(d, results)
    refine(d, results)

    # ---- Step 5: structure ----------------------------------------------
    for phase, flags, idx in (('NAC', 'XU', [1, 2]), ('CaF2', 'U', [2])):
        what = ' and '.join({'X': 'X - coordinates', 'U': 'U - thermal'}[f]
                            for f in flags)
        d.step(f'Atoms tab of {phase}: double-click the "refine" column '
               f'heading; choose {what}')
        d.grid_event(d.atoms_grid(phase), -1, 'refine',
                     wx.grid.wxEVT_GRID_LABEL_LEFT_DCLICK, modal=True)
        dlg = d.dialog('Refinement controls')
        d.highlight(d.find(wx.CheckListBox, parent=dlg))
        d.gui(dlg.SetSelections, idx)
        d.pause(d.delay * 0.5)
        d.click(dlg, wx.ID_OK)
        d.wait_no_modal()
        d.wait_until(lambda: all(a[2] == flags
                                 for a in d.phase_data(phase)['Atoms']),
                     f'{flags} refinement flags on the {phase} atoms')
        d.check(True, f'all {phase} atoms refine {flags}')

    d.step('Refinements 5 and 6: Calculate > Refine (menu), twice')
    refine(d, results)
    rw = refine(d, results)

    # ---- summary --------------------------------------------------------
    print('Rwp (%)   run    tutorial\n' + '\n'.join(
        f'refine {n}  {r:5.2f}  {t or "":>5}'
        for n, (r, t) in enumerate(zip(results, TUTORIAL_RWP), 1)),
        flush=True)
    d.step(f'Final Rwp: {rw:.2f}% (tutorial: 7.822%)')
    d.check(all(b <= a + 0.05 for a, b in zip(results, results[1:])),
            'Rwp never increases: ' + ', '.join(f'{v:.2f}' for v in results))
    d.check(rw < 9, f'final Rwp close to the tutorial ({rw:.2f} < 9)')
    d.select_tree(xray)
