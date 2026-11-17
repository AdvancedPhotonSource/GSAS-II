'''GSAS-II GUI scenario: create a phase, set its space group, add an atom,
save, clear and reopen the project.

Self-test:  python -m GSASII.GSASIIGUIdriver demos/phase_basics.py
Demo:       python -m GSASII.GSASIIGUIdriver --demo demos/phase_basics.py
Step-by-step demo (press Next for each step):
            python -m GSASII.GSASIIGUIdriver --step demos/phase_basics.py

Uses the driver's phase tasks (add_phase, set_space_group, phase_tab,
append_atoms), File menu commands and native file and message dialogs
(answered by the driver's stand-ins).

This is a self-test only (tests/test_gui_scenarios.py runs it); it is not
offered by Help > Run a GUI demo...
'''
import os
import wx

# not listed by Help > Run a GUI demo... (see scenario_in_menu)
GUI_DEMO = False

PHASE = 'Test phase'

def scenario(d):
    gpx = os.path.join(d.tmpdir, 'phase_basics.gpx')

    d.step('Add a new phase: Data menu, "Add new phase"')
    d.add_phase(PHASE)
    d.check(d.tree_children('Phases') == [PHASE], 'phase is in the data tree')

    d.step('Set the space group: click the space group button')
    d.set_space_group('F m -3 m')
    sg = d.gui(lambda: d.phase_data()['General']['SGData']['SpGrp'])
    d.check(sg.replace(' ', '') == 'Fm-3m', f'space group is F m -3 m ({sg!r})')

    d.step('Go to the Atoms tab with the "Select tab" menu')
    d.phase_tab('Atoms')
    d.check(d.menu_enabled('Edit Atoms', 'Append atom'),
            'Atoms-tab menu bar is active ("Edit Atoms" menu)')

    d.step('Add an atom: Edit Atoms menu, "Append atom"')
    d.append_atoms(1)
    atoms = d.phase_data()['Atoms']
    d.check(len(atoms) == 1, f'one atom added ({atoms[0][:4]})')

    d.step('Save the project: File menu, "Save project as..."')
    d.expect_file('Choose GSAS-II project file name', gpx)
    d.menu('File', 'Save project as...')
    d.check(os.path.exists(gpx), f'project saved to {gpx}')

    d.step('Clear the project: File menu, "New project"')
    d.expect_message('Save before Clear?', wx.ID_NO)
    d.menu('File', 'New project')
    d.check(d.tree_data('Phases', PHASE) is None, 'project cleared')

    d.step('Reopen it: File menu, "Open project..."')
    d.expect_file('Choose GSAS-II project file', gpx)
    d.menu('File', 'Open project...')
    d.wait_until(lambda: d._tree_item(('Phases', PHASE)),
                 'phase after reopening', timeout=30)
    phase = d.tree_data('Phases', PHASE)
    sg = phase['General']['SGData']['SpGrp']
    d.check(sg.replace(' ', '') == 'Fm-3m',
            f'space group survived save and reopen ({sg!r})')
    d.check(len(phase['Atoms']) == 1, 'atom survived save and reopen')

    d.step('Show the reopened phase')
    d.select_tree('Phases', PHASE)
    d.select_tab('Atoms')
