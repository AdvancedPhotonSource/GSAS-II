'''
test_gui_scenarios.py
=====================
Runs each GUI scenario in demos/ (e.g. the CWCombined
tutorial) as a self-test: GSAS-II is started in a subprocess with
``python -m GSASII.GSASIIGUIdriver SCENARIO`` and the scenario must pass.
A scenario that sets ``SELF_TEST = False`` is not run (it remains a GUI
demo; see also ``GUI_DEMO = False`` for the opposite).

These tests need a display and take a while, so they are skipped unless
the environment variable G2_GUI_TESTS is set to 1 (or to "require", which
makes a missing display or wxPython an error instead of a skip, as in the
gui-selftest GitHub workflow). On Linux without a display, use a virtual
one::

    G2_GUI_TESTS=1 xvfb-run -a python -m pytest tests/test_gui_scenarios.py

Other settings (environment variables):

* G2_GUI_TIMEOUT_SCALE: multiply the driver's timeouts by this (e.g. 3)
* G2_GUI_TEST_TIMEOUT: seconds allowed for one scenario (default 900)
* G2_GUI_ARTIFACTS: directory for each scenario's log and, if it fails,
  screenshots and a report (default: pytest's temporary directory)
'''
import os
import sys
import glob
import subprocess
import importlib.util

import pytest

home = os.path.dirname(os.path.abspath(__file__))
repo = os.path.dirname(home)

def _in_selftest(path):
    '''False if the scenario file sets ``SELF_TEST = False`` (a demo that
    is not run as a self-test, e.g. the simulation tutorial). The file is
    parsed, not run.'''
    import ast
    try:
        with open(path, encoding='utf-8') as fp:
            tree = ast.parse(fp.read(), path)
    except (OSError, SyntaxError, ValueError):
        return True         # run it, so that the problem is reported
    value = True
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == 'SELF_TEST'
                and isinstance(node.value, ast.Constant)):
            value = bool(node.value.value)
    return value

SCENARIOS = sorted(p for p in glob.glob(os.path.join(repo, 'demos',
                                                     '*.py'))
                   if not os.path.basename(p).startswith('_')
                   and _in_selftest(p))

def _gsas_path():
    '''Directory to put on PYTHONPATH so that the subprocess imports the
    same GSASII package as this one would: an installed GSAS-II if there
    is one, otherwise this source tree'''
    try:
        spec = importlib.util.find_spec('GSASII')
    except (ImportError, ValueError):
        spec = None
    if spec is None or not spec.origin:
        return repo
    return os.path.dirname(os.path.dirname(spec.origin))

def _skip_reason():
    if os.environ.get('G2_GUI_TESTS', '') in ('', '0'):
        return 'GUI scenario tests are off (set G2_GUI_TESTS=1)'
    if importlib.util.find_spec('wx') is None:
        return 'wxPython is not installed'
    if (sys.platform.startswith('linux') and not os.environ.get('DISPLAY')
            and not os.environ.get('WAYLAND_DISPLAY')):
        return 'no display (use xvfb-run)'
    return None

@pytest.mark.parametrize('scenario', SCENARIOS,
                         ids=[os.path.splitext(os.path.basename(p))[0]
                              for p in SCENARIOS])
def test_gui_scenario(scenario, tmp_path):
    'Run one GUI scenario in test mode; it must report PASS'
    reason = _skip_reason()
    if reason and os.environ.get('G2_GUI_TESTS') == 'require':
        pytest.fail(reason)     # CI: a missing display is an error
    elif reason:
        pytest.skip(reason)
    name = os.path.splitext(os.path.basename(scenario))[0]
    artifacts = os.environ.get('G2_GUI_ARTIFACTS') or str(tmp_path)
    os.makedirs(artifacts, exist_ok=True)
    timeout = float(os.environ.get('G2_GUI_TEST_TIMEOUT') or 900)
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(
        [_gsas_path()] + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
    env['PYTHONUNBUFFERED'] = '1'
    cmd = [sys.executable, '-m', 'GSASII.GSASIIGUIdriver',
           '--artifacts', artifacts, scenario]
    try:
        proc = subprocess.run(cmd, cwd=str(tmp_path), env=env,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=timeout)
        out, status = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as err:
        out, status = err.stdout or b'', f'killed after {timeout:g} s'
    out = out.decode('utf-8', 'replace')
    with open(os.path.join(artifacts, name + '.log'), 'w',
              encoding='utf-8') as fp:
        fp.write(out)
    print(out)     # shown by pytest when the test fails
    tail = '\n'.join(out.splitlines()[-60:])
    assert status == 0 and '[GUIdriver] PASS' in out, (
        f'scenario {name} failed (exit status {status}); '
        f'details in {artifacts}\n{tail}')
