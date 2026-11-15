# -*- coding: utf-8 -*-
'''
GSASIIGUIdriver: script the GSAS-II GUI for self-tests and tutorial demos.

A scenario is an ordinary Python function, ``scenario(d)``, that receives a
:class:`GUIDriver` and reads top to bottom::

    def scenario(d):
        d.step('Add a phase')
        d.menu('Data', 'Add new phase', modal=True)  # opens a dialog
        dlg = d.dialog('Phase Name Entry')
        d.type(dlg, 'Test phase')
        d.click(dlg, wx.ID_OK)

The same scenario runs in two modes:

* test mode (default): as fast as possible; steps are logged; a failed
  check ends the run with exit status 1.
* demo mode (``--demo``): each step shows a caption, widgets are outlined
  before they are clicked, text is typed one character at a time, and the
  script pauses between steps so a viewer can follow (or record) it.
* step mode (``--step``, implies ``--demo``): like demo mode, but each
  step's caption has Next and Stop buttons, and the step does not start
  until the viewer presses Next (Enter in the terminal also works). Use
  :meth:`GUIDriver.wait_for_user` for extra stops within a step.

Besides the basic actions (:meth:`GUIDriver.menu`, :meth:`GUIDriver.click`,
:meth:`GUIDriver.type`, ...), the driver has methods for common GSAS-II
tasks, each done through the GUI as a user would do it:
:meth:`~GUIDriver.download`, :meth:`~GUIDriver.import_powder`,
:meth:`~GUIDriver.add_phase`, :meth:`~GUIDriver.phase_tab`,
:meth:`~GUIDriver.set_space_group`, :meth:`~GUIDriver.set_cell`,
:meth:`~GUIDriver.append_atoms`, :meth:`~GUIDriver.set_atom`,
:meth:`~GUIDriver.set_check`, :meth:`~GUIDriver.show_hap`,
:meth:`~GUIDriver.select_hap_histogram`, :meth:`~GUIDriver.add_histograms`
and :meth:`~GUIDriver.refine`. For example::

    def scenario(d):
        d.step('Read the data', wait=False)   # step mode: no stop here
        data, prm = d.download([URL + 'PBSO4.XRA', URL + 'INST_XRY.PRM'])
        hist, = d.import_powder(data, prm)
        d.step('Add the phase')
        d.add_phase('PbSO4')
        d.set_space_group('P n m a')
        d.set_cell(a=8.48, b=5.398, c=6.958)
        ...
        d.add_histograms()
        res = d.refine(gpx=os.path.join(d.tmpdir, 'PbSO4.gpx'))
        d.check(res['Rwp'] < 20, f'Rwp {res["Rwp"]:.2f}')

Run a scenario file (it must define ``scenario(d)``) with::

    python -m GSASII.GSASIIGUIdriver [--demo|--step] [--delay S] [--keep-open] FILE

FILE can also be the name of a scenario in the source tree's demos
directory (``--list`` lists them). For self-tests on slow machines (e.g. CI), ``--timeout-scale
X`` multiplies every timeout by X, and ``--artifacts DIR`` saves
screenshots and a report in DIR if the scenario fails (also settable with
the G2_GUI_TIMEOUT_SCALE and G2_GUI_ARTIFACTS environment variables).
tests/test_gui_scenarios.py runs every scenario this way under pytest.

How it works
------------
GSAS-II is started normally (``GSASIIGUI.main()``). Before that, a
background thread is started; it waits for the main frame and then runs
the scenario. Every GUI action is sent to the GUI thread with
``wx.CallAfter``. While a modal dialog is open, ``ShowModal()`` runs its own
event loop, which still processes those calls, so the script can fill in
and close dialogs that GSAS-II opens.

Native dialogs (``wx.FileDialog``, ``wx.MessageDialog``, ``wx.MessageBox``)
cannot be driven with events, so they are replaced with scripted stand-ins
before GSAS-II starts. GSAS-II's own file browser
(``GSASIIctrlGUI.G2FileBrowser``, used instead of ``wx.FileDialog`` on
Linux) is answered the same way. A scenario queues each answer before the action that
opens the dialog (:meth:`GUIDriver.expect_file`,
:meth:`GUIDriver.expect_message`). A native dialog that was not expected,
or that does not match the queued answer, is recorded as an error and is
answered with Cancel.

So that a run never changes the user's settings, ``G2G.SaveConfigVars`` is
stubbed out, and in test mode the process ends with ``os._exit()`` rather
than ``ExitMain`` (which asks whether to save, and saves config.ini).

Rules for scenarios
-------------------
* Use ``modal=True`` (:meth:`GUIDriver.menu`, :meth:`GUIDriver.click`) for
  any action that opens a wx dialog with ``ShowModal()``; otherwise the call
  waits for the handler, which does not return until the dialog is closed.
  Actions that open only native dialogs (handled by the stand-ins) can be
  called normally.
* Close dialogs by clicking their buttons (:meth:`GUIDriver.click` with a
  button id), not with ``EndModal``: the click runs the dialog's
  validators and OK handler.
* GSAS-II sometimes hides dialogs instead of destroying them, so searches
  are limited to shown windows.
* Do not import GSAS-II GUI modules at the top of a scenario file: they
  must not be imported before the wx.App exists. Import them inside
  ``scenario()`` or use :meth:`GUIDriver.tree_data`.
'''
import os
import sys
import time
import tempfile
import threading
import traceback

import wx
import wx.grid

__all__ = ['GUIDriver', 'ScenarioError', 'ScenarioStopped', 'run', 'main']

class ScenarioError(AssertionError):
    'Raised when a scenario check fails'

class ScenarioStopped(Exception):
    'Raised in step mode when the viewer presses Stop'

def log(*args):
    print('[GUIdriver]', *args, flush=True)

def _env_float(name, default):
    try:
        return float(os.environ.get(name) or default)
    except ValueError:
        return default

#: All waits (:meth:`GUIDriver.wait_until`, :meth:`GUIDriver.gui`, and so
#: every timeout given to a driver method) are multiplied by this factor.
#: Set with ``--timeout-scale`` or the G2_GUI_TIMEOUT_SCALE environment
#: variable, e.g. 3 on a slow CI machine.
TIMEOUT_SCALE = _env_float('G2_GUI_TIMEOUT_SCALE', 1.0)

#==============================================================================
# Scripted stand-ins for native dialogs
#==============================================================================
class Answers:
    '''Queue of answers for native dialogs. Each entry is
    (kind, expected-text, answer, optional); the dialog's message (or
    caption) must contain expected-text, so a wrong dialog is reported, not
    answered. An optional entry is for a dialog that GSAS-II may not show:
    if a different dialog comes first, the optional entry is dropped.
    '''
    def __init__(self):
        self.queue = []
        self.errors = []
        self.lock = threading.Lock()
        self.on_show = None      # called with (kind, text, answer) in demo mode

    def add(self, kind, expect, answer, optional=False):
        '''Queue an answer; returns the entry (for :meth:`discard`)'''
        entry = (kind, expect, answer, optional)
        with self.lock:
            self.queue.append(entry)
        return entry

    def discard(self, entry):
        '''Remove a queued entry if it is still there (an optional dialog
        that was not shown); returns True if it was removed'''
        with self.lock:
            for i, e in enumerate(self.queue):
                if e is entry:
                    del self.queue[i]
                    return True
        return False

    def take(self, kind, text):
        with self.lock:
            while self.queue:
                k, expect, answer, optional = self.queue[0]
                if k == kind and expect in text:
                    break
                if not optional:
                    self.errors.append(f'expected {k} dialog containing '
                                       f'{expect!r}, got {kind} dialog: {text!r}')
                    log(self.errors[-1])
                    return None
                log(f'(optional {k} dialog {expect!r} was not shown)')
                self.queue.pop(0)
            else:
                self.errors.append(f'unexpected {kind} dialog: {text!r}')
                log(self.errors[-1])
                return None
            self.queue.pop(0)
        log(f'{kind} dialog {text.splitlines()[0]!r} -> {_answer_name(answer)}')
        if self.on_show:
            self.on_show(kind, text, answer)
        return answer

def _answer_name(answer):
    'Readable form of a native-dialog answer (button ids by name)'
    if isinstance(answer, int):
        for name in ('ID_YES', 'ID_NO', 'ID_OK', 'ID_CANCEL',
                     'YES', 'NO', 'OK', 'CANCEL'):
            if answer == getattr(wx, name):
                return name.replace('ID_', '').capitalize()
    if isinstance(answer, str):
        return repr(answer)
    return str(answer)

ANSWERS = Answers()
_REAL = {}

class _Inert:
    '''Base for the native-dialog stand-ins: other wx.Dialog methods
    (SetSize, Fit, SetFont, ...) are accepted and do nothing'''
    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return lambda *a, **kw: None

class FakeFileDialog(_Inert):
    'Stand-in for wx.FileDialog; returns a path queued in ANSWERS'
    def __init__(self, parent=None, message='', defaultDir='',
                 defaultFile='', wildcard='*.*', style=0, *args, **kw):
        self.message, self.wildcard, self.style = message, wildcard, style
        self.paths = []
    def ShowModal(self):
        ans = ANSWERS.take('file', self.message)
        if ans is None:
            return wx.ID_CANCEL
        self.paths = [os.path.abspath(p) for p in
                      (ans if isinstance(ans, (list, tuple)) else [ans])]
        return wx.ID_OK
    def GetPath(self): return self.paths[0] if self.paths else ''
    def GetPaths(self): return list(self.paths)
    def GetDirectory(self): return os.path.dirname(self.GetPath())
    def GetFilename(self): return os.path.basename(self.GetPath())
    def GetFilenames(self): return [os.path.basename(p) for p in self.paths]
    def GetWildcard(self): return self.wildcard
    def GetFilterIndex(self): return 0
    def SetFilterIndex(self, i): pass
    def SetDirectory(self, d): pass
    def SetFilename(self, f): pass
    def SetPath(self, p): pass
    def SetWildcard(self, w): self.wildcard = w
    def CenterOnParent(self, *a): pass
    def CentreOnParent(self, *a): pass
    def Raise(self): pass
    def Destroy(self): pass

class FakeMessageDialog(_Inert):
    'Stand-in for wx.MessageDialog; returns a button id queued in ANSWERS'
    def __init__(self, parent=None, message='', caption='', style=0,
                 *args, **kw):
        self.message, self.caption = message, caption
    def ShowModal(self):
        ans = ANSWERS.take('message', f'{self.caption}: {self.message}')
        return wx.ID_CANCEL if ans is None else ans
    def SetYesNoLabels(self, *a): return True
    def SetOKLabel(self, *a): return True
    def SetOKCancelLabels(self, *a): return True
    def SetYesNoCancelLabels(self, *a): return True
    def SetExtendedMessage(self, msg): self.message += '\n' + msg
    def CenterOnParent(self, *a): pass
    def CentreOnParent(self, *a): pass
    def Raise(self): pass
    def Destroy(self): pass

def fake_message_box(message, caption='Message', style=wx.OK, *a, **kw):
    'Stand-in for wx.MessageBox'
    ans = ANSWERS.take('message', f'{caption}: {message}')
    return wx.OK if ans is None else ans

def fake_file_browser(parent, message, defaultDir='', *args,
                      style=wx.FD_OPEN, **kw):
    '''Stand-in for GSASIIctrlGUI.G2FileBrowser, the file dialog GSAS-II
    uses instead of wx.FileDialog on Linux (config G2FileBrowser); answered
    from the same queue as wx.FileDialog'''
    dlg = FakeFileDialog(parent, message, defaultDir, style=style)
    if dlg.ShowModal() != wx.ID_OK:
        return []
    paths = dlg.GetPaths()
    return paths if style & wx.FD_MULTIPLE else paths[:1]

def install_fakes():
    '''Replace the native dialogs. GSAS-II looks these up on the wx module
    at call time, so patching wx is enough. Safe to call before the wx.App
    exists (no GSAS-II module is imported).
    '''
    if _REAL:
        return
    _REAL.update(FileDialog=wx.FileDialog, MessageDialog=wx.MessageDialog,
                 MessageBox=wx.MessageBox)
    wx.FileDialog = FakeFileDialog
    wx.MessageDialog = FakeMessageDialog
    wx.MessageBox = fake_message_box

def uninstall_fakes():
    'Restore the real native dialogs'
    for name, obj in _REAL.items():
        setattr(wx, name, obj)
    _REAL.clear()

_realSaveConfig = None
_realFileBrowser = None
GUI_EXCEPTIONS = []   # uncaught exceptions in GSAS-II event handlers

def install_excepthook():
    '''Log uncaught exceptions in GSAS-II event handlers and remember
    them, so that the scenario fails at its next wait instead of GSAS-II
    (in debug mode) stopping in an IPython console. Run after GSAS-II has
    started (it installs its own hook).'''
    def hook(etype, value, tb):
        text = ''.join(traceback.format_exception(etype, value, tb))
        GUI_EXCEPTIONS.append(f'{etype.__name__}: {value}')
        log('uncaught exception in GSAS-II:\n' + text)
    sys.excepthook = hook

def install_config_stub():
    '''Stop GSAS-II from writing ~/.GSASII/config.ini, and replace GSAS-II's
    own file browser (used on Linux) with :func:`fake_file_browser`. Must
    run on the GUI thread after GSAS-II has started: importing GSAS-II GUI
    modules before the wx.App exists breaks them (e.g. module-level wx
    constants such as WACV in GSASIIphsGUI fail silently, giving a
    NameError later).
    '''
    global _realSaveConfig, _realFileBrowser
    from . import GSASIIctrlGUI as G2G
    if _realSaveConfig is None:
        _realSaveConfig = G2G.SaveConfigVars
        _realFileBrowser = G2G.G2FileBrowser
    G2G.SaveConfigVars = lambda *a, **k: log('(config.ini save skipped)')
    G2G.G2FileBrowser = fake_file_browser

def uninstall_config_stub():
    from . import GSASIIctrlGUI as G2G
    if _realSaveConfig is not None:
        G2G.SaveConfigVars = _realSaveConfig
        G2G.G2FileBrowser = _realFileBrowser

#==============================================================================
# Demo-mode decorations: caption bar and widget outline
#==============================================================================
class _Caption(wx.Frame):
    'A borderless strip of large text shown below the GSAS-II window'
    def __init__(self):
        super().__init__(None, style=wx.FRAME_NO_TASKBAR | wx.STAY_ON_TOP
                         | wx.FRAME_TOOL_WINDOW | wx.BORDER_NONE)
        self.SetBackgroundColour(wx.Colour(255, 245, 170))
        self.text = wx.StaticText(self, style=wx.ALIGN_CENTRE_HORIZONTAL)
        font = self.text.GetFont()
        font.SetPointSize(font.GetPointSize() + 8)
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.text.SetFont(font)
        self.text.SetForegroundColour(wx.BLACK)
        # step mode: buttons that let the viewer continue or stop
        self.nextBtn = wx.Button(self, label='Next  ▸')
        self.stopBtn = wx.Button(self, label='Stop')
        bfont = self.nextBtn.GetFont()
        bfont.SetPointSize(bfont.GetPointSize() + 3)
        for b in (self.nextBtn, self.stopBtn):
            b.SetFont(bfont)
        self.nextBtn.Bind(wx.EVT_BUTTON, lambda e: self._answer('next'))
        self.stopBtn.Bind(wx.EVT_BUTTON, lambda e: self._answer('stop'))
        self.nextBtn.SetToolTip('Continue with the next step '
                                '(or press Enter in the terminal)')
        self.stopBtn.SetToolTip('End the demo here '
                                '(or type q and Enter in the terminal)')
        self.nextBtn.Hide(); self.stopBtn.Hide()
        self.on_choice = None    # called with 'next' or 'stop'

    def _answer(self, choice):
        'GUI thread: a step-mode button was pressed'
        if self.on_choice:
            self.on_choice(choice)

    def show(self, msg, frame, buttons=False):
        self.text.SetLabel(msg)
        self.nextBtn.Show(buttons); self.stopBtn.Show(buttons)
        r = frame.GetScreenRect()
        w = max(r.width, 400)
        bw = 0               # width taken by the buttons on the right
        if buttons:
            nw, nh = self.nextBtn.GetBestSize()
            sw, sh = self.stopBtn.GetBestSize()
            bw = nw + sw + 30
        self.text.Wrap(w - 40 - 2*bw)
        tw, th = self.text.GetBestSize()
        h = max(th + 24, 44)
        self.SetSize(r.x, r.y + r.height + 4, w, h)
        disp = wx.Display(max(0, wx.Display.GetFromWindow(frame)))
        area = disp.GetClientArea()
        if r.y + r.height + 4 + h > area.y + area.height:   # no room below
            self.SetPosition((r.x, r.y + r.height - h - 4))
        self.text.SetPosition(((w - tw) // 2, (h - th) // 2))
        if buttons:
            self.nextBtn.SetSize(w - nw - 10, (h - nh) // 2, nw, nh)
            self.stopBtn.SetSize(w - nw - sw - 20, (h - sh) // 2, sw, sh)
            self.nextBtn.SetDefault()
        if not self.IsShown():
            self.ShowWithoutActivating()
        self.Raise()
        self.Refresh(); self.Update()

class _Outline(wx.Frame):
    'A hollow coloured rectangle drawn around a widget (shaped frame)'
    def __init__(self, colour=wx.Colour(230, 30, 30), thick=4):
        super().__init__(None, style=wx.FRAME_SHAPED | wx.BORDER_NONE
                         | wx.FRAME_NO_TASKBAR | wx.STAY_ON_TOP
                         | wx.FRAME_TOOL_WINDOW)
        self.SetBackgroundColour(colour)
        self.thick = thick

    def show(self, rect):
        t = self.thick
        r = wx.Rect(rect.x - t - 2, rect.y - t - 2,
                    rect.width + 2*t + 4, rect.height + 2*t + 4)
        self.SetSize(r)
        region = wx.Region(0, 0, r.width, r.height)
        region.Subtract(wx.Rect(t, t, r.width - 2*t, r.height - 2*t))
        self.SetShape(region)
        if not self.IsShown():
            self.ShowWithoutActivating()
        self.Raise()
        self.Refresh(); self.Update()

#==============================================================================
# The driver
#==============================================================================
def _set_text(ctrl, text):
    '''GUI thread: put text in a text control as if it had been typed.
    A GSAS-II ValidatedTxtCtrl's own SetValue stores its argument unchanged
    (a str for a number entry), so for those the text is placed with the
    base-class method and checked by the validator, which stores the
    converted value, as happens when typing.'''
    if hasattr(ctrl, '_onLoseFocus') and hasattr(ctrl, 'result'):
        wx.TextCtrl.ChangeValue(ctrl, text)
        ctrl.SetModified(True)
        if getattr(ctrl, 'type', None) is str:
            ctrl.ShowStringValidity()
        else:
            ctrl._TestValidity()
    else:
        ctrl.SetValue(text)

class GUIDriver:
    '''Drives the GSAS-II GUI from a background (scenario) thread.

    All methods are called from the scenario thread; they send work to
    the GUI thread themselves. Code passed to :meth:`gui` runs on the GUI
    thread and may touch any wx object.

    :param G2frame: the GSAS-II main frame
    :param bool demo: True for demo mode (captions, outlines, pauses)
    :param float delay: demo-mode pause after each step, in seconds
    :param float typing: demo-mode delay per typed character, in seconds
    :param bool step_mode: demo mode in which each step waits for the
      viewer to press Next (implies demo=True)
    '''
    def __init__(self, G2frame, demo=False, delay=1.5, typing=0.08,
                 step_mode=False):
        self.frame = G2frame
        self.step_mode = step_mode
        self.demo = demo or step_mode
        self.delay = delay
        self.typing = typing
        self.answers = ANSWERS
        self.nsteps = 0
        self._step_text = ''
        self._caption_text = ''
        self._caption = None
        self._outline = None
        self._waiter = None          # (Event, dict) while waiting for Next
        self._wait_lock = threading.Lock()
        self._reader = None          # thread reading Enter from the terminal
        self._phase = None           # last phase added (default for phase tasks)
        self._tmpdir = None
        if self.demo:
            ANSWERS.on_show = self._show_native

    # ---- running code on the GUI thread ---------------------------------
    def gui(self, fn, *args, timeout=30, **kw):
        '''Run fn(*args, **kw) on the GUI thread, wait for it and return
        its result (an exception there is re-raised here). Not for actions
        that open a modal wx dialog; use :meth:`gui_nowait` for those.
        Called on the GUI thread itself (e.g. inside a wait_until
        condition), fn is simply called.
        '''
        if wx.IsMainThread():
            return fn(*args, **kw)
        done = threading.Event()
        box = {}
        def run():
            try:
                box['result'] = fn(*args, **kw)
            except BaseException as err:
                box['error'] = err
                box['tb'] = traceback.format_exc()
            finally:
                done.set()
        wx.CallAfter(run)
        if not done.wait(timeout * TIMEOUT_SCALE):
            raise TimeoutError(f'GUI call {getattr(fn, "__name__", fn)} '
                               f'did not finish in {timeout * TIMEOUT_SCALE:g}'
                               ' s (did it open a modal dialog? use '
                               'modal=True)')
        if 'error' in box:
            log('error on GUI thread:\n' + box['tb'])
            raise box['error']
        return box.get('result')

    def gui_nowait(self, fn, *args, **kw):
        'Start fn(*args, **kw) on the GUI thread and return at once'
        def run():
            try:
                fn(*args, **kw)
            except Exception:
                log('error on GUI thread:\n' + traceback.format_exc())
        wx.CallAfter(run)

    def wait_until(self, cond, what='condition', timeout=10):
        '''Call cond() on the GUI thread until it returns something truthy;
        return that value. Raises TimeoutError; raises ScenarioError at
        once if GSAS-II raises an exception or opens an unexpected native
        dialog meanwhile (e.g. an error message from a failed refinement).
        The timeout is multiplied by :data:`TIMEOUT_SCALE`.'''
        timeout *= TIMEOUT_SCALE
        end = time.time() + timeout
        nerr = len(self.answers.errors)
        while True:
            if GUI_EXCEPTIONS:
                raise ScenarioError('GSAS-II raised an exception while '
                                    f'waiting for {what}: {GUI_EXCEPTIONS[-1]}')
            if len(self.answers.errors) > nerr:
                raise ScenarioError(f'while waiting for {what}: '
                                    f'{self.answers.errors[-1]}')
            # GSAS-II may be busy (e.g. refining) for a while before the
            # GUI thread gets to the call, so allow at least 30 s for it
            val = self.gui(cond, timeout=max(30, timeout))
            if val:
                return val
            if time.time() > end:
                raise TimeoutError(f'timed out after {timeout:g} s waiting '
                                   f'for {what}')
            time.sleep(0.05)

    # ---- narration and pacing -------------------------------------------
    def step(self, text, wait=True):
        '''Start a step: logged in test mode; shown as a caption in demo
        mode, followed by a pause. In step mode the caption has Next and
        Stop buttons, and the step's actions wait until Next is pressed;
        with wait=False this step starts without waiting (as in demo mode),
        e.g. for a preparation step that has nothing to show.'''
        self.nsteps += 1
        log(f'step {self.nsteps}: {text}')
        self._step_text = text
        if self.step_mode and wait:
            self.wait_for_user(text)
        elif self.demo:
            self.caption(text)
            self.pause()

    def settle(self, seconds=0.25):
        '''Let delayed GUI updates run. GSAS-II often redraws a panel
        with ``wx.CallLater(100, ...)`` after a change (e.g. a check box
        click or a new histogram selection), so widgets found right after
        the change may be about to be replaced.'''
        time.sleep(seconds)
        self.gui(lambda: None)

    def caption(self, text):
        'Show text in the demo caption bar (does nothing in test mode)'
        if not self.demo:
            return
        self._caption_text = text
        def show():
            if self._caption is None:
                self._caption = _Caption()
            self._caption.show(text, self.frame)
        self.gui(show)

    def wait_for_user(self, text=None):
        '''Step mode: show text (default: the current caption) with Next
        and Stop buttons, and wait until the viewer presses one. Pressing
        Enter in the terminal GSAS-II was started from also counts as
        Next ("q" and Enter as Stop), which helps while a modal dialog
        blocks clicks on the caption bar. Stop raises
        :class:`ScenarioStopped`. Does nothing unless in step mode.'''
        if not self.step_mode:
            return
        text = self._caption_text if text is None else text
        self._caption_text = text
        event, box = threading.Event(), {}
        with self._wait_lock:
            self._waiter = (event, box)
        self._start_terminal_reader()
        def show(buttons):
            if self._caption is None:
                self._caption = _Caption()
            self._caption.on_choice = self._choose
            self._caption.show(text, self.frame, buttons=buttons)
            if buttons and not self._caption.IsEnabled():
                self._caption.Enable()   # a modal dialog disabled it
        self.gui(show, True)
        event.wait()
        with self._wait_lock:
            self._waiter = None
        self.gui(show, False)
        if box.get('choice') == 'stop':
            log('stopped by the viewer')
            raise ScenarioStopped('stopped by the viewer')

    def _choose(self, choice):
        'Any thread: the viewer chose "next" or "stop"'
        with self._wait_lock:
            if self._waiter:
                event, box = self._waiter
                box['choice'] = choice
                event.set()

    def _start_terminal_reader(self):
        '''Step mode: read lines typed in the terminal (if stdin is one);
        Enter means Next, "q" means Stop.'''
        if self._reader is not None:
            return
        self._reader = False
        try:
            if not sys.stdin.isatty():
                return
        except (AttributeError, ValueError):
            return
        def read():
            for line in sys.stdin:
                word = line.strip().lower()
                self._choose('stop' if word in ('q', 'quit', 'stop')
                             else 'next')
        self._reader = threading.Thread(target=read, daemon=True,
                                        name='GUIdriver-stdin')
        self._reader.start()
        log('step mode: press Next below the GSAS-II window, or Enter here '
            '(q + Enter to stop)')

    def detail(self, text):
        '''Demo mode: show text as a second caption line, below the
        current step's text'''
        self.caption(f'{self._step_text}\n{text}' if self._step_text else text)

    def pause(self, seconds=None):
        'Pause in demo mode (seconds defaults to the step delay)'
        if self.demo:
            time.sleep(self.delay if seconds is None else seconds)

    def highlight(self, target, seconds=None):
        '''Demo mode: outline a widget (or a wx.Rect in screen coordinates)
        for a moment. Does nothing in test mode.'''
        if not self.demo:
            return
        def show():
            rect = target if isinstance(target, wx.Rect) else \
                target.GetScreenRect()
            if self._outline is None:
                self._outline = _Outline()
            self._outline.show(rect)
        self.gui(show)
        self.pause(self.delay * 0.6 if seconds is None else seconds)
        self.unhighlight()

    def unhighlight(self):
        if self._outline is not None:
            self.gui(self._outline.Hide)

    def _show_native(self, kind, text, answer):
        '''Demo mode: native dialogs are answered by the stand-ins, so tell
        the viewer what would have been chosen. Runs on the GUI thread.'''
        text = text.split('\n')[0]
        if kind == 'file':
            msg = f'{text}:  {os.path.basename(str(answer))}'
        else:
            msg = f'{text}  →  {_answer_name(answer)}'
        if self._step_text:
            msg = f'{self._step_text}\n{msg}'
        if self._caption is None:
            self._caption = _Caption()
        self._caption.show(msg, self.frame)
        # pause on the GUI thread (we are inside the handler that opened
        # the dialog), yielding so that the caption is drawn
        end = time.time() + self.delay
        while time.time() < end:
            wx.GetApp().Yield(True)
            wx.MilliSleep(30)

    # ---- menus ----------------------------------------------------------
    def _find_menu_item(self, path):
        'GUI thread: (menu, item) for a menu path such as ("File","Open...")'
        if len(path) < 2:
            raise ValueError('menu path needs a menu title and an item label')
        mb = self.frame.GetMenuBar()
        if mb is None:
            raise LookupError('the GSAS-II window has no menu bar')
        pos = mb.FindMenu(path[0])
        if pos == wx.NOT_FOUND:
            titles = [wx.MenuItem.GetLabelText(mb.GetMenuLabel(i))
                      for i in range(mb.GetMenuCount())]
            raise LookupError(f'no menu {path[0]!r} in the menu bar {titles}')
        menu = mb.GetMenu(pos)
        for depth, label in enumerate(path[1:], start=1):
            for item in menu.GetMenuItems():
                if item.IsSeparator():
                    continue
                if wx.MenuItem.GetLabelText(item.GetItemLabel()) == label:
                    break
            else:
                labels = [wx.MenuItem.GetLabelText(i.GetItemLabel())
                          for i in menu.GetMenuItems() if not i.IsSeparator()]
                raise LookupError(f'no item {label!r} in menu '
                                  f'{" > ".join(path[:depth])}: {labels}')
            if depth < len(path) - 1:
                if not item.IsSubMenu():
                    raise LookupError(f'{label!r} is not a submenu')
                menu = item.GetSubMenu()
        if item.IsSubMenu():
            raise LookupError(f'{" > ".join(path)} is a submenu, not a command')
        return menu, item

    def _menu_command(self, path, check=None):
        'GUI thread: send the menu event, as a real selection would'
        menu, item = self._find_menu_item(path)
        if not item.IsEnabled():
            raise ScenarioError(f'menu item {" > ".join(path)} is disabled')
        evt = wx.CommandEvent(wx.wxEVT_MENU, item.GetId())
        if item.IsCheckable():
            new = (not item.IsChecked()) if check is None else bool(check)
            if item.GetKind() == wx.ITEM_RADIO:
                new = True
            item.Check(new)
            evt.SetInt(int(new))
        evt.SetEventObject(menu)
        # wx delivers menu events to the menu first, then to the window;
        # GSAS-II binds nearly all its menu commands on the main frame
        if not menu.ProcessEvent(evt):
            self.frame.GetEventHandler().ProcessEvent(evt)

    def menu(self, *path, modal=False, check=None, timeout=60):
        '''Run a menu command, e.g. ``d.menu('File', 'Open project...')`` or,
        for a submenu, ``d.menu('Import', 'Phase', 'from CIF file')``.
        Labels are matched without ``&`` and accelerators. GSAS-II changes
        the menu bar with the selected tree item or tab, so select those
        first.

        :param bool modal: True if the command opens a modal wx dialog
          (the call then returns at once; wait for the dialog with
          :meth:`dialog`)
        :param check: for a check item, the state to set (default: toggle)
        Raises LookupError if the item does not exist and ScenarioError if
        it is disabled.
        '''
        log('menu: ' + ' > '.join(path))
        if self.demo:
            self.detail('  ▸  '.join(path))
            self.pause(self.delay * 0.7)
        self.gui(self._find_menu_item, path)    # fail early, in this thread
        if modal:
            self.gui_nowait(self._menu_command, path, check)
        else:
            self.gui(self._menu_command, path, check, timeout=timeout)
        self.pause(self.delay * 0.3)

    def menu_items(self, title):
        'Labels of the items in a menu of the current menu bar'
        def get():
            mb = self.frame.GetMenuBar()
            pos = mb.FindMenu(title)
            if pos == wx.NOT_FOUND:
                return None
            return [wx.MenuItem.GetLabelText(i.GetItemLabel())
                    for i in mb.GetMenu(pos).GetMenuItems()
                    if not i.IsSeparator()]
        return self.gui(get)

    def menu_enabled(self, *path):
        'True if a menu item exists and is enabled'
        def get():
            try:
                return self._find_menu_item(path)[1].IsEnabled()
            except LookupError:
                return False
        return self.gui(get)

    # ---- dialogs --------------------------------------------------------
    def dialog(self, title, timeout=15, prefix=False):
        '''Wait for a shown modal dialog with this title and return it.
        With prefix=True the title need only start with the given text
        (e.g. ``d.dialog('Refinement results', prefix=True)``).'''
        def find():
            for w in wx.GetTopLevelWindows():
                if (isinstance(w, wx.Dialog) and w.IsShown() and w.IsModal()
                        and (w.GetTitle().startswith(title) if prefix
                             else w.GetTitle() == title)):
                    return w
        try:
            dlg = self.wait_until(find, f'dialog {title!r}', timeout)
        except TimeoutError:
            log('dialogs open:', self.open_dialogs())
            raise
        self.pause(self.delay * 0.5)
        return dlg

    def open_dialogs(self):
        '[(title, shown, modal), ...] for every dialog window'
        return self.gui(lambda: [(w.GetTitle(), w.IsShown(), w.IsModal())
                                 for w in wx.GetTopLevelWindows()
                                 if isinstance(w, wx.Dialog)])

    def wait_no_modal(self, timeout=15):
        'Wait until no modal dialog is showing'
        def none_open():
            return not any(isinstance(w, wx.Dialog) and w.IsShown()
                           and w.IsModal() for w in wx.GetTopLevelWindows())
        try:
            self.wait_until(none_open, 'modal dialogs to close', timeout)
        except TimeoutError:
            log('dialogs still open:', self.open_dialogs())
            raise

    def expect_file(self, message, path, optional=False):
        '''Queue the answer for the next native file dialog. message must
        be part of the dialog's message; path may be a list (multiple
        selection). Call before the action that opens the dialog. With
        optional=True the dialog may not appear (the answer is then
        skipped); remove an unused optional answer with
        ``d.answers.discard(entry)``. Returns the queue entry.'''
        return self.answers.add('file', message, path, optional)

    def expect_message(self, text, answer=wx.ID_OK, optional=False):
        '''Queue the answer (wx.ID_YES, wx.ID_NO, wx.ID_OK, ...) for the next
        native message dialog or wx.MessageBox; text must be part of its
        "caption: message". Call before the action that opens it. See
        :meth:`expect_file` for optional. Returns the queue entry.'''
        return self.answers.add('message', text, answer, optional)

    # ---- finding widgets ------------------------------------------------
    def find(self, cls=wx.Window, label=None, name=None, parent=None,
             timeout=10, wait=True):
        '''Find a shown widget below parent (default: the main frame) by
        class and/or label and/or window name. label may also be a
        function of the label, e.g. ``label=lambda s: s.startswith('Use')``.
        With wait=True, waits for it to appear (raises TimeoutError);
        otherwise returns None if absent.
        '''
        parent = parent or self.frame
        match = label if callable(label) else (
            lambda s: label is None or s == label)
        def walk(w):
            for c in w.GetChildren():
                if (isinstance(c, cls) and c.IsShownOnScreen()
                        and match(c.GetLabel())
                        and (name is None or c.GetName() == name)):
                    return c
                found = walk(c)
                if found:
                    return found
        what = f'{cls.__name__} label={label!r} name={name!r}'
        if wait:
            return self.wait_until(lambda: walk(parent), what, timeout)
        return self.gui(lambda: walk(parent))

    def find_after(self, label, cls=wx.TextCtrl, parent=None, timeout=10):
        '''Find the control that follows a label. GSAS-II usually creates
        a label (wx.StaticText or a labelled check box) and then its
        controls, so the first shown widget of class cls created after the
        widget whose label starts with ``label`` (and in the same parent
        window) is returned. E.g. ``d.find_after(' a = ')`` is the entry
        for the unit cell length a; ``d.find_after(' U (', wx.CheckBox)``
        is the refine box for instrument parameter U. label may also be a
        function of the label text, as for :meth:`find`.'''
        parent = parent or self.frame
        match = label if callable(label) else (
            lambda s: s.startswith(label))
        def look(w):
            kids = list(w.GetChildren())
            for i, c in enumerate(kids):
                if (c.IsShownOnScreen() and not isinstance(c, wx.TextCtrl)
                        and match(c.GetLabel())):
                    for n in kids[i+1:]:
                        if isinstance(n, cls) and n.IsShownOnScreen():
                            return n
                found = look(c)
                if found:
                    return found
        return self.wait_until(lambda: look(parent),
                               f'{cls.__name__} after {label!r}', timeout)

    def find_all(self, cls=wx.Window, parent=None):
        'All shown widgets of a class below parent (default: main frame)'
        parent = parent or self.frame
        def walk(w, out):
            for c in w.GetChildren():
                if isinstance(c, cls) and c.IsShownOnScreen():
                    out.append(c)
                walk(c, out)
            return out
        return self.gui(lambda: walk(parent, []))

    # ---- acting on widgets ----------------------------------------------
    @staticmethod
    def _click(widget, btn_id=None):
        'GUI thread: send a button event'
        if btn_id is not None:      # standard button of a dialog
            # NB dlg.FindWindowById() is static and searches ALL windows
            # (it can return a button of another, hidden, dialog)
            btn = widget.FindWindow(btn_id)
            target = btn or widget
            evt = wx.CommandEvent(wx.wxEVT_BUTTON, btn_id)
        else:
            target = widget
            evtType = wx.wxEVT_BUTTON
            if isinstance(widget, wx.CheckBox):
                widget.SetValue(not widget.GetValue())
                evtType = wx.wxEVT_CHECKBOX
            elif isinstance(widget, wx.ToggleButton):
                widget.SetValue(not widget.GetValue())
                evtType = wx.wxEVT_TOGGLEBUTTON
            elif isinstance(widget, wx.RadioButton):
                widget.SetValue(True)
                evtType = wx.wxEVT_RADIOBUTTON
            evt = wx.CommandEvent(evtType, widget.GetId())
            if evtType != wx.wxEVT_BUTTON:
                evt.SetInt(int(widget.GetValue()))
        evt.SetEventObject(target)
        target.GetEventHandler().ProcessEvent(evt)

    def click(self, widget, btn_id=None, modal=False, timeout=60):
        '''Click a button, check box, toggle or radio button. To press a
        standard dialog button pass the dialog and the id:
        ``d.click(dlg, wx.ID_OK)``.

        :param bool modal: True if the click opens a modal wx dialog
        '''
        if self.demo:
            target = widget
            if btn_id is not None:
                target = self.gui(widget.FindWindow, btn_id) or widget
            self.highlight(target)
        if modal:
            self.gui_nowait(self._click, widget, btn_id)
        else:
            self.gui(self._click, widget, btn_id, timeout=timeout)
        self.pause(self.delay * 0.3)

    def type(self, target, text, enter=False):
        '''Put text into a text control. target is the control, or a dialog
        (its first shown text control is used). In demo mode the text is
        typed one character at a time. With enter=True a text-enter event
        is sent afterwards (needed by GSAS-II validated entries that act on
        Enter).'''
        ctrl = target
        if isinstance(target, (wx.Dialog, wx.Frame, wx.Panel)) and \
                not isinstance(target, wx.TextCtrl):
            ctrl = self.find(wx.TextCtrl, parent=target)
        if self.demo:
            self.highlight(ctrl, self.delay * 0.3)
            self.gui(ctrl.SetFocus)
            for i in range(1, len(text) + 1):
                self.gui(_set_text, ctrl, text[:i])
                time.sleep(self.typing)
        self.gui(_set_text, ctrl, text)
        if enter:
            def send():
                evt = wx.CommandEvent(wx.wxEVT_TEXT_ENTER, ctrl.GetId())
                evt.SetEventObject(ctrl)
                evt.SetString(text)
                ctrl.GetEventHandler().ProcessEvent(evt)
            self.gui(send)
        self.pause(self.delay * 0.3)
        return ctrl

    def choose(self, ctrl, value, modal=False):
        'Select a string in a wx.Choice/wx.ComboBox/wx.ListBox and send the event'
        self.highlight(ctrl)
        def run():
            if not ctrl.SetStringSelection(value):
                raise LookupError(f'{value!r} not in {ctrl.GetStrings()}')
            if isinstance(ctrl, wx.ComboBox):
                kind = wx.wxEVT_COMBOBOX
            elif isinstance(ctrl, wx.ListBox):
                kind = wx.wxEVT_LISTBOX
            else:
                kind = wx.wxEVT_CHOICE
            evt = wx.CommandEvent(kind, ctrl.GetId())
            evt.SetEventObject(ctrl)
            evt.SetInt(ctrl.GetSelection())
            evt.SetString(value)
            ctrl.GetEventHandler().ProcessEvent(evt)
        (self.gui_nowait if modal else self.gui)(run)
        self.pause(self.delay * 0.3)

    # ---- grids ----------------------------------------------------------
    def grid_col(self, grid, label):
        'Index of the grid column with this label (raises LookupError)'
        def get():
            labels = [grid.GetColLabelValue(c)
                      for c in range(grid.GetNumberCols())]
            if label not in labels:
                raise LookupError(f'no column {label!r} in {labels}')
            return labels.index(label)
        return self.gui(get)

    def _cell_rect(self, grid, row, col):
        '''GUI thread: screen rectangle of a cell, or of a column label
        (row < 0) or row label (col < 0)'''
        if row < 0 and col < 0:
            return grid.GetGridCornerLabelWindow().GetScreenRect()
        r = grid.CellToRect(max(row, 0), max(col, 0))
        x, y = grid.CalcScrolledPosition(r.x, r.y)
        if row < 0:
            win = grid.GetGridColLabelWindow()
            return wx.Rect(win.ClientToScreen(wx.Point(x, 0)),
                           wx.Size(r.width, win.GetSize().height))
        if col < 0:
            win = grid.GetGridRowLabelWindow()
            return wx.Rect(win.ClientToScreen(wx.Point(0, y)),
                           wx.Size(win.GetSize().width, r.height))
        return wx.Rect(grid.GetGridWindow().ClientToScreen(wx.Point(x, y)),
                       r.GetSize())

    def edit_cell(self, grid, row, col, text):
        '''Type text into a grid cell with the cell's editor, as a user
        would; the grid then sends EVT_GRID_CELL_CHANGED. col may be a
        column label. GSAS-II often rebuilds a grid after a change, so
        find the grid again before each edit.'''
        if isinstance(col, str):
            col = self.grid_col(grid, col)
        log(f'cell ({row}, {self.gui(grid.GetColLabelValue, col)!r}) '
            f'= {text!r}')
        def _open_editor():
            'GUI thread: open the cell editor; return its control'
            grid.MakeCellVisible(row, col)
            grid.SetGridCursor(row, col)
            if not grid.CanEnableCellControl():
                raise ScenarioError(f'grid cell ({row}, {col}) is read-only')
            grid.EnableCellEditControl()
            editor = grid.GetCellEditor(row, col)
            ctrl = editor.GetControl()
            editor.DecRef()
            if not hasattr(ctrl, 'SetValue'):
                grid.DisableCellEditControl()
                raise ScenarioError(f'cannot type into the editor of cell '
                                    f'({row}, {col}): {ctrl!r}')
            return ctrl
        def _editing():
            'GUI thread: True if the editor is still open on this cell'
            return (grid.IsCellEditControlEnabled()
                    and grid.GetGridCursorRow() == row
                    and grid.GetGridCursorCol() == col)
        if self.demo:     # show the typing; the value is committed below
            self.highlight(self.gui(self._cell_rect, grid, row, col),
                           self.delay * 0.3)
            ctrl = self.gui(_open_editor)
            for i in range(1, len(text) + 1):
                self.gui(lambda: _editing() and ctrl.SetValue(text[:i]))
                time.sleep(self.typing)
        def _commit():
            '''GUI thread: put the text in the editor and close it, which
            stores the value and sends EVT_GRID_CELL_CHANGED. Done in one
            call, so no event can get in between. (wx closes a cell editor
            when it loses focus, e.g. when the window is activated or a
            dialog closes; text put in a closed editor would be lost.)'''
            reopened = False
            if not _editing():
                if grid.IsCellEditControlEnabled():   # open on another cell
                    grid.DisableCellEditControl()
                reopened = self.demo
                ctrl = _open_editor()
            else:
                editor = grid.GetCellEditor(row, col)
                ctrl = editor.GetControl()
                editor.DecRef()
            ctrl.SetValue(text)
            grid.DisableCellEditControl()
            return reopened
        if self.gui(_commit):
            log('(the cell editor had been closed, by a focus change?; '
                'it was reopened)')
        self.pause(self.delay * 0.2)

    def grid_event(self, grid, row, col,
                   kind=wx.grid.wxEVT_GRID_CELL_LEFT_DCLICK, modal=False):
        '''Send a grid mouse event, by default a double click on a cell.
        row=-1 means the column label, col=-1 the row label; col may be a
        column label. E.g. a double click on the "refine" column label:
        ``d.grid_event(grid, -1, 'refine',
        wx.grid.wxEVT_GRID_LABEL_LEFT_DCLICK, modal=True)``.

        :param bool modal: True if the handler opens a modal wx dialog
        '''
        if isinstance(col, str):
            col = self.grid_col(grid, col)
        log(f'grid event at ({row}, {col})')
        if self.demo:
            self.highlight(self.gui(self._cell_rect, grid, row, col))
        def send():
            evt = wx.grid.GridEvent(grid.GetId(), kind, grid, row, col)
            grid.GetEventHandler().ProcessEvent(evt)
        if modal:
            self.gui_nowait(send)
        else:
            self.gui(send)
        self.pause(self.delay * 0.3)

    # ---- data tree and tabs ---------------------------------------------
    def _tree_item(self, path):
        'GUI thread: data tree item for a path of item labels, or None'
        tree = self.frame.GPXtree
        item = tree.GetRootItem()
        for label in path:
            child, cookie = tree.GetFirstChild(item)
            while child and child.IsOk():
                if tree.GetItemText(child) == label:
                    break
                child, cookie = tree.GetNextChild(item, cookie)
            else:
                return None
            if not (child and child.IsOk()):
                return None
            item = child
        return item

    def tree_children(self, *path):
        'Labels of the children of a data tree item (default: top level)'
        def get():
            tree = self.frame.GPXtree
            item = self._tree_item(path)
            if item is None:
                return None
            out = []
            child, cookie = tree.GetFirstChild(item)
            while child and child.IsOk():
                out.append(tree.GetItemText(child))
                child, cookie = tree.GetNextChild(item, cookie)
            return out
        return self.gui(get)

    def tree_data(self, *path):
        '''The data attached to a data tree item, e.g.
        ``d.tree_data('Phases', 'Test phase')``, or None if absent. This is
        the live object; read it inside :meth:`gui` if GSAS-II may be
        changing it.'''
        def get():
            item = self._tree_item(path)
            return None if item is None else \
                self.frame.GPXtree.GetItemPyData(item)
        return self.gui(get)

    def select_tree(self, *path, timeout=10):
        'Select a data tree item, e.g. ``d.select_tree("Phases", "Test phase")``'
        item = self.wait_until(lambda: self._tree_item(path),
                               'tree item ' + ' > '.join(path), timeout)
        log('tree: ' + ' > '.join(path))
        if self.demo:
            def rect():
                tree = self.frame.GPXtree
                tree.EnsureVisible(item)
                r = tree.GetBoundingRect(item, textOnly=True)
                if r is None:
                    return tree.GetScreenRect()
                return wx.Rect(tree.ClientToScreen(r.GetTopLeft()), r.GetSize())
            self.highlight(self.gui(rect))
        self.gui(self.frame.GPXtree.SelectItem, item)
        self.pause(self.delay * 0.3)

    def select_tab(self, name, notebook=None, timeout=10):
        '''Select a tab by its label in a wx notebook (default: the phase
        notebook, G2frame.phaseDisplay). Phase tabs can also be selected
        with the menu, ``d.menu('Select tab', name)``.'''
        def find():
            nb = notebook or getattr(self.frame, 'phaseDisplay', None)
            if not nb:
                return None
            for i in range(nb.GetPageCount()):
                if nb.GetPageText(i) == name:
                    return nb, i
        nb, i = self.wait_until(find, f'tab {name!r}', timeout)
        log(f'tab: {name}')
        self.gui(nb.SetSelection, i)
        self.pause(self.delay * 0.3)

    def current_tab(self, notebook=None):
        'Label of the selected tab (default: phase notebook)'
        def get():
            nb = notebook or getattr(self.frame, 'phaseDisplay', None)
            if not nb or nb.GetSelection() < 0:
                return None
            return nb.GetPageText(nb.GetSelection())
        return self.gui(get)

    def has_menu(self, title):
        'True if the current menu bar has a menu with this title'
        return self.gui(lambda: self.frame.GetMenuBar().FindMenu(title)
                        != wx.NOT_FOUND)

    # ---- GSAS-II tasks --------------------------------------------------
    # Common GUI operations, done as a user would (menus, dialogs, tabs).
    # Phase tasks take an optional phase name; the default is the phase
    # last added with add_phase, or the only phase in the project.
    @property
    def tmpdir(self):
        '''A scratch directory for downloads and saved projects, made on
        first use (it is not removed afterwards)'''
        if self._tmpdir is None:
            self._tmpdir = tempfile.mkdtemp(prefix='G2GUIscenario_')
            log('working directory', self._tmpdir)
        return self._tmpdir

    def download(self, urls, dest=None):
        '''Download one URL (returns the local path) or a list of URLs
        (returns a list of paths) into dest (default :attr:`tmpdir`); each
        file keeps the last part of its URL as its name. The GUI's file
        dialogs need local files, so download input files first and give
        their paths to :meth:`expect_file` or :meth:`import_powder`.'''
        import urllib.request
        import urllib.parse
        single = isinstance(urls, str)
        dest = dest or self.tmpdir
        paths = []
        for url in ([urls] if single else urls):
            name = os.path.basename(urllib.parse.urlparse(url).path)
            path = os.path.join(dest, name)
            try:
                urllib.request.urlretrieve(url, path)
            except Exception:     # e.g. no SSL certificates for urllib
                import requests
                r = requests.get(url, timeout=60)
                r.raise_for_status()
                with open(path, 'wb') as fp:
                    fp.write(r.content)
            log(f'downloaded {url} ({os.path.getsize(path)} bytes)')
            paths.append(path)
        return paths[0] if single else paths

    def get_files(self, names, url, local_dirs=()):
        '''Get input files into :attr:`tmpdir` and return a dict of
        name: local path. Each file is copied from the first directory in
        local_dirs that has it (e.g. the tests/testinp directory of a
        GSAS-II source tree, so a self-test does not need the network),
        otherwise downloaded from url + name. Set the environment variable
        G2_GUI_DOWNLOAD=1 to always download.'''
        import shutil
        if os.environ.get('G2_GUI_DOWNLOAD', '') not in ('', '0'):
            local_dirs = ()
        local_dirs = [os.path.abspath(p) for p in local_dirs if p]
        paths = {}
        for name in names:
            for d in local_dirs:
                src = os.path.join(d, name)
                if os.path.isfile(src):
                    paths[name] = shutil.copy(src, self.tmpdir)
                    log(f'copied {name} from {d}')
                    break
            else:
                paths[name] = self.download(url + name)
        return paths

    def import_powder(self, data, instprm=None,
                      reader='guess format from file', timeout=60):
        '''Import a powder data file with Import > Powder Data > reader
        and return the names of the new histograms (a list). The file
        preview question is answered Yes; the instrument parameter file
        dialog (shown when no matching file is found next to the data
        file) gets instprm. Both are optional: GSAS-II may not ask.

        :param str data: local path of the data file
        :param str instprm: local path of the instrument parameter file
        :param str reader: the reader's item in the Powder Data submenu,
          e.g. 'from GSAS powder data file'
        '''
        before = set(self.tree_children() or [])
        self.expect_file('Choose Powder Data input file', [data])
        optional = [self.expect_message('Is this the file you want?',
                                        wx.ID_YES, optional=True)]
        if instprm:
            optional.append(self.expect_file('Choose inst. param file',
                                             instprm, optional=True))
        try:
            self.menu('Import', 'Powder Data', reader, timeout=timeout)
            new = self.wait_until(
                lambda: [n for n in self.tree_children() if n.startswith(
                    'PWDR') and n not in before],
                f'histogram from {os.path.basename(data)}', timeout=timeout)
        finally:
            for entry in optional:
                self.answers.discard(entry)
        log(f'imported {os.path.basename(data)}: {new}')
        return new

    def phase_names(self):
        'Names of the phases in the data tree'
        return self.tree_children('Phases') or []

    def _phase_name(self, name=None):
        'The phase to work on: name, else the last added, else the only one'
        if name:
            return name
        names = self.phase_names()
        if self._phase in names:
            return self._phase
        if len(names) == 1:
            return names[0]
        raise ScenarioError(f'which phase? (phases: {names})')

    def phase_data(self, phase=None):
        '''The data dict of a phase (live object; see :meth:`tree_data`)'''
        return self.tree_data('Phases', self._phase_name(phase))

    def add_phase(self, name):
        '''Data > Add new phase, named name. Its General tab is then
        shown, and it becomes the default phase for the phase tasks.'''
        self.menu('Data', 'Add new phase', modal=True)
        dlg = self.dialog('Phase Name Entry')
        self.type(dlg, name)
        self.click(dlg, wx.ID_OK)
        self.wait_until(lambda: self._tree_item(('Phases', name)),
                        f'phase {name!r} in the data tree')
        self._phase = name
        self.wait_until(lambda: self.has_menu('Select tab'), 'phase menus')
        return name

    def _selected_phase(self):
        'GUI thread: name of the phase selected in the data tree, or None'
        tree = self.frame.GPXtree
        item = tree.GetSelection()
        if not (item and item.IsOk()):
            return None
        parent = tree.GetItemParent(item)
        if parent and parent.IsOk() and tree.GetItemText(parent) == 'Phases':
            return tree.GetItemText(item)
        return None

    def phase_tab(self, tab, phase=None):
        '''Show a tab of a phase ('General', 'Atoms', 'Data', ...): selects
        the phase in the data tree if needed, then uses the Select tab
        menu unless the tab is already shown'''
        name = self._phase_name(phase)
        if self.gui(self._selected_phase) != name:
            self.select_tree('Phases', name)
        # the phase window (and its menus) is built after the selection
        self.wait_until(lambda: tab in (self.menu_items('Select tab') or []),
                        f'phase menus ({tab} in "Select tab")')
        if self.current_tab() != tab:
            self.menu('Select tab', tab)
            self.wait_until(lambda: self.current_tab() == tab, f'{tab} tab')

    def set_space_group(self, sg, phase=None):
        '''General tab: click the space group button and enter sg, e.g.
        'P n m a' (spaces between the axis symbols)'''
        name = self._phase_name(phase)
        self.phase_tab('General', name)
        self.click(self.find_after(lambda s: s.strip() == 'Space group:',
                                   wx.Button), modal=True)
        dlg = self.dialog('Get Space Group')
        self.type(dlg, sg)
        self.click(dlg, wx.ID_OK)
        dlg = self.dialog('Space Group Information')
        self.click(dlg, wx.ID_OK)
        self.wait_no_modal()
        norm = lambda s: s.replace(' ', '').lower()
        try:
            self.wait_until(lambda: norm(self.phase_data(name)['General']
                                         ['SGData']['SpGrp']) == norm(sg),
                            f'space group {sg}', timeout=5)
        except TimeoutError:
            got = self.gui(lambda: self.phase_data(name)['General']
                           ['SGData']['SpGrp'])
            raise ScenarioError(f'space group is {got!r}, not {sg!r}')

    _CELL = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')

    def set_cell(self, phase=None, **cell):
        '''General tab: type unit cell values, e.g.
        ``d.set_cell(a=8.48, b=5.398, c=6.958)``. Only the entries shown
        for the lattice can be set (e.g. only a for a cubic cell).'''
        name = self._phase_name(phase)
        self.phase_tab('General', name)
        for key, val in cell.items():
            if key not in self._CELL:
                raise ValueError(f'unknown cell parameter {key!r}; '
                                 f'use {self._CELL}')
            # ' a = ', or ' Unit cell: a = ' for a cubic cell
            label = lambda s, k=key: s.strip() in (f'{k} =',
                                                   f'Unit cell: {k} =')
            self.type(self.find_after(label), str(val), enter=True)
        got = self.gui(lambda: list(self.phase_data(name)['General']['Cell']))
        for key, val in cell.items():
            i = self._CELL.index(key) + 1
            if abs(got[i] - float(val)) > 1e-5:
                raise ScenarioError(f'cell {key} is {got[i]}, not {val}')

    def atoms_grid(self, phase=None):
        '''The phase's Atoms grid (shows the Atoms tab). GSAS-II rebuilds
        the grid after changes, so get it again before each use.'''
        self.phase_tab('Atoms', phase)
        return self.find(wx.grid.Grid)

    def append_atoms(self, n=1, phase=None):
        '''Atoms tab: Edit Atoms > Append atom, n times; returns the new
        rows'''
        name = self._phase_name(phase)
        self.phase_tab('Atoms', name)
        start = len(self.phase_data(name)['Atoms'])
        for i in range(n):
            self.menu('Edit Atoms', 'Append atom')
            self.wait_until(
                lambda: len(self.phase_data(name)['Atoms']) >= start + i + 1,
                f'atom {start + i + 1} added')
        return list(range(start, start + n))

    def _atom(self, name, row):
        'GUI thread: (atom record, AtomPtrs) for an atom'
        data = self.phase_data(name)
        return data['Atoms'][row], data['General']['AtomPtrs']

    def set_atom_type(self, row, elem, phase=None):
        '''Atoms tab: double-click an atom's Type cell and pick elem (e.g.
        'Pb' or 'Pb+2') in the periodic table'''
        name = self._phase_name(phase)
        self.grid_event(self.atoms_grid(name), row, 'Type', modal=True)
        dlg = self.dialog('Pick Element')
        # one combo box per element, listing it and its ions
        combo = self.wait_until(
            lambda: [c for c in self.find_all(wx.ComboBox, parent=dlg)
                     if elem in c.GetStrings()],
            f'element {elem} in the periodic table')[0]
        self.choose(combo, elem)
        self.wait_no_modal()
        def done():
            atom, ptrs = self._atom(name, row)
            return atom[ptrs[1]] == elem
        self.wait_until(done, f'atom {row + 1} type {elem}')

    def set_atom(self, row, elem=None, xyz=None, cols=None, phase=None):
        '''Atoms tab: set an atom's type (see :meth:`set_atom_type`), its
        coordinates (three values, typed into the x, y, z cells) and any
        other cells (cols: {column label: value}, e.g. {'Uiso': '0.02'}).
        Values are typed as text.'''
        name = self._phase_name(phase)
        if elem:
            self.set_atom_type(row, elem, name)
        cells = dict(zip('xyz', xyz)) if xyz else {}
        cells.update(cols or {})
        for col, val in cells.items():
            self.edit_cell(self.atoms_grid(name), row, col, str(val))
        if xyz:
            def number(v):
                try:
                    return float(v)
                except (TypeError, ValueError):
                    return None     # e.g. '1/3'; not checked
            want = [number(v) for v in xyz]
            def done():
                atom, ptrs = self._atom(name, row)
                got = atom[ptrs[0]:ptrs[0] + 3]
                return all(w is None or abs(g - w) < 1e-5
                           for g, w in zip(got, want))
            self.wait_until(done, f'atom {row + 1} coordinates')

    def set_check(self, target, value=True, parent=None):
        '''Set a check box (or a label prefix: the first shown check box
        whose label starts with it) to value, clicking it only if needed;
        returns the check box'''
        for attempt in range(5):
            box = target if isinstance(target, wx.CheckBox) else self.find(
                wx.CheckBox, label=lambda s: s.startswith(target),
                parent=parent)
            try:
                if self.gui(box.GetValue) != bool(value):
                    self.click(box)
                    self.settle()
                return box
            except RuntimeError as err:
                # the panel was redrawn under us: find the box again
                if isinstance(target, wx.CheckBox) or 'deleted' not in str(err):
                    raise
                log(f'check box {target!r} was replaced; trying again')
                self.settle()
        raise ScenarioError(f'check box {target!r} keeps being replaced')

    def hap_location(self):
        '''Where the phase/histogram (HAP) settings are shown: 'Hist/Phase
        item' (GSAS-II configured with SeparateHistPhaseTreeItem) or
        'phase Data tab' (for captions)'''
        return ('Hist/Phase item' if 'Hist/Phase' in self.tree_children()
                else 'phase Data tab')

    def show_hap(self, phase=None):
        '''Show the phase/histogram (HAP) settings of a phase: the phase's
        Data tab, or the phase's tab of the Hist/Phase data tree item
        (selected with its Select tab menu) when GSAS-II is configured with
        SeparateHistPhaseTreeItem'''
        name = self._phase_name(phase)
        if 'Hist/Phase' in self.tree_children():
            self.select_tree('Hist/Phase')
            self.wait_until(lambda: name in (self.menu_items('Select tab')
                                             or []), 'Hist/Phase menus')
            if self.gui(lambda: getattr(self.frame, 'HistPhaseLastSel',
                                        None)) != name:
                self.menu('Select tab', name)
                self.wait_until(lambda: self.frame.HistPhaseLastSel == name,
                                f'Hist/Phase tab for {name}')
        else:
            self.phase_tab('Data', name)
        self.wait_until(lambda: self.has_menu('Edit Phase'), 'Edit Phase menu')

    def select_hap_histogram(self, hist):
        '''HAP settings (see :meth:`show_hap`): select a histogram in the
        list at the top'''
        lb = self.wait_until(
            lambda: [c for c in self.find_all(wx.ListBox)
                     if hist in c.GetStrings()],
            'histogram list of the HAP settings')[0]
        if self.gui(lb.GetStringSelection) != hist:
            self.choose(lb, hist)
        self.wait_until(lambda: self.frame.hist == hist, f'{hist} selected')
        self.settle()     # the settings are redrawn 0.1 s later

    def add_histograms(self, hists=None, phase=None):
        '''HAP settings: Edit Phase > Add powder histograms, choosing
        hists (default: all, with Set All). Returns the histograms added.'''
        name = self._phase_name(phase)
        self.show_hap(name)
        used = set(self.gui(lambda: list(self.phase_data(name)['Histograms'])))
        # the dialog lists the PWDR histograms not in the phase, in tree order
        avail = [n for n in self.tree_children()
                 if n.startswith('PWDR') and n not in used]
        want = avail if hists is None else list(hists)
        missing = [h for h in want if h not in avail]
        if missing or not want:
            raise ScenarioError(f'cannot add {missing or "any histogram"} to '
                                f'{name}: available {avail}')
        self.menu('Edit Phase', 'Add powder histograms', modal=True)
        dlg = self.dialog('Use data')
        if hists is None:
            self.click(self.find(wx.Button, label='Set All', parent=dlg))
        else:
            self.highlight(self.find(wx.CheckListBox, parent=dlg))
            self.gui(dlg.SetSelections, [avail.index(h) for h in want])
            self.pause(self.delay * 0.5)
        self.click(dlg, wx.ID_OK)
        self.wait_no_modal()
        self.wait_until(
            lambda: all(h in self.phase_data(name)['Histograms'] for h in want),
            f'histograms added to {name}')
        return want

    def refine(self, gpx=None, timeout=600):
        '''Calculate > Refine, then accept the results (OK in the
        "Refinement results" dialog), and wait for GSAS-II to reload the
        project. A project that was never saved needs gpx, the file to
        save it as. Returns a dict: 'title' (of the results dialog), 'Rwp'
        (total, %), 'GOF', 'wR' ({histogram: wR %}) and 'Rvals' (all of
        the Covariance item's Rvals).'''
        if gpx:
            self.expect_file('Choose GSAS-II project file name', gpx)
        elif not self.gui(lambda: self.frame.GSASprojectfile):
            raise ScenarioError('the project has never been saved; '
                                'pass gpx (a file name) to refine()')
        phases = self.phase_names()
        old = self.tree_data('Covariance')
        self.menu('Calculate', 'Refine', modal=True)
        dlg = self.dialog('Refinement results', prefix=True, timeout=timeout)
        title = self.gui(dlg.GetTitle)
        self.detail(title)
        self.pause()
        self.click(dlg, wx.ID_OK)
        # accepting reloads the project from the .gpx file; wait for that
        def reloaded():
            cov = self.tree_data('Covariance')
            return (cov is not old and 'Rvals' in (cov or {})
                    and all(self._tree_item(('Phases', p)) for p in phases))
        self.wait_until(reloaded, 'project reloaded after refinement',
                        timeout=120)
        self.wait_no_modal(timeout=60)
        def results():
            rv = dict(self.tree_data('Covariance')['Rvals'])
            wR = {}
            for h in self.tree_children():
                if h.startswith(('PWDR', 'HKLF')):
                    try:
                        wR[h] = self.tree_data(h)[0]['wR']
                    except (KeyError, IndexError, TypeError):
                        pass
            return {'title': title, 'Rwp': rv.get('Rwp'),
                    'GOF': rv.get('GOF'), 'wR': wR, 'Rvals': rv}
        res = self.gui(results)
        log(f'refinement: Rwp {res["Rwp"]}, GOF {res["GOF"]}, ' +
            ', '.join(f'{h} {v:.2f}' for h, v in res['wR'].items()))
        return res

    # ---- checks and screenshots -----------------------------------------
    def check(self, cond, msg):
        '''Check that cond is true (if callable, it is called on the GUI
        thread first). Raises ScenarioError if not.'''
        val = self.gui(cond) if callable(cond) else cond
        if not val:
            log('FAILED:', msg)
            if self.demo:
                self.caption('Check failed: ' + msg)
            raise ScenarioError(msg)
        log('ok:', msg)

    def check_native_dialogs(self):
        'Check that every queued native-dialog answer was used, and no other'
        self.check(not self.answers.errors,
                   'no unexpected native dialogs ' + repr(self.answers.errors))
        self.check(not self.answers.queue,
                   'all queued native-dialog answers were used '
                   + repr(self.answers.queue))
        self.check(not GUI_EXCEPTIONS,
                   'no uncaught exceptions in GSAS-II ' + repr(GUI_EXCEPTIONS))

    def screenshot(self, path, window=None):
        '''Save a PNG of a window (default: the main frame), drawn from the
        window's own contents. Native dialogs and decorations are not
        included. Returns the path.'''
        def grab():
            win = window or self.frame
            win.Update()
            w, h = win.GetClientSize()
            bmp = wx.Bitmap(w, h)
            mdc = wx.MemoryDC(bmp)
            mdc.Blit(0, 0, w, h, wx.ClientDC(win), 0, 0)
            mdc.SelectObject(wx.NullBitmap)
            bmp.SaveFile(path, wx.BITMAP_TYPE_PNG)
            return path
        return self.gui(grab)

    def save_failure(self, directory, name='scenario', err=None):
        '''Save what is on screen after a failure into directory:
        name-frame.png (the main frame), name-windowN.png (other shown
        top-level windows, e.g. an open dialog), name-screen.png (the whole
        screen, which also shows native windows; may be blank where the
        screen cannot be read) and name-failure.txt (the error, the step
        and the open windows). Each part is skipped if it fails.'''
        os.makedirs(directory, exist_ok=True)
        base = os.path.join(directory, name)
        saved = []
        def windows():
            return [(w, w.GetTitle(), type(w).__name__)
                    for w in wx.GetTopLevelWindows() if w.IsShown()]
        try:
            wins = self.gui(windows, timeout=10)
        except Exception:
            wins = []
        def tryit(fn, *args):
            try:
                saved.append(fn(*args))
            except Exception as e:
                log(f'could not save {args[0] if args else fn}: {e!r}')
        tryit(self.screenshot, base + '-frame.png')
        n = 0
        for w, title, cls in wins:
            if w is self.frame or w in (self._caption, self._outline):
                continue
            n += 1
            tryit(self.screenshot, f'{base}-window{n}.png', w)
        def screen(path):
            def grab():
                w, h = wx.GetDisplaySize()
                bmp = wx.Bitmap(w, h)
                mdc = wx.MemoryDC(bmp)
                mdc.Blit(0, 0, w, h, wx.ScreenDC(), 0, 0)
                mdc.SelectObject(wx.NullBitmap)
                bmp.SaveFile(path, wx.BITMAP_TYPE_PNG)
                return path
            return self.gui(grab, timeout=10)
        tryit(screen, base + '-screen.png')
        def report(path):
            with open(path, 'w') as fp:
                fp.write(f'error: {err!r}\n')
                fp.write(f'step {self.nsteps}: {self._step_text}\n')
                fp.write('shown top-level windows:\n')
                for w, title, cls in wins:
                    fp.write(f'  {cls}: {title!r}\n')
                fp.write(f'unused native-dialog answers: {self.answers.queue!r}\n')
                fp.write(f'native-dialog errors: {self.answers.errors!r}\n')
                fp.write(f'GSAS-II exceptions: {GUI_EXCEPTIONS!r}\n')
            return path
        tryit(report, base + '-failure.txt')
        log('failure details saved:', ', '.join(saved))
        return saved

    def close_decorations(self):
        'Remove the demo caption and outline windows'
        def run():
            for w in (self._caption, self._outline):
                if w:
                    w.Destroy()
            self._caption = self._outline = None
        self.gui(run)

#==============================================================================
# Runner
#==============================================================================
def _wait_for_main_frame(timeout=180):
    end = time.time() + timeout * TIMEOUT_SCALE
    while wx.GetApp() is None:
        if time.time() > end:
            raise TimeoutError('wx.App was never created')
        time.sleep(0.1)
    while True:
        box = {}
        done = threading.Event()
        def look():
            top = wx.GetApp().GetTopWindow()
            if top is not None and top.IsShown() and hasattr(top, 'GPXtree'):
                box['frame'] = top
            done.set()
        wx.CallAfter(look)
        done.wait(10)
        if 'frame' in box:
            return box['frame']
        if time.time() > end:
            raise TimeoutError('GSAS-II main frame never appeared')
        time.sleep(0.2)

def run(scenario, demo=False, delay=1.5, keep_open=False, startup_wait=1.0,
        step=False, artifacts=None):
    '''Start GSAS-II and run scenario(driver) in a background thread.

    :param scenario: function taking a :class:`GUIDriver`
    :param bool demo: run in demo mode
    :param float delay: demo-mode pause per step (s)
    :param bool step: step mode: demo mode where each step waits for the
      viewer to press Next (Stop ends the run, with exit status 0)
    :param bool keep_open: leave GSAS-II running afterwards, with the real
      dialogs and config saving restored (default: exit; the exit status is
      0 if the scenario passed, 1 if not)
    :param float startup_wait: time (s) to let GSAS-II finish starting
    :param str artifacts: directory for screenshots and a report if the
      scenario fails (default: the G2_GUI_ARTIFACTS environment variable;
      none are saved if neither is set); see :meth:`GUIDriver.save_failure`

    Does not return unless keep_open is set and the user closes GSAS-II.
    '''
    sys.argv = sys.argv[:1]     # GSAS-II imports a .py argument in debug mode
    demo = demo or step
    install_fakes()
    result = {'status': 1}
    mode = '(step mode)' if step else '(demo mode)' if demo else ''
    artifacts = artifacts or os.environ.get('G2_GUI_ARTIFACTS')
    name = getattr(scenario, '__module__', '') or 'scenario'
    name = name.replace('G2scenario_', '')

    def runner():
        status = 1
        d = None
        try:
            G2frame = _wait_for_main_frame()
            d = GUIDriver(G2frame, demo=demo, delay=delay, step_mode=step)
            d.gui(install_config_stub)
            if not keep_open:
                d.gui(install_excepthook)
            time.sleep(startup_wait)    # let startup CallAfters finish
            log('GSAS-II is up; running scenario',
                getattr(scenario, '__module__', ''), mode)
            t0 = time.time()
            scenario(d)
            d.check_native_dialogs()
            log(f'PASS ({d.nsteps} steps, {time.time()-t0:.1f} s)')
            status = 0
            if step:
                d.wait_for_user('Done. Press Next or Stop to finish')
            elif demo:
                d.caption('Done')
                d.pause(2 * delay)
        except ScenarioStopped:
            if status != 0:
                log(f'STOPPED by the viewer at step {d.nsteps}')
            status = 0
        except Exception as err:
            log('FAIL:', repr(err))
            traceback.print_exc()
            if d is not None and artifacts:
                d.save_failure(artifacts, name, err)
            if d is not None and step:
                try:     # leave the failure on screen until the viewer is done
                    d.wait_for_user(None if isinstance(err, ScenarioError)
                                    else f'Failed: {err}')
                except ScenarioStopped:
                    pass
            elif d is not None and demo:
                time.sleep(3 * delay)
        result['status'] = status
        sys.stdout.flush(); sys.stderr.flush()
        if keep_open:
            ANSWERS.on_show = None
            if d is not None:
                d.close_decorations()
            def restore():
                uninstall_fakes()
                uninstall_config_stub()
            wx.CallAfter(restore)
            log('scenario finished; GSAS-II left running')
        else:
            os._exit(status)

    threading.Thread(target=runner, daemon=True, name='GUIdriver').start()
    from . import GSASIIGUI
    GSASIIGUI.main()
    return result['status']

# the demos directory at the top of a GSAS-II source tree (git checkout);
# it is not installed with the package
SCENARIO_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'demos')

def find_scenario(path):
    '''Return the path of a scenario file: path itself if it exists,
    otherwise the file of that name (.py optional) in the demos
    directory of the source tree'''
    if os.path.exists(path):
        return os.path.abspath(path)
    for name in (path, path + '.py'):
        p = os.path.join(SCENARIO_DIR, name)
        if os.path.isfile(p):
            return p
    raise SystemExit(f'scenario {path} not found (also looked in '
                     f'{SCENARIO_DIR})')

def list_scenarios():
    'Paths of the scenario files in the demos directory'
    import glob
    return sorted(p for p in glob.glob(os.path.join(SCENARIO_DIR, '*.py'))
                  if not os.path.basename(p).startswith('_'))

def load_scenario(path):
    'Load a scenario file and return its scenario() function'
    import importlib.util
    path = find_scenario(path)
    spec = importlib.util.spec_from_file_location(
        'G2scenario_' + os.path.splitext(os.path.basename(path))[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not callable(getattr(mod, 'scenario', None)):
        raise SystemExit(f'{path} does not define scenario(d)')
    return mod.scenario

def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(
        prog='python -m GSASII.GSASIIGUIdriver',
        description='Run a GSAS-II GUI scenario as a self-test or a demo')
    p.add_argument('scenario', nargs='?',
                   help='Python file defining scenario(d), or the name of '
                   'one in demos/ (e.g. CWCombined_tutorial)')
    p.add_argument('--list', action='store_true',
                   help='list the scenarios in demos/')
    p.add_argument('--demo', action='store_true',
                   help='demo mode: captions, outlines, typing and pauses')
    p.add_argument('--delay', type=float, default=1.5,
                   help='demo-mode pause per step in seconds (default 1.5)')
    p.add_argument('--step', action='store_true',
                   help='step mode (implies --demo): each step waits until '
                   'the viewer presses Next (or Enter in the terminal)')
    p.add_argument('--keep-open', action='store_true',
                   help='leave GSAS-II running when the scenario ends')
    p.add_argument('--timeout-scale', type=float, metavar='X',
                   help='multiply all timeouts by X, e.g. 3 on a slow '
                   'machine (default: $G2_GUI_TIMEOUT_SCALE or 1)')
    p.add_argument('--artifacts', metavar='DIR',
                   help='if the scenario fails, save screenshots and a '
                   'report in DIR (default: $G2_GUI_ARTIFACTS)')
    args = p.parse_args(argv)
    if args.list:
        for path in list_scenarios():
            print(os.path.splitext(os.path.basename(path))[0])
        return 0
    if not args.scenario:
        p.error('a scenario is required')
    if args.timeout_scale:
        global TIMEOUT_SCALE
        TIMEOUT_SCALE = args.timeout_scale
    scenario = load_scenario(args.scenario)
    return run(scenario, demo=args.demo, delay=args.delay,
               keep_open=args.keep_open, step=args.step,
               artifacts=args.artifacts)

if __name__ == '__main__':
    sys.exit(main())
