"""
GSAS-II Documentation Assistant — wxPython desktop dialog.

Standalone:
    gsas-query --gui
    python -m gsas_query.gui

GSAS-II Help menu integration:
    def OnDocAssistant(self, event):
        from gsas_query.gui import show_assistant
        show_assistant(self)
"""

import atexit
import os
import subprocess
import sys
import threading
import time
import webbrowser

from ._paths import get_chroma_path


def _ollama_url() -> str:
    return os.environ.get("OLLAMA_URL", "http://localhost:11434")


def _ollama_bin() -> str:
    """Return ollama executable path (optionally overridden by env)."""
    return os.environ.get("OLLAMA_BIN", "ollama")


def _is_ollama_running() -> bool:
    try:
        import httpx
        return httpx.get(f"{_ollama_url()}/api/tags", timeout=2).status_code == 200
    except Exception:
        return False


def _launch_ollama() -> "subprocess.Popen | None":
    """Start `ollama serve` if not already running. Returns the process we started, or None."""
    if _is_ollama_running():
        return None
    try:
        proc = subprocess.Popen(
            [_ollama_bin(), "serve"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        return None
    for _ in range(12):          # wait up to 6 s
        time.sleep(0.5)
        if _is_ollama_running():
            # Register a fallback: if GSAS-II exits without closing the dialog
            # (EVT_CLOSE won't fire for child frames), atexit still kills Ollama.
            atexit.register(proc.terminate)
            return proc
    proc.terminate()
    return None

try:
    import wx
    import wx.html
except ImportError:
    print("wxPython is required for the GUI. It is included with GSAS-II.")
    sys.exit(1)


# ── Colour constants ───────────────────────────────────────────────────────────

_BLUE   = wx.Colour(26, 79, 138)     # Argonne blue
_ACCENT = wx.Colour(232, 119, 34)    # APS orange
_BG     = wx.Colour(244, 246, 249)
_BOT_BG = wx.Colour(255, 255, 255)
_BORDER = wx.Colour(209, 217, 230)
_MUTED  = wx.Colour(107, 114, 128)


# ── Background query thread ────────────────────────────────────────────────────

class _QueryThread(threading.Thread):
    def __init__(self, parent, question: str, history: list):
        super().__init__(daemon=True)
        self.parent = parent
        self.question = question
        self.history = list(history)

    def run(self):
        try:
            from .rag import answer_question
            result = answer_question(self.question, self.history)
        except Exception as e:
            result = {"answer": f"Error: {e}", "sources": []}
        wx.CallAfter(self.parent._on_query_done, result)


# ── Source link panel ──────────────────────────────────────────────────────────

class _SourcePanel(wx.Panel):
    def __init__(self, parent, source: dict, number: int):
        super().__init__(parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(parent.GetBackgroundColour())

        url = source.get("url", "")
        title = source.get("title", "Unknown")
        section = source.get("section", "")
        rel = int(source.get("relevance", 0) * 100)

        # Number label — bold, matches inline [N] in answer text
        num_lbl = wx.StaticText(self, label=f"[{number}]")
        num_lbl.SetForegroundColour(_ACCENT)
        nf = num_lbl.GetFont()
        nf.SetWeight(wx.FONTWEIGHT_BOLD)
        nf.SetPointSize(nf.GetPointSize() - 1)
        num_lbl.SetFont(nf)

        # Clickable title + section
        text = title
        if section and section != title:
            text += f"  ›  {section}"
        text += f"  [{rel}%]"

        lnk = wx.StaticText(self, label=text)
        lnk.SetForegroundColour(_BLUE)
        lnk.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        f = lnk.GetFont()
        f.SetPointSize(f.GetPointSize() - 1)
        lnk.SetFont(f)
        lnk.Bind(wx.EVT_LEFT_UP, lambda e: webbrowser.open(url))

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(num_lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        sizer.Add(lnk, 1, wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(sizer)


# ── Main dialog ────────────────────────────────────────────────────────────────

class GSASQueryDialog(wx.Frame):
    """
    Modeless frame — stays open while the user works in GSAS-II.
    Call show_assistant() rather than instantiating directly.
    """

    def __init__(self, parent):
        super().__init__(
            parent,
            title="GSAS-II Documentation Assistant",
            size=(700, 620),
            style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT,
        )
        self._history: list[dict] = []
        self._busy = False
        self._ollama_proc: "subprocess.Popen | None" = None
        self._build_ui()
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)
        self._check_index()
        self._ensure_ollama()

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        self.SetBackgroundColour(_BG)
        outer = wx.BoxSizer(wx.VERTICAL)

        # Header
        header = wx.Panel(self)
        header.SetBackgroundColour(_BLUE)
        h_sizer = wx.BoxSizer(wx.HORIZONTAL)
        title_lbl = wx.StaticText(header, label="GSAS-II Documentation Assistant")
        title_lbl.SetForegroundColour(wx.WHITE)
        f = title_lbl.GetFont()
        f.SetWeight(wx.FONTWEIGHT_BOLD)
        f.SetPointSize(f.GetPointSize() + 2)
        title_lbl.SetFont(f)
        self._index_lbl = wx.StaticText(header, label="")
        self._index_lbl.SetForegroundColour(wx.Colour(200, 220, 255))
        small = self._index_lbl.GetFont()
        small.SetPointSize(small.GetPointSize() - 1)
        self._index_lbl.SetFont(small)
        h_sizer.Add(title_lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 14)
        h_sizer.AddStretchSpacer()
        h_sizer.Add(self._index_lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 14)
        header.SetSizer(h_sizer)
        header.SetMinSize((-1, 46))
        outer.Add(header, 0, wx.EXPAND)

        # Chat transcript
        self._chat = wx.TextCtrl(
            self, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH2 | wx.BORDER_NONE
        )
        self._chat.SetBackgroundColour(_BOT_BG)
        self._chat.SetMinSize((-1, 300))
        outer.Add(self._chat, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)

        # Sources panel
        src_box = wx.StaticBox(self, label="Sources")
        src_box.SetForegroundColour(_MUTED)
        self._src_sizer = wx.StaticBoxSizer(src_box, wx.VERTICAL)
        self._src_panel = wx.ScrolledWindow(self, style=wx.BORDER_NONE)
        self._src_panel.SetScrollRate(0, 12)
        self._src_panel.SetBackgroundColour(_BG)
        self._src_inner = wx.BoxSizer(wx.VERTICAL)
        self._src_panel.SetSizer(self._src_inner)
        self._src_sizer.Add(self._src_panel, 1, wx.EXPAND | wx.ALL, 4)
        self._src_panel.SetMinSize((-1, 90))
        outer.Add(self._src_sizer, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)

        # Input row
        in_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self._input = wx.TextCtrl(
            self, style=wx.TE_PROCESS_ENTER | wx.TE_MULTILINE, size=(-1, 54),
        )
        self._input.SetHint("Ask a question about GSAS-II…")
        self._input.Bind(wx.EVT_TEXT_ENTER, self._on_send)
        self._input.Bind(wx.EVT_KEY_DOWN, self._on_key)

        self._send_btn = wx.Button(self, label="Ask", size=(70, 54))
        self._send_btn.SetBackgroundColour(_BLUE)
        self._send_btn.SetForegroundColour(wx.WHITE)
        self._send_btn.Bind(wx.EVT_BUTTON, self._on_send)

        clear_btn = wx.Button(self, label="Clear", size=(60, 54))
        clear_btn.Bind(wx.EVT_BUTTON, self._on_clear)

        in_sizer.Add(self._input, 1, wx.EXPAND | wx.RIGHT, 6)
        in_sizer.Add(self._send_btn, 0)
        in_sizer.Add(clear_btn, 0, wx.LEFT, 4)
        outer.Add(in_sizer, 0, wx.EXPAND | wx.ALL, 10)

        self._status = self.CreateStatusBar()
        self.SetSizer(outer)
        self.Layout()

    # ── Index check ────────────────────────────────────────────────────────────

    def _check_index(self):
        def _check():
            try:
                import chromadb
                client = chromadb.PersistentClient(path=get_chroma_path())
                count = client.get_or_create_collection("gsasii_docs").count()
            except Exception:
                count = 0
            wx.CallAfter(self._set_index_status, count)
        threading.Thread(target=_check, daemon=True).start()

    def _ensure_ollama(self):
        from .rag import _effective_backend
        if _effective_backend() != "ollama":
            return

        # If an Ollama server is already reachable (possibly started outside this env),
        # do not launch another process from PATH.
        if _is_ollama_running():
            self._status.SetStatusText("Connected to Ollama")
            return

        def _start():
            wx.CallAfter(self._status.SetStatusText, "Checking Ollama…")
            proc = _launch_ollama()
            if proc is not None:
                self._ollama_proc = proc
                wx.CallAfter(self._status.SetStatusText, "Ollama started")
            else:
                wx.CallAfter(
                    self._status.SetStatusText,
                    "Ollama not found. Start it manually: ollama serve",
                )


        threading.Thread(target=_start, daemon=True).start()

    def _on_close(self, event):
        global _instance
        _instance = None
        if self._ollama_proc is not None:
            atexit.unregister(self._ollama_proc.terminate)  # cancel the fallback
            self._ollama_proc.terminate()
            self._ollama_proc = None
        self.Destroy()

    def _set_index_status(self, count: int):
        if count == 0:
            self._index_lbl.SetLabel("Not indexed — run: gsas-query --setup")
            self._append_system(
                "The knowledge base is empty.\n"
                "Run 'gsas-query --setup' to index all GSAS-II documentation (~10 min)."
            )
        else:
            backend = os.environ.get("LLM_BACKEND", "ollama")
            self._index_lbl.SetLabel(f"{count:,} chunks · {backend}")

    # ── Event handlers ─────────────────────────────────────────────────────────

    def _on_key(self, event):
        if event.GetKeyCode() == wx.WXK_RETURN and not event.ShiftDown():
            self._on_send(event)
        else:
            event.Skip()

    def _on_send(self, event):
        question = self._input.GetValue().strip()
        if not question or self._busy:
            return
        self._input.SetValue("")
        self._busy = True
        self._send_btn.Disable()
        self._status.SetStatusText("Thinking…")
        self._append_user(question)
        self._clear_sources()
        _QueryThread(self, question, self._history).start()

    def _on_clear(self, event):
        self._history.clear()
        self._chat.SetValue("")
        self._clear_sources()
        self._status.SetStatusText("History cleared.")

    def _on_query_done(self, result: dict):
        self._busy = False
        self._send_btn.Enable()
        self._status.SetStatusText("")
        answer = result.get("answer", "")
        self._append_assistant(answer)
        self._history.append({"role": "user",
                               "content": self._chat.GetValue().split("\n")[-2]})
        self._history.append({"role": "assistant", "content": answer})
        self._show_sources(result.get("citations", {}))

    # ── Chat helpers ───────────────────────────────────────────────────────────

    def _append_user(self, text: str):
        self._chat.SetDefaultStyle(wx.TextAttr(_BLUE, font=wx.Font(wx.FontInfo(10).Bold())))
        self._chat.AppendText(f"You: {text}\n\n")
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK))

    def _append_assistant(self, text: str):
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK, font=wx.Font(wx.FontInfo(10))))
        self._chat.AppendText(f"Assistant: {text}\n\n")
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK))

    def _append_system(self, text: str):
        self._chat.SetDefaultStyle(wx.TextAttr(_MUTED))
        self._chat.AppendText(f"{text}\n\n")
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK))

    # ── Source helpers ─────────────────────────────────────────────────────────

    def _clear_sources(self):
        self._src_inner.Clear(delete_windows=True)
        self._src_panel.Layout()

    def _show_sources(self, citations: dict):
        self._clear_sources()
        for key in sorted(citations, key=lambda k: int(k)):
            item = _SourcePanel(self._src_panel, citations[key], number=int(key))
            self._src_inner.Add(item, 0, wx.EXPAND | wx.BOTTOM, 4)
        self._src_panel.FitInside()
        self._src_panel.Layout()
        self.Layout()


# ── Public API ─────────────────────────────────────────────────────────────────

_instance: GSASQueryDialog | None = None


def show_assistant(parent=None) -> GSASQueryDialog:
    """
    Show the GSAS-II Documentation Assistant dialog.

    Call from GSAS-II's Help menu:
        from gsas_query.gui import show_assistant
        show_assistant(self)

    If already open, brings the window to front instead of opening a duplicate.
    """
    global _instance
    if _instance is None or not _instance.IsShown():
        _instance = GSASQueryDialog(parent)
        _instance.Show()
    else:
        _instance.Raise()
    return _instance


# ── Standalone entry point ─────────────────────────────────────────────────────

if __name__ == "__main__":
    app = wx.App(False)
    show_assistant()
    app.MainLoop()
