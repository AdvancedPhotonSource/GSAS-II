"""
GSAS-II Documentation Assistant — wxPython desktop dialog.

Standalone:
    gsas-query --gui
    python -m gsas_query.gui

GSAS-II Help menu integration:
    def OnDocAssistant(self, event):
        from gsas_query.gui import show_assistant
        show_assistant(self)

Uses MathJaX to display equations using https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js, if available.
"""

import atexit
import html
import os
import re
import subprocess
import sys
import threading
import time
import webbrowser

from ._paths import get_chroma_path


def _get_gsas_font_size(default: int = 10) -> int:
    """Return the font size to use, in points.

    Priority:
      1. specified font size in show_assistant() or GSASQueryDialog() call
      2. GSAS_QUERY_FONT_SIZE env var (explicit override)
      3. *default* (10 pt)
    """
    env = os.environ.get("GSAS_QUERY_FONT_SIZE")
    if env:
        try:
            return int(env)
        except ValueError:
            pass
    return default


def _ollama_url() -> str:
    return os.environ.get("OLLAMA_URL", "http://localhost:11434")


def _format_llm_markdown(text: str) -> str:
    """Convert markdown/math flavored model output into readable plain text.

    The chat transcript is rendered in a wx.TextCtrl, so we normalize common
    markdown and MathJax delimiters rather than showing raw markup.
    """
    if not text:
        return ""

    out = text.replace("\r\n", "\n")

    # Avoid duplicated speaker prefix in the transcript.
    if out.lstrip().lower().startswith("assistant:"):
        out = out.split(":", 1)[1].lstrip()

    # Convert fenced code blocks into plain indented blocks.
    out = re.sub(r"```(?:[\w+-]+)?\n(.*?)```", lambda m: "\n" + "\n".join(
        f"    {line}" if line else "" for line in m.group(1).split("\n")
    ) + "\n", out, flags=re.DOTALL)

    # Convert display math delimiters into visible equation lines.
    out = out.replace("\\[", "\n")
    out = out.replace("\\]", "\n")
    out = out.replace("$$", "\n")

    # Convert inline math delimiters to simple quoted equation text.
    out = out.replace("\\(", "[")
    out = out.replace("\\)", "]")

    # Simplify common markdown markers for plain-text display.
    out = re.sub(r"^\s*[-*]\s+", "• ", out, flags=re.MULTILINE)
    out = re.sub(r"\*\*(.*?)\*\*", r"\1", out)
    out = re.sub(r"__(.*?)__", r"\1", out)
    out = re.sub(r"`([^`]+)`", r"\1", out)

    # Collapse excess blank space introduced by conversions.
    out = re.sub(r"\n{3,}", "\n\n", out)
    return out.strip()


def _strip_speaker_prefix(text: str) -> str:
    """Drop a leading 'Assistant:' emitted by some models."""
    if text.lstrip().lower().startswith("assistant:"):
        return text.split(":", 1)[1].lstrip()
    return text


def _inline_md_to_html(line: str) -> str:
    """Render a small markdown subset into HTML-safe inline content."""
    out = html.escape(line)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", out)
    out = re.sub(r"__(.+?)__", r"<strong>\1</strong>", out)
    return out


def _assistant_text_to_html(text: str) -> str:
    """Convert model output into safe HTML while preserving TeX delimiters."""
    if not text:
        return ""

    src = _strip_speaker_prefix(text).replace("\r\n", "\n")
    parts = re.split(r"```(?:[\w+-]+)?\n(.*?)```", src, flags=re.DOTALL)
    chunks: list[str] = []

    for i, part in enumerate(parts):
        if i % 2 == 1:
            chunks.append(f"<pre><code>{html.escape(part.strip())}</code></pre>")
            continue

        lines = part.split("\n")
        in_list = False
        para_lines: list[str] = []

        def flush_para():
            if para_lines:
                body = "<br>".join(_inline_md_to_html(x) for x in para_lines)
                chunks.append(f"<p>{body}</p>")
                para_lines.clear()

        for raw in lines:
            line = raw.rstrip()
            m = re.match(r"^\s*[-*]\s+(.*)$", line)
            if m:
                flush_para()
                if not in_list:
                    chunks.append("<ul>")
                    in_list = True
                chunks.append(f"<li>{_inline_md_to_html(m.group(1))}</li>")
                continue

            if line.strip() == "":
                flush_para()
                if in_list:
                    chunks.append("</ul>")
                    in_list = False
                continue

            if in_list:
                chunks.append("</ul>")
                in_list = False
            para_lines.append(line)

        flush_para()
        if in_list:
            chunks.append("</ul>")

    return "\n".join(chunks)


def _chat_document(messages: list[dict], font_size: int) -> str:
    """Build the full chat HTML document with MathJax enabled."""
    msg_html: list[str] = []
    for msg in messages:
        role = msg.get("role", "assistant")
        if role == "user":
            body = f"<p>{_inline_md_to_html(msg.get('content', ''))}</p>"
            msg_html.append(
                f'<div class="msg user"><div class="who">You</div><div class="bubble">{body}</div></div>'
            )
        elif role == "system":
            body = f"<p>{_inline_md_to_html(msg.get('content', ''))}</p>"
            msg_html.append(
                f'<div class="msg system"><div class="bubble">{body}</div></div>'
            )
        else:
            body = _assistant_text_to_html(msg.get("content", ""))
            msg_html.append(
                f'<div class="msg assistant"><div class="who">Assistant</div><div class="bubble">{body}</div></div>'
            )

    if not msg_html:
        msg_html.append('<div class="msg system"><div class="bubble"><p>Ask a question about GSAS-II documentation.</p></div></div>')

    return f"""<!doctype html>
<html>
<head>
<meta charset=\"utf-8\">
<style>
  :root {{
    --blue: rgb(26,79,138);
    --bg: rgb(244,246,249);
    --bot: #ffffff;
    --text: #0f172a;
    --muted: #64748b;
    --border: #d1d9e6;
  }}
  html, body {{ margin: 0; padding: 0; background: var(--bg); color: var(--text); }}
    body {{ font: {font_size}pt system-ui, 'Segoe UI', 'Helvetica Neue', Helvetica, Arial, sans-serif; }}
  #chat {{ padding: 10px; }}
  .msg {{ margin: 0 0 10px 0; max-width: 94%; }}
  .msg.user {{ margin-left: auto; }}
    .who {{ color: var(--muted); font-size: .65em; margin: 0 0 3px 2px; }}
    .bubble {{ background: var(--bot); border: 1px solid var(--border); border-radius: 9px; padding: 8px 10px; line-height: 1.3; font-size: .78em; }}
  .msg.user .bubble {{ border-color: rgba(26,79,138,.3); background: #eef5ff; }}
  .msg.system .bubble {{ color: var(--muted); background: #f8fafc; }}
  p {{ margin: 0 0 7px 0; }}
  p:last-child {{ margin-bottom: 0; }}
  ul {{ margin: 6px 0 8px 20px; padding: 0; }}
  li {{ margin: 2px 0; }}
  code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; background: #eef2f7; padding: 1px 4px; border-radius: 4px; }}
  pre {{ margin: 8px 0; background: #1f2937; color: #e5e7eb; padding: 10px; border-radius: 7px; overflow: auto; }}
  pre code {{ background: transparent; padding: 0; color: inherit; }}
</style>
<script>
window.MathJax = {{
  tex: {{
    inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
    displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']]
  }},
  options: {{ skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code'] }}
}};
</script>
<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>
</head>
<body>
<div id=\"chat\">{''.join(msg_html)}</div>
<script>
window.scrollTo(0, document.body.scrollHeight);
if (window.MathJax && window.MathJax.typesetPromise) {{
  window.MathJax.typesetPromise().then(() => window.scrollTo(0, document.body.scrollHeight));
}}
</script>
</body>
</html>
"""


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
    repeat = True
    cmdarg = "serve"
    sec = 6
    while repeat:
        repeat = False
        try:
            proc = subprocess.Popen(
                [_ollama_bin(), cmdarg],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                encoding='UTF-8'
                )
        except FileNotFoundError:
            return None
        for _ in range(sec*2):          # wait up to `sec` s
            time.sleep(0.5)
            if _is_ollama_running():
                # Register a fallback: if GSAS-II exits without closing the dialog
                # (EVT_CLOSE won't fire for child frames), atexit still kills Ollama.
                atexit.register(proc.terminate)
                return proc
            elif cmdarg == "serve":
                # on Mac provided binary runs w/o server. Longer timeout, as 
                # one may need to respond to GUI questions
                for line in proc.stderr:
                    if "serve command not supported" in line:
                        repeat = True
                        sec = 60
                        break
                else:
                    continue
                break
        print(f'Failed to launch Ollama server with "{_ollama_bin()} {cmdarg}"')
        cmdarg = ''  # if repeat, try again without serve        
        proc.terminate()
    return None

try:
    import wx
    import wx.html
    try:
        import wx.html2
        _HAS_WEBVIEW = True
    except Exception:
        _HAS_WEBVIEW = False
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
    def __init__(self, parent, source: dict, number: int, font_size: int = 10):
        super().__init__(parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(parent.GetBackgroundColour())

        url = source.get("url", "")
        title = source.get("title", "Unknown")
        section = source.get("section", "")
        rel = int(source.get("relevance", 0) * 100)

        small = max(font_size - 1, 8)

        # Number label — bold, matches inline [N] in answer text
        num_lbl = wx.StaticText(self, label=f"[{number}]")
        num_lbl.SetForegroundColour(_ACCENT)
        num_lbl.SetFont(wx.Font(wx.FontInfo(small).Bold()))

        # Clickable title + section
        text = title
        if section and section != title:
            text += f"  ›  {section}"
        text += f"  [{rel}%]"

        lnk = wx.StaticText(self, label=text)
        lnk.SetForegroundColour(_BLUE)
        lnk.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        lnk.SetFont(wx.Font(wx.FontInfo(small)))
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

    def __init__(self, parent, fontsize=None):
        super().__init__(
            parent,
            title="GSAS-II Documentation Assistant",
            size=(700, 620),
            style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT,
        )
        self._history: list[dict] = []
        self._messages: list[dict] = []
        self._pending_question = ""
        self._busy = False
        self._ollama_proc: "subprocess.Popen | None" = None
        if fontsize is None:
            self._font_size = _get_gsas_font_size(default=10)
        else:
            self._font_size = fontsize
        self._build_ui()
        if self._chat_web is None:
            self._append_system(
                "Math rendering unavailable: WebView backend is not available in this wxPython environment."
            )
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

        # Chat transcript (WebView + MathJax when available)
        self._chat_web = None
        self._chat = None
        if _HAS_WEBVIEW:
            try:
                self._chat_web = wx.html2.WebView.New(self)
                outer.Add(self._chat_web, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)
                self._render_chat()
            except Exception:
                # Some wx builds expose wx.html2 but lack a working runtime backend.
                self._chat_web = None

        if self._chat_web is None:
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

        # Progress indicator — shown while the LLM is generating
        self._gauge = wx.Gauge(self, range=50, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
        self._gauge.Hide()
        outer.Add(self._gauge, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)
        self._gauge_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, lambda _: self._gauge.Pulse(), self._gauge_timer)

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
            from .rag import _effective_backend
            self._index_lbl.SetLabel(f"{count:,} chunks · {_effective_backend()}")

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
        self._pending_question = question
        self._busy = True
        self._send_btn.Disable()
        self._status.SetStatusText("Thinking…")
        self._gauge.Show()
        self._gauge_timer.Start(80)   # pulse every 80 ms
        wx.BeginBusyCursor()
        self.Layout()
        self._append_user(question)
        self._clear_sources()
        _QueryThread(self, question, self._history).start()

    def _on_clear(self, event):
        self._history.clear()
        self._messages.clear()
        if self._chat_web is not None:
            self._render_chat()
        else:
            self._chat.SetValue("")
        self._clear_sources()
        self._status.SetStatusText("History cleared.")

    def _on_query_done(self, result: dict):
        self._busy = False
        self._send_btn.Enable()
        self._gauge_timer.Stop()
        self._gauge.Hide()
        wx.EndBusyCursor()
        self.Layout()
        self._status.SetStatusText("")
        answer = result.get("answer", "")
        self._append_assistant(answer)
        self._history.append({"role": "user", "content": self._pending_question})
        self._history.append({"role": "assistant", "content": answer})
        self._show_sources(result.get("citations", {}))

    # ── Chat helpers ───────────────────────────────────────────────────────────

    def _render_chat(self):
        if self._chat_web is not None:
            self._chat_web.SetPage(_chat_document(self._messages, self._font_size), "")

    def _append_user(self, text: str):
        self._messages.append({"role": "user", "content": text})
        if self._chat_web is not None:
            self._render_chat()
            return
        self._chat.SetDefaultStyle(
            wx.TextAttr(_BLUE, font=wx.Font(wx.FontInfo(self._font_size).Bold()))
        )
        self._chat.AppendText(f"You: {text}\n\n")
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK))

    def _append_assistant(self, text: str):
        self._messages.append({"role": "assistant", "content": text})
        if self._chat_web is not None:
            self._render_chat()
            return
        formatted = _format_llm_markdown(text)
        self._chat.SetDefaultStyle(
            wx.TextAttr(wx.BLACK, font=wx.Font(wx.FontInfo(self._font_size)))
        )
        self._chat.AppendText(f"Assistant: {formatted}\n\n")
        self._chat.SetDefaultStyle(wx.TextAttr(wx.BLACK))

    def _append_system(self, text: str):
        self._messages.append({"role": "system", "content": text})
        if self._chat_web is not None:
            self._render_chat()
            return
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
            item = _SourcePanel(self._src_panel, citations[key],
                                number=int(key), font_size=self._font_size)
            self._src_inner.Add(item, 0, wx.EXPAND | wx.BOTTOM, 4)
        self._src_panel.FitInside()
        self._src_panel.Layout()
        self.Layout()


# ── Public API ─────────────────────────────────────────────────────────────────

_instance: GSASQueryDialog | None = None


def show_assistant(parent=None,fontsize=None) -> GSASQueryDialog:
    """
    Show the GSAS-II Documentation Assistant dialog.

    Call from GSAS-II's Help menu:
        from gsas_query.gui import show_assistant
        show_assistant(self)

    If already open, brings the window to front instead of opening a duplicate.
    """
    global _instance
    if _instance is None or not _instance.IsShown():
        _instance = GSASQueryDialog(parent,fontsize=fontsize)
        _instance.Show()
    else:
        _instance.Raise()
    return _instance


# ── Standalone entry point ─────────────────────────────────────────────────────

if __name__ == "__main__":
    app = wx.App(False)
    show_assistant()
    app.MainLoop()
