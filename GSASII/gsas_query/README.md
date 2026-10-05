# gsas2-query — GSAS-II Documentation Assistant

Semantic search + AI answers over the full GSAS-II documentation set:
**129 HTML pages** (home, help, and all 62 tutorials) plus the
**Programmer's Guide** and **Powder Crystallography** book PDFs.

Answers include **inline citations** — every `[N]` in the response is a
clickable link to the exact documentation section that supported that sentence.

All embedding and retrieval runs on your machine — no data is sent externally
unless you explicitly choose the Anthropic API backend.

---

## Installation

### pip

```bash
pip install gsas2-query
```

Or install directly from source:

```bash
pip install git+https://github.com/AdvancedPhotonSource/Query-GSAS-II.git
```

### conda

```bash
conda install -c conda-forge gsas2-query
```

### Into an existing GSAS-II environment

```bash
# Activate the GSAS-II conda environment first, then:
pip install gsas2-query
```

### Development / editable install

```bash
git clone https://github.com/AdvancedPhotonSource/Query-GSAS-II.git
cd Query-GSAS-II
pip install -e ".[dev]"
```

---

## Ollama — free local LLM (recommended when llama.cpp is not installed)

```bash
brew install ollama          # macOS; see https://ollama.com for other platforms
ollama serve &               # start the local server
ollama pull llama3           # ~5 GB one-time download
```

Other model options: `llama3:70b` (better quality, ~40 GB), `mistral` (faster, ~4 GB).

---

## llama-cpp-python — in-process local LLM (auto-selected when installed)

`llama-cpp-python` runs a GGUF model directly inside the Python process — no
separate server required. It is cross-platform and available from conda-forge.

When `llama-cpp-python` is importable, gsas2-query selects it **automatically**
without any `LLM_BACKEND` setting.

### Install

```bash
conda install -c conda-forge llama-cpp-python   # recommended
# or:
pip install "gsas2-query[llama_cpp]"
```

### Download a GGUF model

Download any GGUF-format model from Hugging Face, for example:

```bash
# Llama 3.1 8B (Q4_K_M quantisation, ~5 GB):
wget https://huggingface.co/bartowski/Meta-Llama-3.1-8B-Instruct-GGUF/resolve/main/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf
```

### Configure

Set the model path in your `.env` or environment:

```env
LLAMA_CPP_MODEL=/path/to/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf
# Optional tuning:
# LLAMA_CPP_N_CTX=4096       # context window (default: 4096)
# LLAMA_CPP_MAX_TOKENS=1500  # max tokens to generate (default: 1500)
```

---

## First-time setup — index the documentation

Run once, or again when the docs are updated. Fetches ~130 web pages and 2 PDFs,
embeds everything locally. Takes ~10–20 minutes.

By default the index is stored in `~/.GSASII/gsas_query/chroma_db`
(see [User data directory](#userdatadirectory) for details).

```bash
gsas2-query --setup              # all sources (HTML + PDFs)
gsas2-query --setup --html-only  # skip PDFs, faster (~5 min)
gsas2-query --setup --reset      # drop index and rebuild from scratch
```

HTML-only index is 21 Mb. With GSAS-II Programmer's Guide, 41 Mb; with
textbook as well, 53 Mb.

---

## Usage

### Command-line — single question

```bash
gsas2-query "How do I set up a sequential refinement?"
gsas2-query "What parameters control the background in Rietveld?"
gsas2-query "How do I export a CIF for publication?"
```

### Command-line — interactive REPL

Multi-turn conversation that remembers previous questions in the session.

```bash
gsas2-query
```

```
GSAS-II Documentation Assistant
════════════════════════════════════════════════════════════════════════════════

Knowledge base: 3,847 indexed chunks.
LLM backend: ollama
Type your question and press Enter. 'clear' resets history, 'quit' exits.

────────────────────────────────────────────────────────────────────────────────
You: How do I constrain lattice parameters?
Thinking…
Assistant: To constrain lattice parameters in GSAS-II, open the Constraints
tab in the Phase panel [1]…

Sources:
  [94%] Help: Phase General  ›  Constraints
         https://advancedphotonsource.github.io/GSAS-II-tutorials/help/phasegeneral.html
```

Commands inside the REPL: `clear` (reset history), `quit` / `exit` (exit).

### Desktop GUI

Opens a standalone floating dialog — stays open while you work in GSAS-II.

```bash
gsas2-query --gui
```

### Embed in GSAS-II Help menu

Add one call to the GSAS-II menu handler (e.g. in `GSASIIctrl.py`):

```python
def OnDocAssistant(self, event):
    try:
        from gsas_query.gui import show_assistant
        show_assistant(self)          # self = GSAS-II main frame
    except ImportError:
        wx.MessageBox(
            "GSAS-II Assistant not installed.\n"
            "Run: pip install gsas2-query",
            "Not available"
        )
```

`show_assistant()` is idempotent — calling it a second time raises the existing
window rather than opening a duplicate.

### Web UI

```bash
gsas2-query-web                          # serves on 0.0.0.0:8000
HOST=127.0.0.1 PORT=8765 gsas2-query-web
```

Then open `http://localhost:8000` in a browser. The web UI shows inline citations
as clickable superscript links and source chips below each answer.

---

## CLI flags reference

| Flag | Description |
|---|---|
| `--setup` | Index all documentation sources |
| `--setup --reset` | Drop the existing index and rebuild |
| `--setup --html-only` | Index HTML only, skip PDFs |
| `--gui` | Open the wxPython desktop assistant |
| `--backend ollama\|anthropic\|retrieval\|llama_cpp` | Override `LLM_BACKEND` env var |
| `--model <name>` | Override OLLAMA_MODEL, ANTHROPIC_MODEL, or LLAMA_CPP_MODEL |
| `--stats` | Show chunk count, backend, and DB path |

---

## LLM backend configuration

Backend selection precedence:

1. `LLM_BACKEND` env var (or `--backend` CLI flag) **always takes effect when set**.
2. If `LLM_BACKEND` is **not set** and `llama-cpp-python` is importable, `llama_cpp` is selected automatically.
3. Otherwise the default is `ollama`.

| Backend | Config | Notes |
|---|---|---|
| `llama_cpp` (auto) | `LLAMA_CPP_MODEL=/path/to/model.gguf` | In-process, no daemon needed; auto-selected when llama-cpp-python is installed |
| `ollama` (default) | `OLLAMA_MODEL=llama3` | Free, fully local — no data leaves the network |
| `anthropic` | `ANTHROPIC_API_KEY=sk-ant-…` | Better answers; queries sent to Anthropic |
| `retrieval` | — | No LLM — returns raw matched chunks; useful offline or for testing |

```bash
gsas2-query --backend retrieval "What is Le Bail extraction?"
gsas2-query --backend ollama --model mistral "How do I index peaks?"
gsas2-query --backend llama_cpp --model /path/to/model.gguf "How do I refine a structure?"
```

---

## Inline citations

When using Ollama, llama_cpp, or Anthropic backends, answers contain `[N]` markers inline.
In the **web UI** these render as clickable superscript links opening the exact
source section. In the **CLI**, source URLs are listed below the answer with
relevance scores.

---

## Doc sources

| Category | Count |
|---|---|
| Home / installation pages | 22 |
| Help pages (all sections) | 42 |
| Tutorials | 62 |
| Programmer's Guide (PDF, readthedocs) | 1 |
| Powder Crystallography book (PDF, auto-fetches latest release) | 1 |
| **Total** | **128 sources** |

All HTML sources are fetched from
`https://advancedphotonsource.github.io/GSAS-II-tutorials/`.
The book PDF is fetched from the latest GitHub release of
[briantoby/PowderCrystallography](https://github.com/briantoby/PowderCrystallography/releases).

---

## Re-indexing when docs update

```bash
gsas2-query --setup --reset
```

Or trigger via the web API (requires `ADMIN_KEY` set in `.env`):

```bash
curl -X POST http://localhost:8000/ingest -H "X-Admin-Key: your-key"
```

---

## Security and deployment notes

- All embeddings and vector search run locally (sentence-transformers, ChromaDB).
- Ollama runs entirely on-premises — no queries leave the network.
- The `anthropic` backend sends question text and retrieved doc chunks to Anthropic.
  Do not use it in air-gapped or data-sensitive environments.
- The web server applies per-IP rate limiting (default 30 req/min, configurable
  via `RATE_LIMIT_RPM` in `.env`).
- The `/ingest` endpoint is protected by `X-Admin-Key`; leave `ADMIN_KEY` blank
  to disable remote re-indexing.

---

## User data directory

The ChromaDB index is stored at `~/.GSASII/gsas_query/chroma_db` by default.
Override with the `GSAS_QUERY_DATA_DIR` environment variable:

```bash
GSAS_QUERY_DATA_DIR=/data/gsas_query gsas2-query --setup
```
