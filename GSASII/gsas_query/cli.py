"""
GSAS-II Documentation Assistant — command-line interface.

First-time setup:
    gsas-query --setup
    gsas-query --setup --reset
    gsas-query --setup --html-only

Ask a single question:
    gsas-query "How do I set up a sequential refinement?"

Interactive REPL:
    gsas-query
"""

import os
import sys
import textwrap

from ._paths import get_chroma_path

WIDTH = 80


def _hr(char="─", width=WIDTH):
    print(char * width)


def _wrap(text: str, indent: int = 0) -> str:
    prefix = " " * indent
    return textwrap.fill(text, width=WIDTH, initial_indent=prefix, subsequent_indent=prefix)


def _print_answer(result: dict):
    answer = result.get("answer", "")
    sources = result.get("sources", [])

    print()
    for para in answer.split("\n"):
        if para.strip():
            print(_wrap(para))
        else:
            print()

    if sources:
        print()
        _hr("·")
        print("Sources:")
        seen = set()
        for s in sources:
            key = s["url"]
            if key in seen:
                continue
            seen.add(key)
            rel = int(s.get("relevance", 0) * 100)
            label = f"  [{rel}%] {s['title']}"
            if s.get("section") and s["section"] != s["title"]:
                label += f"  ›  {s['section']}"
            print(label)
            print(f"         {s['url']}")
    print()


def _collection_count() -> int:
    try:
        import chromadb
        client = chromadb.PersistentClient(path=get_chroma_path())
        return client.get_or_create_collection("gsasii_docs").count()
    except Exception:
        return 0


def run_setup(reset: bool = False, html_only: bool = False):
    from . import ingest
    import argparse

    print("Starting ingestion. This may take 10–20 minutes.")
    print(f"Index will be stored at: {get_chroma_path()}\n")

    argv_backup = sys.argv
    sys.argv = ["gsas-query"]
    if reset:
        sys.argv.append("--reset")
    if html_only:
        sys.argv.append("--html-only")
    try:
        ingest.main()
    finally:
        sys.argv = argv_backup


def ask(question: str, history: list | None = None) -> dict:
    from .rag import answer_question
    return answer_question(question, history or [])


def interactive():
    from .rag import _effective_backend
    count = _collection_count()
    if count == 0:
        print("\nWarning: the knowledge base is empty.")
        print("Run:  gsas-query --setup\n")
    else:
        print(f"\nKnowledge base: {count:,} indexed chunks.")

    backend = _effective_backend()
    print(f"LLM backend: {backend}")
    print("Type your question and press Enter. 'clear' resets history, 'quit' exits.\n")
    _hr()

    history: list[dict] = []

    while True:
        try:
            question = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye.")
            break

        if not question:
            continue
        if question.lower() in {"quit", "exit", "q"}:
            print("Goodbye.")
            break
        if question.lower() in {"clear", "reset"}:
            history.clear()
            print("(History cleared)\n")
            continue

        print("Thinking…", end="", flush=True)
        try:
            result = ask(question, history)
        except Exception as e:
            print(f"\rError: {e}\n")
            continue

        print("\r" + " " * 12 + "\r", end="")
        print("Assistant:", end="")
        _print_answer(result)

        history.append({"role": "user", "content": question})
        history.append({"role": "assistant", "content": result.get("answer", "")})


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="GSAS-II Documentation Assistant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("question", nargs="?", help="Question to ask (omit for interactive mode)")
    parser.add_argument("--setup", action="store_true", help="Index documentation (first-time setup)")
    parser.add_argument("--reset", action="store_true", help="Drop and rebuild the index")
    parser.add_argument("--html-only", action="store_true", help="Skip PDFs during setup")
    parser.add_argument("--backend", choices=["ollama", "anthropic", "retrieval", "llama_cpp"],
                        help="Override LLM_BACKEND env var")
    parser.add_argument("--model", help="Override OLLAMA_MODEL or ANTHROPIC_MODEL")
    parser.add_argument("--stats", action="store_true", help="Show index statistics and exit")
    parser.add_argument("--gui", action="store_true", help="Open the wxPython desktop assistant")
    args = parser.parse_args()

    if args.backend:
        os.environ["LLM_BACKEND"] = args.backend
    if args.model:
        effective = os.environ.get("LLM_BACKEND", "")
        if effective == "anthropic":
            os.environ["ANTHROPIC_MODEL"] = args.model
        elif effective == "llama_cpp":
            os.environ["LLAMA_CPP_MODEL"] = args.model
        else:
            os.environ["OLLAMA_MODEL"] = args.model

    print("GSAS-II Documentation Assistant")
    _hr("═")

    if args.stats:
        from .rag import _effective_backend
        count = _collection_count()
        print(f"Indexed chunks : {count:,}")
        print(f"LLM backend    : {_effective_backend()}")
        print(f"Chroma DB      : {get_chroma_path()}")
        return

    if args.setup:
        run_setup(reset=args.reset, html_only=args.html_only)
        return

    if args.gui:
        from .gui import show_assistant
        try:
            import wx
        except ImportError:
            print("wxPython is required for the GUI. Install it or use GSAS-II's Python.")
            sys.exit(1)
        app = wx.App(False)
        show_assistant()
        app.MainLoop()
        return

    if args.question:
        count = _collection_count()
        if count == 0:
            print("Knowledge base is empty. Run:  gsas-query --setup")
            sys.exit(1)
        print(f"({count:,} chunks indexed)\n")
        result = ask(args.question)
        _print_answer(result)
    else:
        interactive()


if __name__ == "__main__":
    main()
