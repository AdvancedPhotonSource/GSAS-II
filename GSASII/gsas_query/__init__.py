"""
GSAS-II Documentation Assistant

Semantic search + AI answers over GSAS-II tutorials, help pages, and PDFs.
All embedding and retrieval runs locally; the LLM backend is configurable.

Quick start:
    from gsas_query.rag import answer_question
    result = answer_question("How do I set up a sequential refinement?", [])

wxPython GUI (for GSAS-II Help menu integration):
    from gsas_query.gui import show_assistant
    show_assistant(parent_wx_window)
"""

from gsas_query._paths import get_chroma_path, get_data_dir, get_static_dir
from gsas_query.rag import answer_question

__all__ = ["answer_question", "get_chroma_path", "get_data_dir", "get_static_dir"]
__version__ = "0.2.0"
