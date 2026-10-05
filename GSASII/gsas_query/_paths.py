"""Central path management for user data (chroma index) and package assets."""

import os
from pathlib import Path

# Suppress HuggingFace tokenizer fork warning and OpenMP threading —
# both must be set before sentence-transformers / loky are imported.
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
os.environ.setdefault("OMP_NUM_THREADS", "1")


def get_data_dir() -> Path:
    """Return (and create) the user data directory for gsas_query."""
    d = Path(os.environ.get("GSAS_QUERY_DATA_DIR",
                            Path.home() / ".GSASII" / "query_gsas2"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def get_chroma_path() -> str:
    return str(get_data_dir() / "chroma_db")


def get_static_dir() -> Path:
    return Path(__file__).parent / "static"
