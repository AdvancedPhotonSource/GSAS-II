"""Embedding model selection for Query-GSAS.

Default: BAAI/bge-base-en-v1.5 via fastembed (ONNX, no PyTorch required).
Override: set GSAS_QUERY_EMBED_MODEL=minilm to use all-MiniLM-L6-v2 instead.
"""

import os
from functools import lru_cache

EMBED_MODEL_ENV = os.environ.get("GSAS_QUERY_EMBED_MODEL", "bge-base")

BGE_BASE_MODEL = "BAAI/bge-base-en-v1.5"


class _FastEmbedWrapper:
    """ChromaDB-compatible wrapper around a fastembed TextEmbedding model."""

    def __init__(self, model_name: str):
        from fastembed import TextEmbedding  # type: ignore
        self._model = TextEmbedding(model_name)
        self.model_name = model_name

    def __call__(self, input: list[str]) -> list[list[float]]:
        return [emb.tolist() for emb in self._model.embed(input)]


@lru_cache(maxsize=1)
def get_embedding_function():
    """Return the configured embedding function (cached after first call)."""
    if EMBED_MODEL_ENV == "minilm":
        from chromadb.utils.embedding_functions import DefaultEmbeddingFunction
        return DefaultEmbeddingFunction()

    try:
        return _FastEmbedWrapper(BGE_BASE_MODEL)
    except ImportError:
        print(
            "Warning: fastembed is not installed. "
            "Install with: pip install fastembed\n"
            "Falling back to all-MiniLM-L6-v2 (lower retrieval quality)."
        )
        from chromadb.utils.embedding_functions import DefaultEmbeddingFunction
        return DefaultEmbeddingFunction()
