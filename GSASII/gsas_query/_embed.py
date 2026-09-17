"""Embedding model for Query-GSAS-II.

Default: BAAI/bge-base-en-v1.5 via onnxruntime + tokenizers (no fastembed,
no PyTorch). The quantized ONNX model (~110 MB) is downloaded once to the
GSAS-II data directory and reused on every subsequent launch.

Override: set GSAS_QUERY_EMBED_MODEL=minilm to use chromadb's built-in
all-MiniLM-L6-v2 (no download required, lower retrieval quality).
"""

import os
import urllib.request
from functools import lru_cache
from pathlib import Path

from ._paths import get_data_dir

_BGE_MODEL = "BAAI/bge-base-en-v1.5"

# Files fetched on first use; keyed by local filename → remote URL.
# Xenova repo hosts the quantized ONNX (~110 MB); tokenizer from BAAI origin.
_MODEL_FILES = {
    "model.onnx": (
        "https://huggingface.co/Xenova/bge-base-en-v1.5"
        "/resolve/main/onnx/model_quantized.onnx"
    ),
    "tokenizer.json": (
        "https://huggingface.co/BAAI/bge-base-en-v1.5"
        "/resolve/main/tokenizer.json"
    ),
}


def get_model_dir() -> Path:
    d = get_data_dir() / "models" / "bge-base-en-v1.5"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_model(model_dir: Path) -> None:
    print(f"Downloading {_BGE_MODEL} ONNX model (~110 MB, one-time download)...")
    for filename, url in _MODEL_FILES.items():
        dest = model_dir / filename
        if dest.exists():
            continue
        print(f"  {filename} ...", end="", flush=True)
        tmp = dest.with_suffix(".tmp")
        try:
            urllib.request.urlretrieve(url, tmp)
            tmp.rename(dest)
            print(" done")
        except Exception as exc:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"Could not download {filename} from {url}: {exc}\n"
                "Check your internet connection or pre-populate "
                f"{model_dir} for offline deployments."
            ) from exc


class _OnnxBgeWrapper:
    """onnxruntime + tokenizers embedding wrapper, ChromaDB-compatible."""

    def __init__(self, model_dir: Path):
        import numpy as np
        import onnxruntime as rt
        from tokenizers import Tokenizer

        self._np = np

        tok = Tokenizer.from_file(str(model_dir / "tokenizer.json"))
        tok.enable_padding()          # uses [PAD] token defined in tokenizer.json
        tok.enable_truncation(max_length=512)
        self._tok = tok

        opts = rt.SessionOptions()
        opts.inter_op_num_threads = 1
        opts.intra_op_num_threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
        self._sess = rt.InferenceSession(
            str(model_dir / "model.onnx"),
            sess_options=opts,
            providers=["CPUExecutionProvider"],
        )
        self._has_token_type = any(
            inp.name == "token_type_ids" for inp in self._sess.get_inputs()
        )

    def __call__(self, input: list[str]) -> list[list[float]]:
        np = self._np
        encodings = self._tok.encode_batch(input)
        ids  = np.array([e.ids            for e in encodings], dtype=np.int64)
        mask = np.array([e.attention_mask for e in encodings], dtype=np.int64)

        feeds: dict = {"input_ids": ids, "attention_mask": mask}
        if self._has_token_type:
            feeds["token_type_ids"] = np.zeros_like(ids)

        last_hidden = self._sess.run(None, feeds)[0]   # (B, L, 768)

        # Mean pool over non-padding tokens, then L2-normalise.
        fmask = mask[:, :, np.newaxis].astype(np.float32)
        emb   = (last_hidden * fmask).sum(axis=1) / fmask.sum(axis=1).clip(1e-12)
        emb  /= np.linalg.norm(emb, axis=1, keepdims=True).clip(1e-12)

        return emb.tolist()


@lru_cache(maxsize=1)
def get_embedding_function():
    """Return the configured embedding function (cached after first call)."""
    if os.environ.get("GSAS_QUERY_EMBED_MODEL") == "minilm":
        from chromadb.utils.embedding_functions import DefaultEmbeddingFunction
        return DefaultEmbeddingFunction()

    model_dir = get_model_dir()
    missing = [f for f in _MODEL_FILES if not (model_dir / f).exists()]
    if missing:
        _download_model(model_dir)
    return _OnnxBgeWrapper(model_dir)
