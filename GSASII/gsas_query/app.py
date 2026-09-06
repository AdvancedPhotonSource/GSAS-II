"""
FastAPI web server for the GSAS-II documentation chatbot.

  GET  /         -> chat UI (static/index.html)
  GET  /health   -> {"status": "ok"}
  GET  /stats    -> {"chunks_indexed": N, "llm_backend": "..."}
  POST /chat     -> RAG query, returns {answer, sources, citations}
  POST /ingest   -> trigger re-indexing (requires X-Admin-Key header)
"""

import os
import subprocess
import sys
import time
from typing import Optional

from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field

from ._paths import get_static_dir

STATIC_DIR = get_static_dir()

app = FastAPI(
    title="GSAS-II Documentation Assistant",
    description="RAG chatbot over GSAS-II tutorials and help documentation",
    version="0.2.0",
)

origins = os.environ.get("ALLOWED_ORIGINS", "*").split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


@app.exception_handler(Exception)
async def generic_exception_handler(request: Request, exc: Exception):
    return JSONResponse(status_code=500, content={"detail": str(exc)})


# --------------------------------------------------------------------------- #
# Request / response models                                                    #
# --------------------------------------------------------------------------- #

class HistoryTurn(BaseModel):
    role: str
    content: str


class ChatRequest(BaseModel):
    message: str = Field(..., max_length=2000)
    history: list[HistoryTurn] = Field(default_factory=list, max_length=20)


class Source(BaseModel):
    title: str
    section: str
    url: str
    category: str
    relevance: float


class Citation(BaseModel):
    title: str
    section: str
    url: str
    relevance: float = 0.0


class ChatResponse(BaseModel):
    answer: str
    sources: list[Source]
    citations: dict[str, Citation]
    elapsed_ms: int


# --------------------------------------------------------------------------- #
# Simple rate limiting (in-memory, per IP)                                     #
# --------------------------------------------------------------------------- #

_rate_store: dict[str, list[float]] = {}
RATE_LIMIT = int(os.environ.get("RATE_LIMIT_RPM", "30"))


def check_rate_limit(request: Request):
    if RATE_LIMIT <= 0:
        return
    ip = request.client.host if request.client else "unknown"
    now = time.time()
    timestamps = [t for t in _rate_store.get(ip, []) if now - t < 60]
    if len(timestamps) >= RATE_LIMIT:
        raise HTTPException(status_code=429, detail="Rate limit exceeded. Please slow down.")
    timestamps.append(now)
    _rate_store[ip] = timestamps


# --------------------------------------------------------------------------- #
# Routes                                                                       #
# --------------------------------------------------------------------------- #

@app.get("/", include_in_schema=False)
async def serve_ui():
    return FileResponse(str(STATIC_DIR / "index.html"))


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.get("/stats")
async def stats():
    from .rag import _effective_backend, _get_collection
    try:
        col = _get_collection()
        count = col.count()
    except Exception:
        count = -1
    return {"chunks_indexed": count, "llm_backend": _effective_backend()}


@app.post("/chat", response_model=ChatResponse)
async def chat(
    request: ChatRequest,
    http_request: Request,
    _: None = Depends(check_rate_limit),
):
    from .rag import answer_question

    t0 = time.monotonic()
    history = [h.model_dump() for h in request.history]
    result = answer_question(request.message, history)
    elapsed = int((time.monotonic() - t0) * 1000)

    return ChatResponse(
        answer=result["answer"],
        sources=result["sources"],
        citations=result.get("citations", {}),
        elapsed_ms=elapsed,
    )


@app.post("/ingest")
async def trigger_ingest(http_request: Request):
    admin_key = os.environ.get("ADMIN_KEY", "")
    if admin_key and http_request.headers.get("X-Admin-Key") != admin_key:
        raise HTTPException(status_code=403, detail="Invalid admin key.")
    subprocess.Popen([sys.executable, "-m", "gsas_query.ingest"])
    return {"status": "ingestion started in background"}
