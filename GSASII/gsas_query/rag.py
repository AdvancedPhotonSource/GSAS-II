"""
RAG engine: query ChromaDB, retrieve relevant chunks, generate answer.

Backend selection (in priority order):
  1. If ``LLM_BACKEND`` env var is set explicitly, it is always honoured.
  2. If ``LLM_BACKEND`` is not set and ``llama_cpp`` (llama-cpp-python) is
     importable, the llama_cpp backend is selected automatically.
  3. Otherwise the default is ``ollama``.

Supported LLM_BACKEND values:
  "ollama"     — local Ollama server, fully on-premises
  "anthropic"  — Anthropic Claude API (requires ANTHROPIC_API_KEY)
  "llama_cpp"  — in-process llama-cpp-python (requires LLAMA_CPP_MODEL path)
  "retrieval"  — no LLM; returns raw matched chunks (offline / testing)
"""

import os
from functools import lru_cache

import chromadb

from ._embed import get_embedding_function
from ._paths import get_chroma_path

COLLECTION_NAME = "gsasii_docs"
TOP_K = 6
FETCH_K = 30  # over-fetch then deduplicate by URL to maximise source diversity

SYSTEM_PROMPT = """\
You are an expert assistant for GSAS-II (General Structure Analysis System-II), \
crystallographic analysis software developed at Argonne National Laboratory.

You answer questions using the GSAS-II tutorials, help manual, and documentation \
provided as context. Your users are crystallographers, materials scientists, and \
physicists who need precise, actionable answers.

Guidelines:
- Be specific and technical. Use proper crystallographic terminology.
- When a question involves a step-by-step procedure, preserve the numbered steps.
- Each context section is numbered [1], [2], etc. Cite sources inline as you write \
  by inserting the number in brackets (e.g. "Open the Phase tab [1] and select..."). \
  Place the citation immediately after the sentence or clause it supports. Do NOT add \
  a separate references section at the end. Do NOT reproduce the raw "[Source: ...]" \
  labels from the context — use only the numeric [N] markers.
- If the provided context does not contain enough information to answer, say so clearly \
  and suggest which tutorial might cover the topic.
- Do not fabricate parameter names, menu paths, or file formats.
"""


@lru_cache(maxsize=1)
def _get_collection() -> chromadb.Collection:
    client = chromadb.PersistentClient(path=get_chroma_path())
    return client.get_or_create_collection(COLLECTION_NAME)


def _retrieve(question: str) -> tuple[str, list[dict], dict[str, dict], list[str]]:
    collection = _get_collection()

    if collection.count() == 0:
        return "", [], {}, []

    ef = get_embedding_function()
    embedding = ef([question])[0]
    fetch_k = min(FETCH_K, collection.count())
    results = collection.query(
        query_embeddings=[embedding],
        n_results=fetch_k,
        include=["documents", "metadatas", "distances"],
    )

    # URL-diversity dedup: keep only the best-scoring chunk per unique URL,
    # then take top TOP_K. Prevents a single page (e.g. tutorials.html) from
    # filling multiple result slots with similar content.
    raw_docs = results["documents"][0]
    raw_metas = results["metadatas"][0]
    raw_dists = results["distances"][0]

    seen_urls: dict[str, int] = {}  # url -> index of best (lowest dist) chunk
    ordered_indices: list[int] = []
    for idx, (meta, dist) in enumerate(zip(raw_metas, raw_dists)):
        url = meta.get("url", "")
        if url not in seen_urls:
            seen_urls[url] = idx
            ordered_indices.append(idx)
        # lower distance = better; first occurrence is always the best since
        # ChromaDB returns results sorted by ascending distance
    top_indices = ordered_indices[:TOP_K]

    # Normalise distances within the selected set
    selected_dists = [raw_dists[i] for i in top_indices]
    min_d = min(selected_dists) if selected_dists else 0.0
    max_d = max(selected_dists) if selected_dists else 1.0
    span = (max_d - min_d) or 1.0

    context_parts = []
    citations: dict[str, dict] = {}
    sources = []
    seen_sources: set = set()
    chunk_texts: list[str] = []

    for i, idx in enumerate(top_indices, start=1):
        doc = raw_docs[idx]
        meta = raw_metas[idx]
        dist = raw_dists[idx]
        context_parts.append(
            f"[{i}] [Source: {meta['title']} | Section: {meta['section']}]\n{doc}"
        )
        relevance = round((max_d - dist) / span, 3)   # 1.0 = best, 0.0 = weakest
        citations[str(i)] = {
            "title": meta["title"],
            "section": meta["section"],
            "url": meta["url"],
            "relevance": relevance,
        }
        chunk_texts.append(doc)
        source_key = (meta["url"], meta["section"])
        if source_key not in seen_sources:
            seen_sources.add(source_key)
            sources.append({
                "title": meta["title"],
                "section": meta["section"],
                "url": meta["url"],
                "category": meta.get("category", ""),
                "relevance": relevance,
            })

    context = "\n\n---\n\n".join(context_parts)
    return context, sources, citations, chunk_texts


def _build_messages(question: str, context: str, history: list[dict]) -> list[dict]:
    messages = []
    for turn in history[-6:]:
        if turn.get("role") in {"user", "assistant"}:
            messages.append({"role": turn["role"], "content": turn["content"]})

    user_content = (
        f"Context from GSAS-II documentation:\n\n{context}\n\n"
        f"---\n\nQuestion: {question}"
        if context
        else question
    )
    messages.append({"role": "user", "content": user_content})
    return messages


def _effective_backend() -> str:
    """Return the backend to use.

    Priority:
      1. ``LLM_BACKEND`` env var when explicitly set.
      2. ``llama_cpp`` when llama-cpp-python is importable and LLM_BACKEND unset.
      3. ``ollama`` as the final default.
    """
    env_backend = os.environ.get("LLM_BACKEND", "").strip().lower()
    if env_backend:
        return env_backend
    try:
        import llama_cpp  # noqa: F401
        return "llama_cpp"
    except ImportError:
        pass
    return "ollama"


def _answer_anthropic(messages: list[dict]) -> str:
    import anthropic
    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
    response = client.messages.create(
        model=os.environ.get("ANTHROPIC_MODEL", "claude-sonnet-4-6"),
        max_tokens=1500,
        system=SYSTEM_PROMPT,
        messages=messages,
    )
    return response.content[0].text


@lru_cache(maxsize=1)
def _get_llama(model_path: str, n_ctx: int):
    from llama_cpp import Llama
    return Llama(model_path=model_path, n_ctx=n_ctx, verbose=False)


def _answer_llama_cpp(messages: list[dict]) -> str:
    model_path = os.environ.get("LLAMA_CPP_MODEL", "").strip()
    if not model_path:
        raise RuntimeError(
            "LLAMA_CPP_MODEL is not set. "
            "Set it to the path of a GGUF model file, e.g. "
            "LLAMA_CPP_MODEL=/path/to/model.gguf"
        )
    n_ctx = int(os.environ.get("LLAMA_CPP_N_CTX", "4096"))
    # Cap at 800 tokens for small models to prevent repetition loops.
    max_tokens = int(os.environ.get("LLAMA_CPP_MAX_TOKENS", "800"))
    llm = _get_llama(model_path, n_ctx)
    full_messages = [{"role": "system", "content": SYSTEM_PROMPT}] + messages
    response = llm.create_chat_completion(messages=full_messages, max_tokens=max_tokens)
    return response["choices"][0]["message"]["content"]


def _ollama_url() -> str:
    return os.environ.get("OLLAMA_URL", "http://localhost:11434")


def _installed_ollama_models() -> list[str]:
    import httpx
    resp = httpx.get(f"{_ollama_url()}/api/tags", timeout=5)
    resp.raise_for_status()
    data = resp.json() or {}
    return [m.get("name", "") for m in data.get("models", []) if m.get("name")]


def _choose_ollama_model() -> str:
    preferred = os.environ.get("OLLAMA_MODEL", "").strip()
    models = _installed_ollama_models()

    if not models:
        raise RuntimeError(
            "No Ollama models are installed. Run e.g. "
            "`ollama pull llama3.1:8b` or `ollama pull qwen2.5:3b`."
        )

    if preferred:
        if preferred in models:
            return preferred
        raise RuntimeError(
            f"OLLAMA_MODEL='{preferred}' is not installed. "
            f"Available models: {', '.join(models)}"
        )

    for candidate in ("llama3.1:8b", "llama3", "qwen2.5:3b"):
        if candidate in models:
            return candidate

    return models[0]


def _answer_ollama(messages: list[dict]) -> str:
    import httpx
    ollama_url = _ollama_url()
    ollama_model = _choose_ollama_model()

    full_messages = [{"role": "system", "content": SYSTEM_PROMPT}] + messages
    try:
        resp = httpx.post(
            f"{ollama_url}/api/chat",
            json={"model": ollama_model, "messages": full_messages, "stream": False},
            timeout=120,
        )
        resp.raise_for_status()
    except httpx.HTTPStatusError as e:
        detail = ""
        try:
            detail = f" | Response: {e.response.text}"
        except Exception:
            pass
        raise RuntimeError(f"Ollama API error: {e}{detail}") from e
    return resp.json()["message"]["content"]


def answer_question(question: str, history: list[dict]) -> dict:
    """Main entry point: returns {answer, sources, citations, backend}."""
    if not question.strip():
        return {"answer": "Please enter a question.", "sources": [], "citations": {}, "backend": ""}

    context, sources, citations, chunk_texts = _retrieve(question)

    if not context:
        return {
            "answer": (
                "The knowledge base is empty. Run `gsas2-query --setup` "
                "to index the GSAS-II documentation."
            ),
            "sources": [],
            "citations": {},
            "backend": "",
        }

    backend = _effective_backend()

    if backend == "retrieval":
        answer = (
            "Most relevant sections (no LLM synthesis — install llama-cpp-python, "
            "Ollama, or set LLM_BACKEND=anthropic for generated answers):\n\n" + context
        )
        return {"answer": answer, "sources": sources, "citations": citations, "backend": backend}

    messages = _build_messages(question, context, history)

    if backend == "llama_cpp":
        answer = _answer_llama_cpp(messages)
    elif backend == "ollama":
        answer = _answer_ollama(messages)
    else:
        answer = _answer_anthropic(messages)

    answer = _inject_citations(answer, chunk_texts)
    answer, citations = _renumber_by_appearance(answer, citations)
    return {"answer": answer, "sources": sources, "citations": citations, "backend": backend}


def _renumber_by_appearance(answer: str, citations: dict) -> tuple[str, dict]:
    """Renumber [N] markers so citations are sequential in order of first appearance.

    The injector assigns numbers by relevance rank, so [4] may appear before [2]
    in the text. This makes citations appear in order: first cited = [1], etc.
    Sources not cited inline are appended at the end of the new citations dict
    (for the bibliography) but get no inline marker.
    """
    import re

    # Walk the text to get first-appearance order of each cited number
    appeared: list[str] = []
    seen: set[str] = set()
    for m in re.finditer(r'\[(\d+)\]', answer):
        k = m.group(1)
        if k in citations and k not in seen:
            appeared.append(k)
            seen.add(k)

    # Append any retrieved-but-uncited chunks (for full bibliography)
    for k in sorted(citations, key=int):
        if k not in seen:
            appeared.append(k)

    # old key → new sequential key
    old_to_new = {old: str(i + 1) for i, old in enumerate(appeared)}

    # Rewrite [N] markers in the answer text
    new_answer = re.sub(
        r'\[(\d+)\]',
        lambda m: f'[{old_to_new.get(m.group(1), m.group(1))}]',
        answer,
    )

    # Rebuild citations dict in new order
    new_citations = {old_to_new[k]: citations[k] for k in appeared}
    return new_answer, new_citations


def _inject_citations(answer: str, chunk_texts: list[str]) -> str:
    """Inject [N] citation markers inline using keyword overlap with context chunks.

    Small models (3B) rarely cite more than one source. This post-processor
    scans each substantive line, finds the best-matching context chunk by
    content-word overlap, and appends [N] if not already present. Requires at
    least 4 shared content words (≥4 chars) to avoid spurious matches.
    """
    import re

    # Pre-compute content-word sets for each chunk (1-indexed)
    _STOP = {"gsas", "gsasii", "that", "with", "this", "from", "will", "have",
             "your", "then", "they", "each", "also", "been", "used", "when"}

    def content_words(text: str) -> set[str]:
        return {w for w in re.findall(r'\b[a-zA-Z]{4,}\b', text.lower())
                if w not in _STOP}

    chunk_word_sets = [content_words(t) for t in chunk_texts]

    def best_chunk_for(line: str) -> int | None:
        lwords = content_words(line)
        if len(lwords) < 5:
            return None
        best_score = 3  # minimum overlap threshold
        best_n = None
        for n, cwords in enumerate(chunk_word_sets, start=1):
            score = len(lwords & cwords)
            if score > best_score:
                best_score = score
                best_n = n
        return best_n

    lines = answer.split('\n')
    result = []
    for line in lines:
        stripped = line.strip()
        # Skip: already cited, blank, bullet markers alone, or very short
        if (not stripped
                or re.search(r'\[\d+\]', stripped)
                or len(stripped) < 40):
            result.append(line)
            continue
        n = best_chunk_for(stripped)
        if n is not None:
            result.append(line.rstrip() + f' [{n}]')
        else:
            result.append(line)
    return '\n'.join(result)
