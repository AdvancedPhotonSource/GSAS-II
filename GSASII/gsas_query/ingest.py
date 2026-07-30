"""
Ingestion pipeline: fetch GSAS-II tutorials (HTML + PDFs),
chunk by section, embed with sentence-transformers, store in ChromaDB.
"""

import argparse
import hashlib
import io
import multiprocessing
import re
import sys
import time

from ._paths import get_chroma_path

COLLECTION_NAME = "gsasii_docs"
MAX_CHUNK_CHARS = 1200
OVERLAP_CHARS = 150
REQUEST_DELAY = 0.5


def get_collection(reset: bool = False):
    import chromadb
    client = chromadb.PersistentClient(path=get_chroma_path())
    if reset:
        try:
            client.delete_collection(COLLECTION_NAME)
            print("Dropped existing collection.")
        except Exception:
            pass
    return client.get_or_create_collection(
        name=COLLECTION_NAME,
        metadata={"hnsw:space": "cosine"},
    )


def chunk_text(text: str, max_chars: int = MAX_CHUNK_CHARS, overlap: int = OVERLAP_CHARS) -> list[str]:
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) <= max_chars:
        return [text] if text else []

    chunks = []
    start = 0
    while start < len(text):
        end = min(start + max_chars, len(text))
        if end < len(text):
            boundary = max(
                text.rfind(". ", start, end),
                text.rfind(".\n", start, end),
                text.rfind("! ", start, end),
                text.rfind("? ", start, end),
            )
            if boundary > start + max_chars // 2:
                end = boundary + 1
        chunk = text[start:end].strip()
        if chunk:
            chunks.append(chunk)
        if end >= len(text):
            break
        start = end - overlap
    return chunks


def extract_html_sections(html: str, source_title: str) -> list[dict]:
    from bs4 import BeautifulSoup, Tag
    soup = BeautifulSoup(html, "html.parser")

    for tag in soup(["script", "style", "nav", "footer", "img"]):
        tag.decompose()

    body = soup.find("body") or soup
    sections = []
    current_heading = source_title
    current_text_parts: list[str] = []
    heading_tags = {"h1", "h2", "h3", "h4"}

    def flush():
        text = " ".join(current_text_parts).strip()
        if text and len(text) > 80:
            sections.append({"heading": current_heading, "text": text})

    for elem in body.descendants:
        if not isinstance(elem, Tag):
            continue
        if elem.name in heading_tags:
            flush()
            current_heading = elem.get_text(separator=" ", strip=True)
            current_text_parts = []
        elif elem.name in {"p", "li", "td", "th", "pre", "blockquote", "dd", "dt"}:
            txt = elem.get_text(separator=" ", strip=True)
            if txt:
                current_text_parts.append(txt)

    flush()
    return sections


def ingest_html_source(source: dict, collection, model):
    import requests
    url = source["url"]
    title = source["title"]
    category = source["category"]

    print(f"  Fetching: {title} ({url})")
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        print(f"  ERROR fetching {url}: {e}")
        return 0

    sections = extract_html_sections(resp.text, title)
    # Use a dict keyed by doc_id to deduplicate identical chunks within the source.
    records: dict[str, tuple] = {}

    for section in sections:
        chunks = chunk_text(section["text"])
        for chunk in chunks:
            doc_id = hashlib.md5(f"{url}|{section['heading']}|{chunk}".encode()).hexdigest()
            if doc_id not in records:
                records[doc_id] = (chunk, model([chunk])[0], {
                    "url": url,
                    "title": title,
                    "section": section["heading"],
                    "category": category,
                    "source_type": "html",
                })

    if records:
        ids = list(records)
        docs = [v[0] for v in records.values()]
        embeddings = [v[1] for v in records.values()]
        metadatas = [v[2] for v in records.values()]
        collection.upsert(ids=ids, documents=docs, embeddings=embeddings, metadatas=metadatas)
        print(f"    -> {len(ids)} chunks stored")
        time.sleep(REQUEST_DELAY)
        return len(ids)

    time.sleep(REQUEST_DELAY)
    return 0


def ingest_pdf_source(source: dict, collection, model):
    try:
        from pypdf import PdfReader
    except ImportError:
        print("  pypdf not installed, skipping PDF ingestion.")
        return 0

    import requests
    url = source["url"]
    title = source["title"]
    category = source["category"]

    print(f"  Fetching PDF: {title}")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except Exception as e:
        print(f"  ERROR fetching PDF {url}: {e}")
        return 0

    try:
        reader = PdfReader(io.BytesIO(resp.content))
    except Exception as e:
        print(f"  ERROR parsing PDF: {e}")
        return 0

    ids, docs, embeddings, metadatas = [], [], [], []
    total_pages = len(reader.pages)
    print(f"    {total_pages} pages")

    for page_start in range(0, total_pages, 3):
        page_end = min(page_start + 3, total_pages)
        combined_text = ""
        for p in range(page_start, page_end):
            combined_text += (reader.pages[p].extract_text() or "") + "\n"

        section_label = f"Pages {page_start + 1}-{page_end}"
        for i, chunk in enumerate(chunk_text(combined_text)):
            doc_id = hashlib.md5(f"{url}|{section_label}|{i}".encode()).hexdigest()
            embedding = model([chunk])[0]
            ids.append(doc_id)
            docs.append(chunk)
            embeddings.append(embedding)
            metadatas.append({
                "url": url,
                "title": title,
                "section": section_label,
                "category": category,
                "source_type": "pdf",
            })

    if ids:
        collection.upsert(ids=ids, documents=docs, embeddings=embeddings, metadatas=metadatas)
        print(f"    -> {len(ids)} chunks stored")

    return len(ids)


def main():
    parser = argparse.ArgumentParser(description="Ingest GSAS-II docs into ChromaDB")
    parser.add_argument("--html-only", action="store_true", help="Skip PDF ingestion")
    parser.add_argument("--reset", action="store_true", help="Drop and rebuild collection")
    parser.add_argument("--book", action="store_true",
                        help="Include Powder Diffraction Crystallography book (185 HTML pages)")
    parser.add_argument("--manual", action="store_true",
                        help="Include Programmers' Manual (24 HTML chapters)")
    args = parser.parse_args()

    from chromadb.utils.embedding_functions import DefaultEmbeddingFunction
    from .sources import get_tutorial_sources
    from .sources import HOME_SOURCES, HELP_SOURCES, TUTORIAL_SOURCES
    from .sources import READTHEDOCS_SOURCES, BOOK_HTML_SOURCES

    print("Loading embedding model (ONNX all-MiniLM-L6-v2)...")
    model = DefaultEmbeddingFunction()

    print(f"ChromaDB path: {get_chroma_path()}")
    collection = get_collection(reset=args.reset)

    print("Fetching tutorial list from GSAS-II repository…")
    tutorial_sources = get_tutorial_sources()
    print(f"  {len(tutorial_sources)} tutorials found. ({len(TUTORIAL_SOURCES)} in hard-coded list)")

    total_chunks = 0

    # add a notes to distinguish information sources
    # W: Web; T: Tutorial; M: Manual; (Help pages & Book already tagged.)
    for i in HOME_SOURCES:
        if 'title' in i:
            i['title'] += ' (W)'

    for i in tutorial_sources:
        if 'title' in i:
            i['title'] += ' (T)'
    
    html_sources = [HOME_SOURCES, HELP_SOURCES, tutorial_sources]
    if args.manual:
        for i in READTHEDOCS_SOURCES:
            if 'title' in i:
                i['title'] += ' (M)'

        print(f"  Adding Programmers' manual ({len(READTHEDOCS_SOURCES)} HTML pages)")
        html_sources = html_sources + [READTHEDOCS_SOURCES]
    if args.book:
        html_sources = html_sources + [BOOK_HTML_SOURCES]

    print("\n=== Ingesting HTML pages ===")
    for pagelist in html_sources:
        for source in pagelist:
            if type(source) is str:
                print(f"\n*** processing {source}")
                continue
            total_chunks += ingest_html_source(source, collection, model)
    if not args.html_only:
        print("\n=== Ingesting PDFs ===")
        from .sources import PDF_SOURCES
        for source in PDF_SOURCES:
            total_chunks += ingest_pdf_source(source, collection, model)

    print(f"\nDone. Total chunks in collection: {collection.count()}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
