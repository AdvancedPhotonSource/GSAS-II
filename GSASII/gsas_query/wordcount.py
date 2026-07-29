"""
word count pipeline: fetch GSAS-II documents and count words
"""

#import argparse
#import hashlib
#import io
#import multiprocessing
#import re
#import sys
#import time
import requests
from bs4 import BeautifulSoup

def wordcount_html_source(source: dict):
    url = source["url"]
    title = source["title"]
    category = source["category"]

    print(f"  Fetching: {title} ({url})")
    try:
        response = requests.get(url, timeout=10)
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Extracts text only from the body, ignoring raw HTML code
        text = soup.body.get_text() if soup.body else soup.get_text()
        
        # Split text by whitespace to count individual words
        word_count = len(text.split())
        print(f"{url}: {word_count} words")
        return word_count
    except Exception as e:
        print(f"Could not process {url}: {e}")
        return 0

def main():
    from .sources import get_tutorial_sources
    from .sources import HOME_SOURCES, HELP_SOURCES, TUTORIAL_SOURCES
    from .sources import READTHEDOCS_SOURCES, BOOK_HTML_SOURCES

    print("Fetching tutorial list from GSAS-II repository…")
    tutorial_sources = get_tutorial_sources()
    print(f"  {len(tutorial_sources)} tutorials found. ({len(TUTORIAL_SOURCES)} in hard-coded list)")

    html_sources = [HOME_SOURCES, HELP_SOURCES, tutorial_sources]
    html_sources = html_sources + [READTHEDOCS_SOURCES]
    html_sources = html_sources + [BOOK_HTML_SOURCES]

    print("\n=== Scanning HTML pages ===")
    total_words = 0
    wordsDict = {}
    filesDict = {}
    for pagelist in html_sources:
        words = 0
        lbl = '?'
        i = 0
        for source in pagelist:
            if type(source) is str:
                print(f"\n*** processing {source}")
                lbl = source
                continue
            i += 1
            
            words += wordcount_html_source(source)
            total_words += wordcount_html_source(source)
        print(f"words: {words} files: {i}")
        wordsDict[lbl] = words
        filesDict[lbl] = i
            #if i >= 5: break  # DEBUG

    print(f"\nDone. Total words: {total_words}")
    for lbl in wordsDict:
        print(lbl,'words',wordsDict[lbl],'files',filesDict[lbl])


if __name__ == "__main__":
#    multiprocessing.freeze_support()
    main()
