"""LitAgent — PubMed literature retrieval and abstract summarisation."""

from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional


ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENTREZ_EMAIL = "bio.research@example.com"   # required by NCBI; update to your email
RATE_LIMIT_SLEEP = 0.35   # respect NCBI's 3 req/s limit for unauthenticated requests


def _get(url: str, retries: int = 3) -> Optional[bytes]:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                return resp.read()
        except urllib.error.URLError as e:
            if attempt < retries - 1:
                time.sleep(1.0)
            else:
                return None
    return None


def search_pubmed(query: str, max_results: int = 5) -> List[str]:
    """Return a list of PubMed IDs for the query."""
    params = urllib.parse.urlencode({
        "db":      "pubmed",
        "term":    query,
        "retmax":  max_results,
        "rettype": "json",
        "tool":    "biomultiagent",
        "email":   ENTREZ_EMAIL,
    })
    url  = f"{ENTREZ_BASE}/esearch.fcgi?{params}"
    data = _get(url)
    if not data:
        return []
    try:
        root = ET.fromstring(data)
        return [id_el.text for id_el in root.findall(".//Id")]
    except ET.ParseError:
        return []


def fetch_abstracts(pmids: List[str]) -> List[Dict[str, str]]:
    """Fetch title + abstract for a list of PMIDs."""
    if not pmids:
        return []
    time.sleep(RATE_LIMIT_SLEEP)
    params = urllib.parse.urlencode({
        "db":      "pubmed",
        "id":      ",".join(pmids),
        "rettype": "abstract",
        "retmode": "xml",
        "tool":    "biomultiagent",
        "email":   ENTREZ_EMAIL,
    })
    url  = f"{ENTREZ_BASE}/efetch.fcgi?{params}"
    data = _get(url)
    if not data:
        return []

    results = []
    try:
        root = ET.fromstring(data)
        for article in root.findall(".//PubmedArticle"):
            pmid_el    = article.find(".//PMID")
            title_el   = article.find(".//ArticleTitle")
            abstract_el = article.find(".//AbstractText")
            year_el    = article.find(".//PubDate/Year")
            journal_el = article.find(".//Journal/Title")

            pmid      = pmid_el.text if pmid_el is not None else "?"
            title     = title_el.text if title_el is not None else "No title"
            abstract  = abstract_el.text if abstract_el is not None else "No abstract available."
            year      = year_el.text if year_el is not None else ""
            journal   = journal_el.text if journal_el is not None else ""

            if title:
                title = title.rstrip(".")

            results.append({
                "pmid":     pmid,
                "title":    title,
                "abstract": (abstract or "")[:600],
                "year":     year,
                "journal":  journal,
                "url":      f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })
    except ET.ParseError:
        pass

    return results


def run(
    task: str,
    query: str = "",
    max_results: int = 5,
    **kwargs,
) -> Dict[str, Any]:
    search_query = query or task
    if not search_query:
        return {"result": "No search query provided.", "details": {}}

    pmids   = search_pubmed(search_query, max_results)
    if not pmids:
        return {
            "result": f"No PubMed results found for: '{search_query}'",
            "details": {"query": search_query},
        }

    papers = fetch_abstracts(pmids)
    if not papers:
        return {
            "result": f"Retrieved {len(pmids)} PMIDs but could not fetch abstracts.",
            "details": {"pmids": pmids},
        }

    summary_lines = []
    citations = []
    for i, p in enumerate(papers, 1):
        year_str = f" ({p['year']})" if p.get("year") else ""
        summary_lines.append(
            f"[{i}] {p['title']}{year_str} — {p.get('journal', '')} "
            f"[PMID {p['pmid']}]"
        )
        citations.append(f"PMID:{p['pmid']} — {p['title']}")

    result_text = f"Found {len(papers)} papers on '{search_query}':\n" + "\n".join(summary_lines)
    return {
        "result":    result_text,
        "papers":    papers,
        "citations": citations,
        "details":   {"query": search_query, "n_results": len(papers)},
    }
