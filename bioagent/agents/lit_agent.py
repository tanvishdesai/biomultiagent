"""LitAgent — PubMed literature retrieval via Bio.Entrez."""

from __future__ import annotations

import os
import re
import time
from typing import Any, Dict, List, Optional

from Bio import Entrez

ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "bio.research@example.com")
ENTREZ_TOOL  = "biomultiagent"
RATE_LIMIT_SLEEP = 0.34

Entrez.email = ENTREZ_EMAIL
Entrez.tool  = ENTREZ_TOOL

_NL_BOILERPLATE = re.compile(
    r"\b(?:find|search|get|list|show|recent|latest|pubmed|literature|"
    r"papers?|articles?|studies|me|the|a|an|\d+)\b",
    re.I,
)


def _strip_nl_boilerplate(text: str) -> str:
    cleaned = _NL_BOILERPLATE.sub(" ", text)
    return re.sub(r"\s+", " ", cleaned).strip(" ,.")


def _extract_pubmed_query(query: str) -> str:
    """Turn natural-language queries into PubMed search terms."""
    q = query.strip()
    if not q:
        return q

    lit_patterns = [
        r"find\s+(?:\d+\s+)?(?:recent\s+)?(?:pubmed\s+)?(?:papers?\s+)?"
        r"(?:on|about|regarding|for)\s+(.+)",
        r"(?:papers?\s+)(?:on|about|regarding|for)\s+(.+)",
        r"search\s+(?:pubmed\s+)?(?:for\s+)?(.+)",
    ]
    for pat in lit_patterns:
        m = re.search(pat, q, re.I)
        if m:
            return _strip_nl_boilerplate(m.group(1))

    m = re.search(r"\b(?:on|about|regarding|for)\s+(.+)$", q, re.I)
    if m:
        return _strip_nl_boilerplate(m.group(1))

    return _strip_nl_boilerplate(q)


def search_pubmed(
    query: str,
    max_results: int = 5,
    sort_recent: bool = False,
) -> tuple[List[str], Optional[str]]:
    """
    Return PubMed IDs for the query.

    Returns (pmids, error_message). error_message is set on network/API failure.
    """
    try:
        sort_mode = "pub date" if sort_recent else "relevance"
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort=sort_mode,
        )
        record = Entrez.read(handle)
        handle.close()
        pmids = record.get("IdList", [])
        return pmids, None
    except Exception as exc:
        return [], f"PubMed search failed: {exc}"


def _abstract_text(article: Dict) -> str:
    """Concatenate all AbstractText elements (handles structured abstracts)."""
    abstract = article.get("Abstract", {})
    texts = abstract.get("AbstractText", [])
    if not texts:
        return "No abstract available."
    parts = []
    for item in texts:
        if isinstance(item, str):
            parts.append(item)
        else:
            label = item.attributes.get("Label", "") if hasattr(item, "attributes") else ""
            parts.append(f"{label}: {item}" if label else str(item))
    return " ".join(parts)[:600]


def fetch_abstracts(pmids: List[str]) -> tuple[List[Dict[str, str]], Optional[str]]:
    """Fetch title + abstract for a list of PMIDs."""
    if not pmids:
        return [], None
    time.sleep(RATE_LIMIT_SLEEP)
    try:
        handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="xml", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as exc:
        return [], f"PubMed fetch failed: {exc}"

    results = []
    for article in records.get("PubmedArticle", []):
        med = article.get("MedlineCitation", {})
        art = med.get("Article", {})
        pmid = str(med.get("PMID", "?"))
        title = str(art.get("ArticleTitle", "No title")).rstrip(".")
        abstract = _abstract_text(art)
        journal = art.get("Journal", {})
        journal_title = str(journal.get("Title", ""))
        year = ""
        pub_date = journal.get("JournalIssue", {}).get("PubDate", {})
        if "Year" in pub_date:
            year = str(pub_date["Year"])

        results.append({
            "pmid":     pmid,
            "title":    title,
            "abstract": abstract,
            "year":     year,
            "journal":  journal_title,
            "url":      f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })
    return results, None


def run(
    task: str = "",
    query: str = "",
    max_results: int = 5,
    **kwargs,
) -> Dict[str, Any]:
    raw_query = (query or task).strip()
    if not raw_query:
        return {"result": "No search query provided.", "details": {}}

    search_query = _extract_pubmed_query(raw_query)
    if not search_query:
        search_query = raw_query

    sort_recent = bool(re.search(r"\b(?:recent|latest|new)\b", raw_query, re.I))
    pmids, search_err = search_pubmed(search_query, max_results, sort_recent=sort_recent)
    if search_err:
        return {
            "result": search_err,
            "error":  search_err,
            "details": {"query": search_query, "raw_query": raw_query},
        }
    if not pmids:
        return {
            "result": f"No PubMed hits for: '{search_query}'",
            "details": {"query": search_query, "raw_query": raw_query},
        }

    papers, fetch_err = fetch_abstracts(pmids)
    if fetch_err:
        return {
            "result": fetch_err,
            "error":  fetch_err,
            "details": {"pmids": pmids},
        }
    if not papers:
        return {
            "result": f"Retrieved {len(pmids)} PMIDs but no parseable abstracts.",
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
        if p.get("abstract") and p["abstract"] != "No abstract available.":
            summary_lines.append(f"    {p['abstract'][:200]}…")
        citations.append(f"PMID:{p['pmid']} — {p['title']}")

    result_text = f"Found {len(papers)} papers on '{search_query}':\n" + "\n".join(summary_lines)
    return {
        "result":    result_text,
        "papers":    papers,
        "citations": citations,
        "details":   {
            "query": search_query,
            "raw_query": raw_query,
            "n_results": len(papers),
        },
    }
