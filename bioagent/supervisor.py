"""
BioMultiAgent Supervisor
========================
LangGraph-orchestrated multi-agent pipeline for bioinformatics compound queries.

Graph nodes:
  classify → route → [seq/align/annot/phylo/lit] → fuse → END

Falls back to keyword-based intent classification when no LLM is configured,
so the system works fully offline for demos.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import uuid
from typing import Any, Dict, List

from langgraph.graph import END, StateGraph

from bioagent.agents import AGENTS
from bioagent.memory.chroma_store import get_memory
from bioagent.state import BioAgentState


# ─── Intent classification ────────────────────────────────────────────────────

INTENT_KEYWORDS: Dict[str, List[str]] = {
    "seq":        ["translate", "gc content", "gc%", "motif", "sequence analysis", "reverse complement"],
    "align":      ["align", "alignment", "msa", "multiple sequence", "pairwise"],
    "annot":      ["orf", "annotate", "annotation", "open reading frame", "codon"],
    "phylo":      ["phylo", "tree", "phylogen", "neighbor", "evolution", "dendrogram"],
    "literature": ["pubmed", "paper", "literature", "article", "recent", "find papers"],
}


def _keyword_classify(query: str) -> tuple[str, list[str]]:
    q = query.lower()
    matched = [k for k, kws in INTENT_KEYWORDS.items() if any(kw in q for kw in kws)]
    # Compound = 2+ intents, or query contains "and", or all/multiple explicit tasks
    if len(matched) >= 2 or " and " in q or "compound" in q:
        sub = [i for i in ["seq", "annot", "align", "phylo", "literature"] if i in matched]
        return "compound", sub or list(matched)
    if matched:
        return matched[0], [matched[0]]
    return "seq", ["seq"]   # default to sequence analysis


def _llm_classify(query: str) -> tuple[str, list[str]]:
    """Try to use Ollama for more accurate intent classification."""
    prompt = (
        'Classify this bioinformatics query. Respond with ONLY valid JSON:\n'
        '{"intent": "<seq|align|annot|phylo|literature|compound>", '
        '"sub_tasks": ["<intent1>", ...]}\n\n'
        f'Query: "{query}"'
    )
    try:
        r = subprocess.run(
            ["ollama", "run", "mistral", prompt],
            capture_output=True, text=True, timeout=20,
        )
        if r.returncode == 0:
            m = re.search(r'\{.*?\}', r.stdout, re.DOTALL)
            if m:
                data = json.loads(m.group())
                intent    = data.get("intent", "seq")
                sub_tasks = data.get("sub_tasks", [intent])
                return intent, sub_tasks
    except Exception:
        pass
    return _keyword_classify(query)


def _extract_sequences(query: str) -> List[str]:
    """Extract DNA/RNA sequences from query text (heuristic)."""
    return re.findall(r'\b([ACGTUN]{15,})\b', query.upper())


# ─── LangGraph nodes ──────────────────────────────────────────────────────────

def classify_node(state: BioAgentState) -> BioAgentState:
    """Classify query intent and extract sequences."""
    query     = state["query"]
    intent, sub_tasks = _llm_classify(query)
    sequences = _extract_sequences(query)

    # Load session memory context
    memory   = get_memory()
    session  = state.get("session_id") or str(uuid.uuid4())
    context  = memory.retrieve(session, query, n=3)

    return {
        **state,
        "intent":         intent,
        "sub_tasks":      sub_tasks,
        "sequences":      sequences,
        "memory_context": [{"text": c} for c in context],
        "session_id":     session,
        "agent_results":  {},
        "citations":      [],
        "error":          None,
    }


def route_node(state: BioAgentState) -> str:
    """LangGraph conditional edge — returns next node name."""
    intent = state["intent"]
    if intent == "compound":
        return "run_agents"
    return f"run_{intent}"


def _run_single_agent(intent: str, state: BioAgentState) -> Dict[str, Any]:
    """Dispatch to the appropriate specialist agent."""
    runner = AGENTS.get(intent)
    if runner is None:
        return {"result": f"No agent available for intent '{intent}'."}
    try:
        sequences = state.get("sequences", [])
        query     = state["query"]
        if intent == "literature":
            return runner(task=intent, query=query)
        elif intent in ("align", "phylo"):
            return runner(task=intent, sequences=sequences, query=query)
        else:
            seq = sequences[0] if sequences else ""
            return runner(task=intent, sequence=seq, query=query)
    except Exception as exc:
        return {"result": f"Agent error: {exc}", "error": str(exc)}


def agent_node_factory(intent: str):
    """Create a LangGraph node function for a specific intent."""
    def node(state: BioAgentState) -> BioAgentState:
        result = _run_single_agent(intent, state)
        return {
            **state,
            "agent_results": {**state.get("agent_results", {}), intent: result},
        }
    node.__name__ = f"run_{intent}"
    return node


def run_agents_node(state: BioAgentState) -> BioAgentState:
    """Execute all sub-tasks for a compound query."""
    results = dict(state.get("agent_results", {}))
    for intent in state.get("sub_tasks", []):
        results[intent] = _run_single_agent(intent, state)
    return {**state, "agent_results": results}


def fuse_node(state: BioAgentState) -> BioAgentState:
    """Assemble the final response from all agent results."""
    parts: List[str] = []
    citations: List[str] = []

    for intent, result in state["agent_results"].items():
        label = intent.upper()
        text  = result.get("result", "No result")
        parts.append(f"**{label}**\n{text}")

        # Collect PubMed citations from LitAgent
        if intent == "literature" and "citations" in result:
            citations += result["citations"]

    final_response = "\n\n".join(parts) if parts else "No results generated."

    # Store in session memory
    memory = get_memory()
    memory.store(
        state.get("session_id", "default"),
        state["query"],
        final_response[:500],
    )

    return {
        **state,
        "final_response": final_response,
        "citations":      citations,
    }


# ─── Graph construction ───────────────────────────────────────────────────────

def build_graph() -> StateGraph:
    g = StateGraph(BioAgentState)

    g.add_node("classify",     classify_node)
    g.add_node("run_seq",       agent_node_factory("seq"))
    g.add_node("run_align",     agent_node_factory("align"))
    g.add_node("run_annot",     agent_node_factory("annot"))
    g.add_node("run_phylo",     agent_node_factory("phylo"))
    g.add_node("run_literature",agent_node_factory("literature"))
    g.add_node("run_agents",    run_agents_node)
    g.add_node("fuse",          fuse_node)

    g.set_entry_point("classify")

    g.add_conditional_edges(
        "classify",
        route_node,
        {
            "run_seq":        "run_seq",
            "run_align":      "run_align",
            "run_annot":      "run_annot",
            "run_phylo":      "run_phylo",
            "run_literature": "run_literature",
            "run_agents":     "run_agents",
        },
    )

    for agent_node in ("run_seq", "run_align", "run_annot", "run_phylo", "run_literature", "run_agents"):
        g.add_edge(agent_node, "fuse")

    g.add_edge("fuse", END)
    return g.compile()


# ─── Singleton graph ──────────────────────────────────────────────────────────

_graph = None


def get_graph():
    global _graph
    if _graph is None:
        _graph = build_graph()
    return _graph


def run_bio_agent(query: str, session_id: str = "") -> BioAgentState:
    """Main entry point — run the full multi-agent pipeline."""
    initial: BioAgentState = {
        "query":          query,
        "intent":         "",
        "sub_tasks":      [],
        "sequences":      [],
        "agent_results":  {},
        "memory_context": [],
        "final_response": "",
        "citations":      [],
        "error":          None,
        "session_id":     session_id or str(uuid.uuid4()),
    }
    return get_graph().invoke(initial)
