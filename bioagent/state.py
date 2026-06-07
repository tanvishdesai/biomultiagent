"""BioAgentState — LangGraph typed state schema for the multi-agent system."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from typing_extensions import TypedDict


class BioAgentState(TypedDict):
    """
    Central shared state passed between all nodes in the LangGraph graph.

    Fields
    ------
    query         : Original user natural language query.
    intent        : Classified intent — seq | align | annot | phylo | literature | compound.
    sub_tasks     : List of intents to execute for compound queries.
    sequences     : Extracted DNA/RNA sequences from the query.
    agent_results : Dict mapping intent name → result dict from each specialist agent.
    memory_context: Retrieved session history (previous Q&A pairs).
    final_response: Assembled natural language answer for the user.
    citations     : List of source strings (PubMed PMIDs, database names, etc.).
    error         : Optional error message if a node failed.
    session_id    : Session identifier for memory persistence across queries.
    """

    query         : str
    intent        : str
    sub_tasks     : List[str]
    sequences     : List[str]
    agent_results : Dict[str, Any]
    memory_context: List[Dict[str, str]]
    final_response: str
    citations     : List[str]
    error         : Optional[str]
    session_id    : str
