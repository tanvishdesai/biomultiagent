"""Agent registry — maps intent names to agent run() functions."""

from . import align_agent, annot_agent, lit_agent, phylo_agent, seq_agent

AGENTS = {
    "seq":        seq_agent.run,
    "align":      align_agent.run,
    "annot":      annot_agent.run,
    "phylo":      phylo_agent.run,
    "literature": lit_agent.run,
}
