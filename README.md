# BioMultiAgent · Multi-Agent Bioinformatics AI System

> **Part of the Bioinformatics AI Portfolio** — Project 2 of 4.
> See [PLANNING.md](../PLANNING.md) for the full roadmap.

A LangGraph-orchestrated multi-agent system where a Supervisor routes compound bioinformatics queries to specialist agents, each owning one analysis domain.

Directly extends [Suru1496/Bio_NLP](https://github.com/Suru1496/Bio_NLP) with multi-agent orchestration, compound task decomposition, and session memory — addressing the "future scope" items in Surbhi Pawar's BioNLP Platform report.

---

## Architecture

```
User query (natural language)
      │
      ▼
Supervisor (LangGraph) — classify intent → route
      │
      ├──▶ SeqAgent     translation, GC%, motif search, reverse complement
      ├──▶ AlignAgent   pairwise + star MSA (pure Python fallback, no MUSCLE needed)
      ├──▶ AnnotAgent   ORF prediction in all 6 reading frames, codon analysis
      ├──▶ PhyloAgent   UPGMA tree, Newick output, ASCII dendrogram
      └──▶ LitAgent     PubMed Entrez search + abstract retrieval
                │
                ▼
       Session Memory (ChromaDB)
       FuseNode → natural language response + citations
```

**Compound queries** fire multiple agents in parallel and fuse their results into one coherent response.

---

## Kaggle Setup (Recommended — No Data Download Required)

Seq/align/annot/phylo agents run locally. LitAgent requires network access to NCBI PubMed (Bio.Entrez). Set `ENTREZ_EMAIL` to your real email before running.

```python
!pip install langgraph langchain-community biopython flask chromadb -q

# Run integration demo (no data needed)
!python integration_demo.py

# Or start the web app
!python app.py
```

---

## Local Setup

```bash
# 1. Install dependencies
pip install -r requirements.txt

# Optional: install Ollama for LLM-based intent classification
# https://ollama.ai
ollama pull mistral

# 2. Run demo
python integration_demo.py

# 3. Start web app
python app.py
# Open: http://localhost:5001
```

---

## Run Order

```
1. pip install -r requirements.txt      ← install deps
2. python integration_demo.py           ← verify all agents work
3. python app.py                        ← start web interface
```

---

## Example Compound Query

```
"Translate ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG,
 predict ORFs, and find 5 recent PubMed papers on KRAS mutations"
```

**Output:**
- `SeqAgent` → translation in 3 frames, GC content
- `AnnotAgent` → ORF list with positions and protein sequences
- `LitAgent` → 5 PubMed abstracts with citations

One query. Three agents. One coherent response.

---

## Project Structure

```
biomultiagent/
  app.py                    Flask web API
  integration_demo.py       End-to-end compound query demo
  requirements.txt
  bioagent/
    supervisor.py           LangGraph state machine + routing
    state.py                BioAgentState TypedDict schema
    agents/
      seq_agent.py          SeqAgent
      align_agent.py        AlignAgent
      annot_agent.py        AnnotAgent
      phylo_agent.py        PhyloAgent
      lit_agent.py          LitAgent (PubMed Entrez)
    memory/
      chroma_store.py       ChromaDB session memory (with fallback)
  templates/
    index.html              Web UI
```

---

## Collaboration Context

Surbhi Pawar's BioNLP Platform report explicitly lists "multi-agent AI systems" as future scope.
This project is the direct implementation of that roadmap using LangGraph.

**PR strategy**: Fork `Suru1496/Bio_NLP` → add `bioagent/` module → open PR with description referencing her Future Scope section.
