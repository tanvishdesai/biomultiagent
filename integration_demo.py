"""
BioMultiAgent Integration Demo
================================
End-to-end compound query demonstration.
Run: python integration_demo.py

Sequence/annotation agents run locally; LitAgent requires network access to PubMed (NCBI Entrez).
"""

from __future__ import annotations

import json

from bioagent.supervisor import run_bio_agent

DEMO_QUERIES = [
    # Single-agent
    "Translate ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG to protein",
    "Find 5 recent PubMed papers on RNA secondary structure prediction",
    # Compound
    (
        "Translate ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG, "
        "predict ORFs, and find papers on KRAS mutations in cancer"
    ),
]


def main() -> None:
    print("=" * 70)
    print("BioMultiAgent Integration Demo")
    print("=" * 70)

    for i, query in enumerate(DEMO_QUERIES, 1):
        print(f"\n{'-'*70}")
        print(f"Query {i}: {query[:80]}{'...' if len(query) > 80 else ''}")
        print("-" * 70)

        state = run_bio_agent(query)

        print(f"Intent   : {state['intent']}")
        print(f"Sub-tasks: {state['sub_tasks']}")
        print(f"\nResponse:\n{state['final_response'][:800]}")

        if state["citations"]:
            print(f"\nCitations ({len(state['citations'])}):")
            for c in state["citations"][:3]:
                print(f"  · {c}")

    print("\n" + "=" * 70)
    print("Demo complete. Run `python app.py` for the web interface.")


if __name__ == "__main__":
    main()
