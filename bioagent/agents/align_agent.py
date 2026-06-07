"""AlignAgent — pairwise and multiple sequence alignment."""

from __future__ import annotations

from typing import Any, Dict, List, Optional


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def pairwise_align(seq1: str, seq2: str) -> Dict[str, Any]:
    s1, s2 = _to_dna(seq1), _to_dna(seq2)
    try:
        from Bio import pairwise2
        alns = pairwise2.align.globalxx(s1, s2)
        if not alns:
            return {"result": "No alignment found.", "score": 0}
        best = alns[0]
        identity = sum(
            a == b for a, b in zip(best[0], best[1]) if a != "-" or b != "-"
        ) / max(len(s1), len(s2)) * 100
        return {
            "result": (
                f"Pairwise alignment score={best[2]:.0f}, "
                f"identity={identity:.1f}%\n"
                f"  {best[0][:60]}…\n"
                f"  {best[1][:60]}…"
            ),
            "score":      float(best[2]),
            "identity":   round(identity, 2),
            "aligned_a":  str(best[0]),
            "aligned_b":  str(best[1]),
        }
    except ImportError:
        # Naive needleman-wunsch with match=1, mismatch/gap=0
        score = sum(a == b for a, b in zip(s1, s2))
        return {
            "result": f"Simple alignment (BioPython unavailable): {score} matching bases",
            "score":  float(score),
        }


def simple_star_msa(sequences: List[str]) -> str:
    """Progressive star MSA against the first (reference) sequence."""
    if not sequences:
        return ""
    if len(sequences) == 1:
        return sequences[0]
    ref = _to_dna(sequences[0])
    aligned = [ref]
    for seq in sequences[1:]:
        result = pairwise_align(ref, seq)
        aligned.append(result.get("aligned_b", _to_dna(seq)))
    max_len = max(len(s) for s in aligned)
    return "\n".join(s + "-" * (max_len - len(s)) for s in aligned)


def run(
    task: str,
    sequences: Optional[List[str]] = None,
    sequence: str = "",
    **kwargs,
) -> Dict[str, Any]:
    seqs = sequences or kwargs.get("seqs", [])
    if sequence and not seqs:
        seqs = [sequence]

    if len(seqs) < 2:
        return {
            "result": "Provide at least 2 sequences for alignment.",
            "details": {},
        }

    task_l = task.lower()
    if len(seqs) > 2 or "msa" in task_l or "multiple" in task_l:
        msa = simple_star_msa(seqs)
        lines = msa.split("\n")
        preview = "\n".join(lines[:4]) + (f"\n… (+{len(lines)-4} more)" if len(lines) > 4 else "")
        return {
            "result": f"MSA ({len(seqs)} sequences):\n{preview}",
            "msa":    msa,
            "n_seqs": len(seqs),
        }

    return pairwise_align(seqs[0], seqs[1])
