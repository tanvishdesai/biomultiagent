"""AlignAgent — pairwise and multiple sequence alignment (BioPython / MUSCLE)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
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
    except ImportError as exc:
        raise RuntimeError("BioPython is required for sequence alignment.") from exc


def _write_fasta(path: Path, sequences: List[str], labels: List[str]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for label, seq in zip(labels, sequences):
            fh.write(f">{label}\n{_to_dna(seq)}\n")


def muscle_msa(sequences: List[str], labels: Optional[List[str]] = None) -> str:
    """Run MUSCLE v3+ if installed; raises if unavailable."""
    muscle = shutil.which("muscle")
    if not muscle:
        raise RuntimeError(
            "MUSCLE not found on PATH. Install MUSCLE or use ≤2 sequences for pairwise alignment."
        )
    lbls = labels or [f"seq{i+1}" for i in range(len(sequences))]
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "in.fasta"
        out = Path(tmp) / "out.fasta"
        _write_fasta(inp, sequences, lbls)
        subprocess.run(
            [muscle, "-align", str(inp), "-output", str(out)],
            check=True,
            capture_output=True,
            text=True,
        )
        return out.read_text(encoding="utf-8")


def run(
    task: str,
    sequences: Optional[List[str]] = None,
    sequence: str = "",
    query: str = "",
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

    task_l = f"{task} {query}".lower()
    if len(seqs) > 2 or "msa" in task_l or "multiple" in task_l:
        try:
            msa = muscle_msa(seqs)
            lines = [ln for ln in msa.splitlines() if not ln.startswith(">")]
            preview = "\n".join(lines[:4])
            if len(lines) > 4:
                preview += f"\n… (+{len(lines)-4} more)"
            return {
                "result": f"MUSCLE MSA ({len(seqs)} sequences):\n{preview}",
                "msa":    msa,
                "n_seqs": len(seqs),
            }
        except RuntimeError as exc:
            return {"result": str(exc), "error": str(exc)}

    return pairwise_align(seqs[0], seqs[1])
