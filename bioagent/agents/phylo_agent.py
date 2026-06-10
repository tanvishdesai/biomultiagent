"""PhyloAgent — phylogenetic tree construction via Bio.Phylo (Neighbor-Joining)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _hamming_distance(s1: str, s2: str) -> float:
    n = min(len(s1), len(s2))
    if n == 0:
        return 0.0
    mismatches = sum(a != b for a, b in zip(s1[:n], s2[:n]))
    return mismatches / n


def _distance_matrix(seqs: List[str]) -> List[List[float]]:
    n = len(seqs)
    d = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist = _hamming_distance(seqs[i], seqs[j])
            d[i][j] = dist
            d[j][i] = dist
    return d


def _nj_tree(seqs: List[str], labels: List[str]) -> tuple[str, str]:
    """Build NJ tree with Bio.Phylo; return (newick, ascii)."""
    import sys
    from io import StringIO

    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    dm = DistanceMatrix(labels, _distance_matrix([_to_dna(s) for s in seqs]))
    tree = DistanceTreeConstructor().nj(dm)

    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    newick = buf.getvalue().strip()

    ascii_buf = StringIO()
    old_stdout = sys.stdout
    sys.stdout = ascii_buf
    try:
        Phylo.draw_ascii(tree)
    finally:
        sys.stdout = old_stdout
    ascii_art = ascii_buf.getvalue().strip() or newick
    return newick, ascii_art


def run(
    task: str,
    sequences: Optional[List[str]] = None,
    labels: Optional[List[str]] = None,
    query: str = "",
    **kwargs,
) -> Dict[str, Any]:
    seqs = sequences or kwargs.get("seqs", [])
    if len(seqs) < 2:
        return {
            "result": "Provide at least 2 sequences to build a phylogenetic tree.",
            "details": {},
        }

    lbls = labels or [f"Seq{i+1}" for i in range(len(seqs))]
    try:
        newick, ascii_art = _nj_tree(seqs, lbls)
    except Exception as exc:
        return {"result": f"Phylogeny failed: {exc}", "error": str(exc)}

    dist_matrix = _distance_matrix([_to_dna(s) for s in seqs])
    return {
        "result": f"Neighbor-joining tree ({len(seqs)} sequences):\n{ascii_art}",
        "newick":       newick,
        "ascii_tree":   ascii_art,
        "dist_matrix":  [[round(d, 4) for d in row] for row in dist_matrix],
        "labels":       lbls,
    }
