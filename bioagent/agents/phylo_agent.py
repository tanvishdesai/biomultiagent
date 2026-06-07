"""PhyloAgent — Neighbor-joining phylogenetic tree construction and ASCII rendering."""

from __future__ import annotations

from typing import Any, Dict, List, Optional


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _hamming_distance(s1: str, s2: str) -> float:
    """Normalised Hamming distance between two sequences (same length)."""
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


def _upgma_tree(seqs: List[str], labels: List[str]) -> str:
    """Build a simple UPGMA-style Newick string from pairwise distances."""
    if len(seqs) < 2:
        return labels[0] if labels else ""
    dist = _distance_matrix([_to_dna(s) for s in seqs])
    n = len(seqs)
    nodes = list(labels)

    while len(nodes) > 1:
        # Find minimum distance pair
        min_d = float("inf")
        mi, mj = 0, 1
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if dist[i][j] < min_d:
                    min_d, mi, mj = dist[i][j], i, j
        merged = f"({nodes[mi]}:{min_d/2:.3f},{nodes[mj]}:{min_d/2:.3f})"
        nodes = [merged] + [nodes[k] for k in range(len(nodes)) if k not in (mi, mj)]
        # Update distances (average linkage)
        new_dist = [[0.0] * len(nodes) for _ in range(len(nodes))]
        old_indices = [k for k in range(len(dist)) if k not in (mi, mj)]
        for ni, oi in enumerate(old_indices):
            avg = (dist[mi][oi] + dist[mj][oi]) / 2
            new_dist[0][ni + 1] = avg
            new_dist[ni + 1][0] = avg
        for ni, oi in enumerate(old_indices):
            for nj, oj in enumerate(old_indices):
                new_dist[ni + 1][nj + 1] = dist[oi][oj]
        dist = new_dist

    return nodes[0] + ";"


def _ascii_tree(seqs: List[str], labels: List[str]) -> str:
    """Minimal ASCII dendrogram from pairwise distances."""
    if len(seqs) < 2:
        return labels[0] if labels else ""
    dist = _distance_matrix([_to_dna(s) for s in seqs])

    lines = []
    used = [False] * len(seqs)
    for i in range(len(seqs)):
        if used[i]:
            continue
        # Find closest neighbour
        min_d = float("inf")
        partner = i
        for j in range(len(seqs)):
            if j != i and not used[j] and dist[i][j] < min_d:
                min_d = dist[i][j]
                partner = j
        if partner != i and not used[partner]:
            lines.append(f"  ┌─ {labels[i]}")
            lines.append(f"  └─ {labels[partner]}  (d={min_d:.3f})")
            used[i] = used[partner] = True
        else:
            lines.append(f"  ── {labels[i]}")
    return "\n".join(lines)


def run(
    task: str,
    sequences: Optional[List[str]] = None,
    labels: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Any]:
    seqs = sequences or kwargs.get("seqs", [])
    if len(seqs) < 2:
        return {
            "result": "Provide at least 2 sequences to build a phylogenetic tree.",
            "details": {},
        }

    lbls = labels or [f"Seq{i+1}" for i in range(len(seqs))]
    newick = _upgma_tree(seqs, lbls)
    ascii_art = _ascii_tree(seqs, lbls)

    dist_matrix = _distance_matrix([_to_dna(s) for s in seqs])
    return {
        "result": f"Phylogenetic tree ({len(seqs)} sequences, UPGMA):\n{ascii_art}",
        "newick":       newick,
        "ascii_tree":   ascii_art,
        "dist_matrix":  [[round(d, 4) for d in row] for row in dist_matrix],
        "labels":       lbls,
    }
