"""AnnotAgent — ORF prediction, codon analysis, gene annotation."""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def find_orfs(seq: str, min_len: int = 30) -> List[Dict]:
    """
    Find all ORFs in all 6 reading frames (3 forward + 3 reverse complement).
    Returns list of {frame, start, stop, length, protein}.
    """
    def _rc(s: str) -> str:
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(comp.get(c, c) for c in s[::-1])

    def _orfs_in_frame(dna: str, frame: int, reverse: bool) -> List[Dict]:
        orfs = []
        seq = dna[frame:]
        i = 0
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            if codon == "ATG":
                start = i + frame
                protein = ["M"]
                j = i + 3
                while j < len(seq) - 2:
                    c = seq[j:j+3]
                    if c in STOP_CODONS:
                        stop = j + frame + 3
                        aa = "".join(protein)
                        if len(aa) * 3 >= min_len:
                            orfs.append({
                                "frame":   (-frame - 1) if reverse else frame,
                                "start":   start,
                                "stop":    stop,
                                "length":  stop - start,
                                "protein": aa[:50] + ("…" if len(aa) > 50 else ""),
                            })
                        break
                    protein.append(CODON_TABLE.get(c, "X"))
                    j += 3
                i = j + 3 if j < len(seq) - 2 else i + 3
            else:
                i += 3
        return orfs

    dna = _to_dna(seq)
    rc  = _rc(dna)
    orfs = []
    for f in range(3):
        orfs += _orfs_in_frame(dna, f, False)
        orfs += _orfs_in_frame(rc, f, True)
    return sorted(orfs, key=lambda o: -o["length"])


def codon_usage(seq: str) -> Dict[str, float]:
    dna = _to_dna(seq)
    counts: Dict[str, int] = {}
    total = 0
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i+3]
        if len(c) == 3 and c not in STOP_CODONS:
            counts[c] = counts.get(c, 0) + 1
            total += 1
    if total == 0:
        return {}
    return {c: round(n / total * 100, 2) for c, n in sorted(counts.items())}


def run(
    task: str,
    sequence: str = "",
    sequences: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Any]:
    seq = sequence or (sequences[0] if sequences else kwargs.get("seq", ""))
    if not seq:
        return {"result": "No sequence provided.", "details": {}}

    q = re.sub(r"\b[ACGTUN]{10,}\b", " ", kwargs.get("query", ""), flags=re.I)
    task_l = f"{task} {q}".lower()
    orfs = find_orfs(seq)

    if "codon" in task_l:
        usage = codon_usage(seq)
        top_codons = dict(list(usage.items())[:10])
        return {
            "result": f"Top 10 codons by frequency: {top_codons}",
            "details": {"codon_usage": usage},
        }

    if not orfs:
        return {"result": "No ORFs found (min 30 bp).", "details": {"n_orfs": 0}}

    summary = f"Found {len(orfs)} ORFs. Longest: {orfs[0]['length']} bp, protein: {orfs[0]['protein']}"
    orf_list = orfs[:5]
    return {
        "result": summary,
        "details": {
            "n_orfs":    len(orfs),
            "top_orfs":  orf_list,
        },
    }
