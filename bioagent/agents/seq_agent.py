"""SeqAgent — sequence translation, GC content, motif search, reverse complement."""

from __future__ import annotations

import re
from typing import Any, Dict, List


def _to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def gc_content(seq: str) -> float:
    s = _to_dna(seq)
    if not s:
        return 0.0
    gc = sum(1 for c in s if c in "GC")
    return round(100.0 * gc / len(s), 2)


def reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(c, c) for c in _to_dna(seq)[::-1])


def translate(seq: str, frame: int = 0) -> str:
    """Translate DNA/RNA sequence in the given reading frame (0, 1, or 2)."""
    try:
        from Bio.Seq import Seq
        dna = _to_dna(seq)[frame:]
        return str(Seq(dna).translate(to_stop=True))
    except ImportError:
        # Manual codon table (standard genetic code)
        codon_table = {
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
        dna = _to_dna(seq)[frame:]
        protein = []
        for i in range(0, len(dna) - 2, 3):
            aa = codon_table.get(dna[i:i+3], "X")
            if aa == "*":
                break
            protein.append(aa)
        return "".join(protein)


def find_motifs(seq: str, pattern: str = "ATG") -> List[int]:
    return [m.start() for m in re.finditer(pattern, seq.upper())]


def run(task: str, sequence: str = "", sequences: List[str] | None = None, **kwargs) -> Dict[str, Any]:
    seq = sequence or (sequences[0] if sequences else kwargs.get("seq", ""))
    if not seq:
        return {"result": "No sequence provided.", "details": {}}

    task_l = task.lower()
    details: Dict[str, Any] = {"length": len(seq), "gc_percent": gc_content(seq)}

    if "gc" in task_l or "content" in task_l:
        return {
            "result": f"GC content: {details['gc_percent']}%",
            "details": details,
        }

    if "reverse" in task_l or "complement" in task_l:
        rc = reverse_complement(seq)
        details["reverse_complement"] = rc[:80] + ("…" if len(rc) > 80 else "")
        return {"result": f"Reverse complement: {rc[:80]}…", "details": details}

    if "translate" in task_l or "protein" in task_l:
        proteins = {f: translate(seq, f) for f in (0, 1, 2)}
        best_frame = max(proteins, key=lambda f: len(proteins[f]))
        details["translations"] = {f"frame{f}": p[:100] for f, p in proteins.items()}
        return {
            "result": (
                f"Best translation (frame {best_frame}): {proteins[best_frame][:120]}"
                + ("…" if len(proteins[best_frame]) > 120 else "")
            ),
            "details": details,
        }

    if "motif" in task_l:
        pattern = kwargs.get("pattern", "ATG")
        positions = find_motifs(seq, pattern)
        details["motif_positions"] = positions
        return {
            "result": f"Motif '{pattern}' found at positions: {positions}",
            "details": details,
        }

    # Default: full analysis
    translations = {f"frame{f}": translate(seq, f)[:80] for f in (0, 1, 2)}
    details["translations"] = translations
    return {
        "result": (
            f"Sequence: {len(seq)} bp, GC={details['gc_percent']}%, "
            f"ATG starts={len(find_motifs(seq, 'ATG'))}"
        ),
        "details": details,
    }
