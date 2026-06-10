"""Unit tests for BioMultiAgent specialist agents."""

from __future__ import annotations

import unittest
from unittest.mock import patch

from bioagent.agents import annot_agent, seq_agent


class TestSeqAgent(unittest.TestCase):
    def test_translate_from_query(self):
        seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        out = seq_agent.run(task="seq", sequence=seq, query="Translate to protein")
        self.assertIn("translation", out["result"].lower())
        self.assertNotIn("GC=", out["result"])

    def test_gc_content(self):
        out = seq_agent.run(task="seq", sequence="GCGC", query="GC content")
        self.assertIn("100.0%", out["result"])


class TestAnnotAgent(unittest.TestCase):
    def test_orf_on_long_sequence(self):
        seq = "ATG" + "GCC" * 20 + "TAA"
        out = annot_agent.run(task="annot", sequence=seq, query="predict ORFs")
        self.assertIn("ORF", out["result"])


class TestLitAgent(unittest.TestCase):
    @patch("bioagent.agents.lit_agent.search_pubmed")
    @patch("bioagent.agents.lit_agent.fetch_abstracts")
    def test_pubmed_success(self, mock_fetch, mock_search):
        from bioagent.agents import lit_agent
        mock_search.return_value = (["12345"], None)
        mock_fetch.return_value = ([{
            "pmid": "12345", "title": "RNA structure paper",
            "abstract": "Secondary structure prediction methods.",
            "year": "2024", "journal": "Bioinformatics", "url": "http://x",
        }], None)
        out = lit_agent.run(query="Find papers on RNA structure")
        self.assertIn("12345", out["result"])
        self.assertEqual(len(out["citations"]), 1)

    @patch("bioagent.agents.lit_agent.search_pubmed")
    def test_pubmed_network_error(self, mock_search):
        from bioagent.agents import lit_agent
        mock_search.return_value = ([], "PubMed search failed: network down")
        out = lit_agent.run(query="Find papers on RNA")
        self.assertIn("failed", out["result"].lower())
        self.assertIn("error", out)


if __name__ == "__main__":
    unittest.main()
