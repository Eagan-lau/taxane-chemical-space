import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.provenance import (
    summarize_rule_library_provenance,
)


class ProvenanceTests(unittest.TestCase):
    def test_source_support_is_overlapping_and_scaffold_free(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "01_sources").mkdir()
            (root / "04_release").mkdir()
            source = pd.DataFrame(
                {
                    "source_database": ["Rhea", "BioNavi", "Rhea"],
                }
            )
            source.to_csv(
                root / "01_sources/source_reactions.normalized.tsv",
                sep="\t",
                index=False,
            )
            source.to_csv(
                root / "01_sources/source_reactions.main_pair.tsv",
                sep="\t",
                index=False,
            )
            release = pd.DataFrame(
                {"template_sources": ["Rhea;BioNavi", "Rhea"]}
            )
            for tier in ("T1", "T2", "T3"):
                release.to_csv(
                    root
                    / f"04_release/reaction_smarts_library.{tier}_only.tsv",
                    sep="\t",
                    index=False,
                )
            (root / "build_summary.json").write_text(
                json.dumps(
                    {
                        "source_records": 3,
                        "main_pair_records": 3,
                        "raw_templates_or_anchors": 3,
                        "deduplicated_templates_or_anchors": 2,
                        "predictive_generalized_rules": 6,
                        "reaction_smarts_library_T1_only": 2,
                        "reaction_smarts_library_T2_only": 2,
                        "reaction_smarts_library_T3_only": 2,
                    }
                ),
                encoding="utf-8",
            )
            grammar = pd.DataFrame(
                {
                    "final_grammar_rule_id": ["F1"],
                    "reaction_smarts": ["[C:1]>>[C:1]O"],
                    "reaction_type": ["hydroxylation"],
                    "template_sources": ["Rhea;BioNavi"],
                }
            )
            grammar_path = root / "final.tsv"
            grammar.to_csv(grammar_path, sep="\t", index=False)
            output = root / "out"
            paths = summarize_rule_library_provenance(
                root, grammar_path, output
            )
            contributions = pd.read_csv(
                paths["source_contributions"], sep="\t"
            ).set_index("source_database")
            self.assertEqual(contributions.loc["Rhea", "T1_rule_rows_supported"], 2)
            self.assertEqual(
                contributions.loc["BioNavi", "T1_rule_rows_supported"], 1
            )
            summary = json.loads(paths["summary"].read_text())
            self.assertTrue(summary["source_support_counts_are_overlapping"])
            self.assertNotIn("scaffold_fields_used", summary)


if __name__ == "__main__":
    unittest.main()
