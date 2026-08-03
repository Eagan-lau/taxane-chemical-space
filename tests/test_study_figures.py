import tempfile
import unittest
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.study_figures import (
    SourceRecorder,
    _format_element_delta,
    _portable_source_reference,
    _rule_label,
)


class StudyFigureTests(unittest.TestCase):
    def test_element_delta_uses_publication_notation(self):
        self.assertEqual(
            _format_element_delta('{"C": 2, "H": -2, "O": 1}'),
            "ΔC +2, ΔH -2, ΔO +1",
        )
        self.assertEqual(_format_element_delta("{}"), "No formula change")

    def test_rule_label_distinguishes_domain_and_global_rules(self):
        domain = pd.Series(
            {
                "grammar_rule_id": "SMRT_TAXANE_DOMAIN_CONSENSUS_000000004",
                "reaction_type": "O_acylation",
            }
        )
        global_rule = pd.Series(
            {
                "grammar_rule_id": "SMRT_T1_CORE_000192836",
                "reaction_type": "O_acylation",
            }
        )
        self.assertTrue(_rule_label(domain).startswith("TD-4:"))
        self.assertTrue(_rule_label(global_rule).startswith("T1-192836:"))

    def test_source_manifest_declares_scaffold_independence(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory)
            recorder = SourceRecorder(output)
            recorder.add(
                "FigX",
                "A",
                pd.DataFrame([{"generation": 1, "count": 2}]),
                source_path="analysis.tsv",
                description="Synthetic source table.",
            )
            manifest_path = recorder.write(output)
            manifest = pd.read_csv(manifest_path, sep="\t")
            self.assertNotIn("scaffold_fields_used", manifest.columns)
            source_table = pd.read_csv(
                output / manifest.loc[0, "source_data_file"], sep="\t"
            )
            self.assertFalse(
                any(
                    "scaffold" in column.lower()
                    for column in source_table.columns
                )
            )

    def test_figure_sources_use_portable_release_references(self):
        reference = _portable_source_reference(
            "/tmp/taxane_space_study_outputs/06_analysis_final/table.tsv;"
            "/tmp/taxane_space_study_outputs/01_provenance/source.tsv"
        )
        self.assertEqual(
            reference,
            "06_analysis_final/table.tsv;01_provenance/source.tsv",
        )


if __name__ == "__main__":
    unittest.main()
