import tempfile
import unittest
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.manuscript import render_manuscript


class ManuscriptTests(unittest.TestCase):
    def test_final_context_is_rendered_without_scaffold_fields(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            analysis = root / "analysis"
            output = root / "output"
            analysis.mkdir()

            pd.DataFrame(
                {
                    "generation": [0, 1, 2, 3],
                    "unique_nodes_first_observed": [2, 3, 4, 5],
                    "cumulative_unique_nodes": [2, 5, 9, 14],
                    "derivation_events": [0, 4, 5, 6],
                    "unique_node_yield_per_raw_product": [
                        0.0,
                        0.5,
                        0.4,
                        0.3,
                    ],
                }
            ).to_csv(
                analysis / "generation_expansion_summary.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "generation_scope": [
                        "all_generations",
                        "G1",
                        "G2",
                        "G3",
                    ],
                    "top1_rule_event_fraction": [0.5, 0.6, 0.5, 0.4],
                    "top5_rule_event_fraction": [0.9, 0.9, 0.8, 0.7],
                    "effective_rule_number_exp_shannon": [
                        3.0,
                        2.0,
                        3.0,
                        4.0,
                    ],
                }
            ).to_csv(
                analysis / "reaction_grammar_use_concentration.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "grammar_rule_id": ["R1", "R2"],
                    "semantic_group_id": ["S1", "S1"],
                }
            ).to_csv(
                analysis / "reaction_grammar_usage.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "generation": [1, 2, 3],
                    "nearest_G0_tanimoto": [0.8, 0.7, 0.6],
                }
            ).to_csv(
                analysis / "nearest_G0_similarity.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "generation": [0, 3],
                    "descriptor": ["exact_mass", "exact_mass"],
                    "mean": [500.0, 510.0],
                    "standard_deviation": [10.0, 11.0],
                }
            ).to_csv(
                analysis / "physicochemical_descriptor_summary.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "generation": [1, 1, 2, 2, 3, 3],
                    "is_convergent": [
                        True,
                        False,
                        True,
                        True,
                        False,
                        False,
                    ],
                }
            ).to_csv(
                analysis / "convergence_and_route_multiplicity.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "latent_bridge_candidate": [True, False],
                    "distinct_G0_pair_bridge_count": [2, 0],
                }
            ).to_csv(
                analysis / "latent_bridge_candidates.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {"known_G0_source_space_id": ["G0_A"]}
            ).to_csv(
                analysis / "known_G0_pair_bridge_summary.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame(
                {
                    "observed_changed_source_atoms": [1, 4],
                    "derivation_event_count": [9, 1],
                }
            ).to_csv(
                analysis / "reaction_edit_landscape.tsv",
                sep="\t",
                index=False,
            )
            template = root / "template.md"
            template.write_text(
                "G3={{FINAL_G3_UNIQUE_STRUCTURES}}; "
                "bridges={{FINAL_LATENT_BRIDGE_COUNT}}",
                encoding="utf-8",
            )

            paths = render_manuscript(template, analysis, output)
            manuscript = paths["manuscript"].read_text(encoding="utf-8")
            self.assertEqual(manuscript, "G3=5; bridges=1")
            summary = paths["summary"].read_text(encoding="utf-8")
            self.assertNotIn("scaffold_fields_used", summary)


if __name__ == "__main__":
    unittest.main()
