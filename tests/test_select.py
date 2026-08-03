import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from taxane_reaction_grammar_space.select import (
    _select_tier,
    assemble_open_grammar,
    augment_open_grammar_with_domain,
)


class SelectTests(unittest.TestCase):
    def test_primary_rejects_unassigned_large_context_rule(self):
        frame = pd.DataFrame(
            [
                {
                    "reaction_smarts": "[O:1]C(=O)C>>[O:1]C(=O)CCCCCCCCCCCCCCCCCCCC",
                    "reaction_type": "unassigned_ec_supported",
                    "biochemical_step_granularity": "uncertain",
                    "single_center_evidence_mode": "smarts_recomputed",
                    "consensus_generation_mode": "",
                    "consensus_qc_status": "",
                    "consensus_source_database_count": "0",
                    "template_sources": "RetroRules",
                    "structural_element_delta": "C+18",
                    "inferred_independent_reaction_centers": "1",
                    "mapped_atom_retention": "1",
                    "reactant_atoms": "4",
                    "grammar_selection_score": "0.9",
                    "g0_match_count": "600",
                    "semantic_group_id": "bad",
                    "reaction_smarts_hash": "x",
                }
            ]
        )
        selected, audit, reasons = _select_tier(
            frame,
            tier="primary",
            representatives_per_group=3,
            max_heavy_atom_gain=24,
            max_product_pattern_growth=32,
        )
        self.assertTrue(selected.empty)
        self.assertIn("reaction_semantics_unassigned", reasons)

    def test_primary_rejects_unassigned_structural_delta(self):
        frame = pd.DataFrame(
            [
                {
                    "reaction_smarts": "[C;H1:1]>>[C:1]O",
                    "reaction_type": "unassigned_structural_delta",
                    "biochemical_step_granularity": "likely_single_enzyme_step",
                    "single_center_evidence_mode": "source_annotation",
                    "consensus_generation_mode": "",
                    "consensus_qc_status": "",
                    "consensus_source_database_count": "1",
                    "template_sources": "BioNaviNP_BioChem",
                    "structural_element_delta": "O+1",
                    "inferred_independent_reaction_centers": "1",
                    "mapped_atom_retention": "1",
                    "reactant_atoms": "1",
                    "grammar_selection_score": "0.9",
                    "g0_match_count": "1",
                    "semantic_group_id": "unassigned",
                    "reaction_smarts_hash": "unassigned",
                }
            ]
        )
        selected, _audit, reasons = _select_tier(
            frame,
            tier="primary",
            representatives_per_group=3,
            max_heavy_atom_gain=24,
            max_product_pattern_growth=32,
        )
        self.assertTrue(selected.empty)
        self.assertIn("reaction_semantics_unassigned", reasons)

    def test_primary_accepts_validated_single_step_consensus(self):
        frame = pd.DataFrame(
            [
                {
                    "reaction_smarts": "[O;H1:1]>>[O:1]C(C)=O",
                    "reaction_type": "O_acetylation",
                    "biochemical_step_granularity": "likely_single_enzyme_step",
                    "single_center_evidence_mode": "source_annotation",
                    "consensus_generation_mode": "transfer_consensus",
                    "consensus_qc_status": "passed",
                    "consensus_source_database_count": "2",
                    "template_sources": "Rhea;BioNaviNP_BioChem",
                    "structural_element_delta": "C+2|O+1",
                    "inferred_independent_reaction_centers": "1",
                    "mapped_atom_retention": "1",
                    "reactant_atoms": "1",
                    "grammar_selection_score": "0.9",
                    "g0_match_count": "400",
                    "semantic_group_id": "good",
                    "reaction_smarts_hash": "y",
                }
            ]
        )
        selected, _audit, reasons = _select_tier(
            frame,
            tier="primary",
            representatives_per_group=3,
            max_heavy_atom_gain=24,
            max_product_pattern_growth=32,
        )
        self.assertEqual(len(selected), 1)
        self.assertFalse(reasons)

    def test_open_grammar_tracks_domain_and_global_scope(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            global_frame = pd.DataFrame(
                [
                    {
                        "reaction_smarts_hash": "shared",
                        "semantic_group_id": "group_a",
                        "taxane_selection_score": 1.0,
                        "reactant_atoms": 6,
                    },
                    {
                        "reaction_smarts_hash": "global_only",
                        "semantic_group_id": "group_b",
                        "taxane_selection_score": 1.0,
                        "reactant_atoms": 6,
                    },
                ]
            )
            domain_frame = pd.DataFrame(
                [
                    {
                        "reaction_smarts_hash": "shared",
                        "semantic_group_id": "group_a",
                        "taxane_selection_score": 0.5,
                        "reactant_atoms": 5,
                    },
                    {
                        "reaction_smarts_hash": "domain_only",
                        "semantic_group_id": "group_a",
                        "taxane_selection_score": 0.4,
                        "reactant_atoms": 5,
                    },
                ]
            )
            global_path = root / "global.tsv"
            domain_path = root / "domain.tsv"
            global_frame.to_csv(global_path, sep="\t", index=False)
            domain_frame.to_csv(domain_path, sep="\t", index=False)
            outputs = assemble_open_grammar(
                global_path,
                domain_path,
                root / "output",
                representatives_per_group=4,
            )
            result = pd.read_csv(outputs["open_grammar"], sep="\t")
            scopes = dict(
                zip(result["reaction_smarts_hash"], result["open_grammar_scope"])
            )
            self.assertEqual(scopes["shared"], "both")
            self.assertEqual(scopes["domain_only"], "G0_domain")
            self.assertEqual(scopes["global_only"], "global")

    def test_domain_augmentation_tracks_provenance(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            external_path = root / "external.tsv"
            domain_path = root / "domain.tsv"
            pd.DataFrame(
                [
                    {
                        "reaction_smarts_hash": "shared",
                        "reaction_smarts": "[C:1]>>[C:1]O",
                        "semantic_group_id": "a",
                    },
                    {
                        "reaction_smarts_hash": "external",
                        "reaction_smarts": "[O:1]>>[O:1]C",
                        "semantic_group_id": "b",
                    },
                ]
            ).to_csv(external_path, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "reaction_smarts_hash": "shared",
                        "reaction_smarts": "[C:1]>>[C:1]O",
                        "semantic_group_id": "a",
                    },
                    {
                        "reaction_smarts_hash": "domain",
                        "reaction_smarts": "[N:1]>>[N:1]C",
                        "semantic_group_id": "c",
                    },
                ]
            ).to_csv(domain_path, sep="\t", index=False)
            outputs = augment_open_grammar_with_domain(
                external_path,
                domain_path,
                root / "output",
            )
            result = pd.read_csv(outputs["grammar"], sep="\t")
            scopes = dict(
                zip(
                    result["reaction_smarts_hash"],
                    result["grammar_provenance_scope"],
                )
            )
            self.assertEqual(scopes["shared"], "external_and_taxane_domain")
            self.assertEqual(scopes["external"], "external")
            self.assertEqual(scopes["domain"], "taxane_domain")


if __name__ == "__main__":
    unittest.main()
