import unittest

import pandas as pd
from rdkit import Chem

from taxane_reaction_grammar_space.analyze import (
    _chemical_space_projection,
    _compile_functional_queries,
    _convergence_and_paths,
    _functional_state,
    _functional_transition_label,
    _grammar_use_concentration,
    _interpretation_layer,
)


class AnalyzeTests(unittest.TestCase):
    def test_generation_boundaries_and_molecular_fingerprint_axis_names(self):
        self.assertEqual(_interpretation_layer(0), "known_taxane_seed_space")
        self.assertEqual(
            _interpretation_layer(1), "primary_near_seed_chemical_space"
        )
        self.assertEqual(
            _interpretation_layer(2), "primary_near_seed_chemical_space"
        )
        self.assertEqual(_interpretation_layer(3), "exploratory_frontier")

        nodes = pd.DataFrame(
            {
                "space_id": ["G0_A", "G1_A", "G2_A"],
                "generation_first": [0, 1, 2],
                "smiles": ["CCO", "CCOC", "CCOC(C)=O"],
                "full_inchikey": ["A", "B", "C"],
                "connectivity_key": ["A14", "B14", "C14"],
                "formula": ["C2H6O", "C3H8O", "C4H8O2"],
            }
        )
        projection = _chemical_space_projection(
            nodes, max_nodes_per_generation=10, random_seed=1729
        )
        self.assertIn("molecular_fp_axis_1", projection.columns)
        self.assertIn("molecular_fp_axis_2", projection.columns)
        self.assertNotIn("reaction_fp_axis_1", projection.columns)

    def test_functional_transition_detects_ester_gain(self):
        queries = _compile_functional_queries()
        alcohol = _functional_state(Chem.MolFromSmiles("CCO"), queries)
        ester = _functional_state(Chem.MolFromSmiles("CCOC(C)=O"), queries)
        transition = _functional_transition_label(alcohol, ester)
        self.assertIn("ester:+1", transition)
        self.assertIn("free_hydroxyl:-1", transition)

    def test_grammar_concentration_reports_dominant_rule(self):
        events = pd.DataFrame(
            {
                "generation": [1, 1, 1, 1],
                "grammar_rule_id": ["R1", "R1", "R1", "R2"],
                "semantic_group_id": ["S1", "S1", "S1", "S2"],
            }
        )
        result = _grammar_use_concentration(events).set_index(
            "generation_scope"
        )
        self.assertAlmostEqual(
            result.loc["G1", "top1_rule_event_fraction"], 0.75
        )
        self.assertEqual(result.loc["G1", "rules_used"], 2)

    def test_latent_bridge_uses_known_seed_lineage_without_scaffold(self):
        nodes = pd.DataFrame(
            {
                "space_id": ["G0_A", "G0_B", "G1_X"],
                "generation_first": [0, 0, 1],
                "smiles": ["C", "CO", "CC"],
                "formula": ["CH4", "CH4O", "C2H6"],
                "full_inchikey": ["A", "B", "X"],
            }
        )
        events = pd.DataFrame(
            {
                "event_id": [1, 2],
                "generation": [1, 2],
                "source_space_id": ["G0_A", "G1_X"],
                "target_space_id": ["G1_X", "G0_B"],
                "grammar_rule_id": ["R1", "R2"],
                "semantic_group_id": ["S1", "S2"],
                "target_is_new": [1, 0],
                "target_generation_first": [1, 0],
            }
        )
        _convergence, bridges, bridge_pairs = _convergence_and_paths(
            nodes, events
        )
        bridge = bridges.set_index("space_id").loc["G1_X"]
        self.assertTrue(bool(bridge["latent_bridge_candidate"]))
        self.assertEqual(bridge["distinct_G0_pair_bridge_count"], 1)
        self.assertEqual(len(bridge_pairs), 1)

    def test_convergence_and_paths_do_not_count_redundant_rules_as_routes(self):
        nodes = pd.DataFrame(
            {
                "space_id": ["G0_A", "G1_B", "G2_C"],
                "generation_first": [0, 1, 2],
                "smiles": ["C", "CO", "COC"],
                "formula": ["CH4", "CH4O", "C2H6O"],
                "full_inchikey": ["A", "B", "C"],
            }
        )
        events = pd.DataFrame(
            {
                "event_id": [1, 2, 3],
                "generation": [1, 1, 2],
                "source_space_id": ["G0_A", "G0_A", "G1_B"],
                "target_space_id": ["G1_B", "G1_B", "G2_C"],
                "grammar_rule_id": ["R1", "R2", "R3"],
                "semantic_group_id": ["S1", "S1", "S2"],
                "target_is_new": [1, 0, 1],
                "target_generation_first": [1, 1, 2],
            }
        )
        convergence, _bridges, _bridge_pairs = _convergence_and_paths(
            nodes, events
        )
        indexed = convergence.set_index("space_id")
        g1 = indexed.loc["G1_B"]
        g2 = indexed.loc["G2_C"]
        self.assertEqual(g1["unique_parent_count"], 1)
        self.assertEqual(g1["unique_rule_count"], 2)
        self.assertEqual(g1["unique_semantic_group_count"], 1)
        self.assertFalse(bool(g1["is_convergent"]))
        self.assertEqual(g1["structural_path_count"], 1)
        self.assertEqual(g1["semantic_edge_path_count"], 1)
        self.assertEqual(g1["raw_rule_event_path_count"], 2)
        self.assertEqual(g2["structural_path_count"], 1)
        self.assertEqual(g2["semantic_edge_path_count"], 1)
        self.assertEqual(g2["raw_rule_event_path_count"], 2)

    def test_path_counts_include_later_events_reaching_same_new_target(self):
        nodes = pd.DataFrame(
            {
                "space_id": ["G0_A", "G0_D", "G1_B", "G2_C"],
                "generation_first": [0, 0, 1, 2],
                "smiles": ["C", "N", "CO", "COC"],
                "formula": ["CH4", "NH3", "CH4O", "C2H6O"],
                "full_inchikey": ["A", "D", "B", "C"],
            }
        )
        events = pd.DataFrame(
            {
                "event_id": [1, 2, 3, 4],
                "generation": [1, 1, 1, 2],
                "source_space_id": ["G0_A", "G0_A", "G0_D", "G1_B"],
                "target_space_id": ["G1_B", "G1_B", "G1_B", "G2_C"],
                "grammar_rule_id": ["R1", "R2", "R4", "R3"],
                "semantic_group_id": ["S1", "S1", "S3", "S2"],
                "target_is_new": [1, 0, 0, 1],
                "target_generation_first": [1, 1, 1, 2],
            }
        )
        convergence, _bridges, _bridge_pairs = _convergence_and_paths(
            nodes, events
        )
        indexed = convergence.set_index("space_id")
        g1 = indexed.loc["G1_B"]
        g2 = indexed.loc["G2_C"]
        self.assertTrue(bool(g1["is_convergent"]))
        self.assertEqual(g1["structural_path_count"], 2)
        self.assertEqual(g1["semantic_edge_path_count"], 2)
        self.assertEqual(g1["raw_rule_event_path_count"], 3)
        self.assertEqual(g2["structural_path_count"], 2)
        self.assertEqual(g2["semantic_edge_path_count"], 2)
        self.assertEqual(g2["raw_rule_event_path_count"], 3)

    def test_later_generation_rediscovery_does_not_define_convergence(self):
        nodes = pd.DataFrame(
            {
                "space_id": ["G0_A", "G1_B", "G2_C"],
                "generation_first": [0, 1, 2],
                "smiles": ["C", "CO", "COC"],
                "formula": ["CH4", "CH4O", "C2H6O"],
                "full_inchikey": ["A", "B", "C"],
            }
        )
        events = pd.DataFrame(
            {
                "event_id": [1, 2, 3],
                "generation": [1, 2, 2],
                "source_space_id": ["G0_A", "G1_B", "G2_C"],
                "target_space_id": ["G1_B", "G2_C", "G1_B"],
                "grammar_rule_id": ["R1", "R2", "R3"],
                "semantic_group_id": ["S1", "S2", "S3"],
                "target_is_new": [1, 1, 0],
                "target_generation_first": [1, 2, 1],
            }
        )
        convergence, _bridges, _bridge_pairs = _convergence_and_paths(
            nodes, events
        )
        g1 = convergence.set_index("space_id").loc["G1_B"]
        self.assertEqual(g1["unique_parent_count"], 1)
        self.assertEqual(g1["all_unique_parent_count"], 2)
        self.assertEqual(g1["later_rediscovery_event_count"], 1)
        self.assertFalse(bool(g1["is_convergent"]))


if __name__ == "__main__":
    unittest.main()
