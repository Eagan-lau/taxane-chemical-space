import tempfile
import unittest
import json
import sqlite3
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.analyze import analyze_chemical_space
from taxane_reaction_grammar_space.generate import (
    _compile_rules,
    _molecule_identity,
    _origin_retention,
    generate_chemical_space,
)


class GenerateTests(unittest.TestCase):
    def test_rule_compile_and_product_identity(self):
        from rdkit import Chem

        grammar = pd.DataFrame(
            [
                {
                    "grammar_rule_id": "g1",
                    "smarts_rule_id": "s1",
                    "reaction_smarts": "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]",
                    "semantic_group_id": "sem1",
                    "reaction_type": "oxidation",
                    "structural_element_delta": "",
                }
            ]
        )
        rules, failures = _compile_rules(grammar)
        self.assertFalse(failures)
        source = Chem.MolFromSmiles("CCO")
        outcomes = rules[0].reaction.RunReactants((source,))
        self.assertTrue(outcomes)
        product = outcomes[0][0]
        self.assertGreaterEqual(_origin_retention(product, source), 0.99)
        Chem.SanitizeMol(product)
        for atom in product.GetAtoms():
            atom.SetAtomMapNum(0)
        identity = _molecule_identity(product)
        self.assertEqual(identity["formula"], "C2H4O")
        self.assertEqual(len(identity["connectivity_key"]), 14)

    def test_end_to_end_one_generation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            grammar_path = root / "grammar.tsv"
            nodes_path = root / "nodes.tsv"
            output_dir = root / "output"
            pd.DataFrame(
                [
                    {
                        "grammar_rule_id": "g1",
                        "smarts_rule_id": "s1",
                        "reaction_smarts": (
                            "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]"
                        ),
                        "semantic_group_id": "sem1",
                        "reaction_type": "oxidation",
                        "structural_element_delta": "",
                        "evidence_layer_best": "T1_Bio_Core",
                        "final_rule_confidence": "0.9",
                    }
                ]
            ).to_csv(grammar_path, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "molecule_id": "known_ethanol",
                        "molecule_name": "ethanol",
                        "standardized_smiles": "CCO",
                    }
                ]
            ).to_csv(nodes_path, sep="\t", index=False)
            outputs = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=1,
            )
            nodes = pd.read_csv(outputs["nodes"], sep="\t")
            events = pd.read_csv(outputs["events"], sep="\t")
            self.assertEqual(set(nodes["generation_first"]), {0, 1})
            self.assertEqual(len(events), 1)
            self.assertEqual(events.iloc[0]["source_space_id"], "G0_00001")
            analysis = analyze_chemical_space(
                outputs["nodes"],
                outputs["events"],
                outputs["application_audit"],
                outputs["rejections"],
                root / "analysis",
            )
            generation_summary = pd.read_csv(
                analysis["generation_summary"], sep="\t"
            )
            self.assertEqual(
                generation_summary.loc[
                    generation_summary["generation"] == 1,
                    "unique_nodes_first_observed",
                ].iloc[0],
                1,
            )

            with sqlite3.connect(outputs["database"]) as connection:
                connection.execute("DELETE FROM generation_parent_progress")
                connection.commit()

            resumed = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=2,
                resume=True,
            )
            resumed_events = pd.read_csv(resumed["events"], sep="\t")
            self.assertEqual(len(resumed_events), 1)
            summary = json.loads(
                resumed["summary"].read_text(encoding="utf-8")
            )
            self.assertTrue(summary["resumed_from_existing_checkpoint"])
            self.assertEqual(summary["completed_generation_at_start"], 1)
            parent_progress = pd.read_csv(
                resumed["parent_progress"], sep="\t"
            )
            self.assertEqual(
                set(parent_progress["generation"].astype(int)), {1, 2}
            )

    def test_resume_preserves_an_incomplete_generation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            grammar_path = root / "grammar.tsv"
            nodes_path = root / "nodes.tsv"
            output_dir = root / "output"
            pd.DataFrame(
                [
                    {
                        "grammar_rule_id": "g1",
                        "smarts_rule_id": "s1",
                        "reaction_smarts": (
                            "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]"
                        ),
                        "semantic_group_id": "sem1",
                        "reaction_type": "oxidation",
                        "structural_element_delta": "",
                        "evidence_layer_best": "T1_Bio_Core",
                        "final_rule_confidence": "0.9",
                    }
                ]
            ).to_csv(grammar_path, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "molecule_id": "known_ethanol",
                        "molecule_name": "ethanol",
                        "standardized_smiles": "CCO",
                    }
                ]
            ).to_csv(nodes_path, sep="\t", index=False)
            first = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=2,
            )
            initial_nodes = pd.read_csv(first["nodes"], sep="\t")
            initial_events = pd.read_csv(first["events"], sep="\t")
            (output_dir / "G2_generation_summary.json").unlink()

            resumed = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=2,
                resume=True,
            )
            resumed_nodes = pd.read_csv(resumed["nodes"], sep="\t")
            resumed_events = pd.read_csv(resumed["events"], sep="\t")
            self.assertEqual(len(resumed_nodes), len(initial_nodes))
            self.assertEqual(len(resumed_events), len(initial_events))
            summary = json.loads(
                resumed["summary"].read_text(encoding="utf-8")
            )
            self.assertEqual(summary["completed_generation_at_start"], 1)
            self.assertTrue(
                summary["generations"][-1]["partial_generation_resumed"]
            )

    def test_resume_migrates_a_legacy_database_without_parent_ledger(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            grammar_path = root / "grammar.tsv"
            nodes_path = root / "nodes.tsv"
            output_dir = root / "output"
            pd.DataFrame(
                [
                    {
                        "grammar_rule_id": "g1",
                        "smarts_rule_id": "s1",
                        "reaction_smarts": (
                            "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]"
                        ),
                        "semantic_group_id": "sem1",
                        "reaction_type": "oxidation",
                        "structural_element_delta": "",
                        "evidence_layer_best": "T1_Bio_Core",
                        "final_rule_confidence": "0.9",
                    }
                ]
            ).to_csv(grammar_path, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "molecule_id": "known_ethanol",
                        "molecule_name": "ethanol",
                        "standardized_smiles": "CCO",
                    }
                ]
            ).to_csv(nodes_path, sep="\t", index=False)

            initial = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=1,
            )
            with sqlite3.connect(initial["database"]) as connection:
                connection.execute("DROP TABLE generation_parent_progress")
                connection.commit()

            resumed = generate_chemical_space(
                grammar_path,
                nodes_path,
                output_dir,
                max_generation=1,
                resume=True,
            )
            ledger = pd.read_csv(resumed["parent_progress"], sep="\t")
            self.assertEqual(len(ledger), 1)
            self.assertEqual(int(ledger.iloc[0]["generation"]), 1)
            self.assertEqual(ledger.iloc[0]["status"], "complete")


if __name__ == "__main__":
    unittest.main()
