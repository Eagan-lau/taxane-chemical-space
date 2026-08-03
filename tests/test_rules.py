import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from taxane_reaction_grammar_space.rules import (
    _semantic_group_key,
    prepare_taxane_domain_grammar,
)


class RuleTests(unittest.TestCase):
    def test_semantic_group_is_stable(self):
        row = {
            "reaction_delta_fingerprint": "O+1|H+0",
            "reaction_type": "hydroxylation_or_oxygenation",
            "normalized_direction": "substrate_to_product",
            "transferred_group_class": "oxygen",
            "donor_class": "",
            "acceptor_atom_class": "carbon",
        }
        self.assertEqual(_semantic_group_key(row), _semantic_group_key(dict(row)))

    def test_domain_harmonization_drops_scaffold_metadata(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "domain.tsv"
            pd.DataFrame(
                [
                    {
                        "smarts_rule_id": "domain_1",
                        "reaction_smarts": "[C;H1:1]>>[C:1][OH]",
                        "reaction_smarts_hash": "hash_1",
                        "reaction_delta_fingerprint": "O+1",
                        "reaction_type": "hydroxylation",
                        "taxane_domain_core_use": "true",
                        "predictive_rule_use": "true",
                        "template_qc_status": "ok",
                        "biochemical_step_granularity": (
                            "likely_single_enzyme_step"
                        ),
                        "composite_rule_flag": "false",
                        "taxane_domain_rule_id": "domain_1",
                        "taxane_consensus_type": "C_H_hydroxylation",
                        "taxane_scaffold_constraint": "must_not_survive",
                    }
                ]
            ).to_csv(source, sep="\t", index=False)
            outputs = prepare_taxane_domain_grammar(source, root / "output")
            grammar = pd.read_csv(outputs["grammar"], sep="\t")
            self.assertEqual(len(grammar), 1)
            self.assertFalse(
                any("scaffold" in column.lower() for column in grammar.columns)
            )


if __name__ == "__main__":
    unittest.main()
