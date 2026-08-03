import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.sensitivity import (
    compare_g1_sensitivity_spaces,
)


class SensitivityTests(unittest.TestCase):
    def test_full_inchikey_overlap_is_reported(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            specifications = {
                "primary": ["A", "B"],
                "t2": ["B", "C"],
                "t3": ["B", "D"],
            }
            paths = {}
            for label, keys in specifications.items():
                frame = pd.DataFrame(
                    {
                        "generation_first": [0] + [1] * len(keys),
                        "full_inchikey": ["G0"] + keys,
                        "connectivity_key": ["G0"] + keys,
                        "smiles": ["C"] * (len(keys) + 1),
                        "formula": ["CH4"] * (len(keys) + 1),
                    }
                )
                nodes = root / f"{label}.tsv"
                frame.to_csv(nodes, sep="\t", index=False)
                summary = root / f"{label}.json"
                summary.write_text(
                    json.dumps(
                        {
                            "compiled_grammar_rules": 1,
                            "grammar_compile_failures": 0,
                            "G0": {"unique_full_stereo_structures": 1},
                            "generations": [
                                {
                                    "generation": 1,
                                    "activated_rules": 1,
                                    "accepted_derivation_events": len(keys),
                                    "raw_product_tuples": len(keys),
                                    "rejection_reason_counts": {},
                                }
                            ],
                        }
                    ),
                    encoding="utf-8",
                )
                paths[label] = (nodes, summary)
            outputs = compare_g1_sensitivity_spaces(
                paths["primary"][0],
                paths["primary"][1],
                paths["t2"][0],
                paths["t2"][1],
                paths["t3"][0],
                paths["t3"][1],
                root / "out",
            )
            overlap = pd.read_csv(outputs["pairwise_overlap"], sep="\t")
            row = overlap.iloc[0]
            self.assertEqual(row["intersection_G1_structures"], 1)
            summary = json.loads(outputs["summary"].read_text())
            self.assertEqual(summary["G1_shared_by_all_three_layers"], 1)
            self.assertNotIn("scaffold_fields_used", summary)


if __name__ == "__main__":
    unittest.main()
