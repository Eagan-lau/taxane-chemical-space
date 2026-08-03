import unittest

from taxane_reaction_grammar_space.chemistry import (
    atom_delta_matches,
    formula_delta,
    parse_reaction_delta_fingerprint,
    reaction_template_metrics,
)


class ChemistryTests(unittest.TestCase):
    def test_formula_delta(self):
        self.assertEqual(
            formula_delta("C20H30O3", "C22H32O4"),
            {"C": 2, "H": 2, "O": 1},
        )

    def test_parse_reaction_delta(self):
        self.assertEqual(
            parse_reaction_delta_fingerprint(
                "C+2|H+2|O+1|HA+3|R+0|M+42.011"
            ),
            {"C": 2, "H": 2, "O": 1},
        )

    def test_hydrogen_tolerance(self):
        self.assertTrue(
            atom_delta_matches(
                {"C": 2, "H": 1, "O": 1},
                {"C": 2, "H": 2, "O": 1},
                hydrogen_tolerance=1,
            )
        )

    def test_smarts_recomputed_single_center(self):
        metrics = reaction_template_metrics(
            "[C:1]-[O:2]>>[C:1]=[O:2]"
        )
        self.assertEqual(metrics.compile_status, "ok")
        self.assertEqual(metrics.inferred_independent_reaction_centers, 1)
        self.assertGreater(metrics.changed_mapped_atom_count, 0)
        self.assertTrue(metrics.reaction_edit_signature)


if __name__ == "__main__":
    unittest.main()
