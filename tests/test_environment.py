import tempfile
import unittest
from pathlib import Path

import pandas as pd

from taxane_reaction_grammar_space.environment import record_environment


class EnvironmentTests(unittest.TestCase):
    def test_environment_snapshot_without_inputs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = record_environment(Path(temporary_directory))
            environment = pd.read_csv(
                paths["software_environment"], sep="\t"
            )
            self.assertIn("python", set(environment["component"]))
            hashes = pd.read_csv(paths["input_file_hashes"], sep="\t")
            self.assertEqual(len(hashes), 0)


if __name__ == "__main__":
    unittest.main()
