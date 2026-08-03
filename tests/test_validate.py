import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from taxane_reaction_grammar_space.generate import _initialize_database
from taxane_reaction_grammar_space.validate import validate_generated_space


class ValidateTests(unittest.TestCase):
    def test_empty_schema_passes_relational_checks(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            database = root / "space.sqlite"
            connection = _initialize_database(database)
            connection.close()
            outputs = validate_generated_space(database, root / "validation")
            summary = outputs["summary"].read_text(encoding="utf-8")
            self.assertIn('"validation_status": "pass"', summary)


if __name__ == "__main__":
    unittest.main()
