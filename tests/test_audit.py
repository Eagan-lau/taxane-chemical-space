import unittest

from taxane_reaction_grammar_space.audit import _contains_cjk


class AuditTests(unittest.TestCase):
    def test_cjk_detection(self):
        self.assertFalse(
            _contains_cjk(
                "G3 is reported as an exploratory chemical-space frontier."
            )
        )
        self.assertTrue(_contains_cjk("不得包含中文字符"))


if __name__ == "__main__":
    unittest.main()
