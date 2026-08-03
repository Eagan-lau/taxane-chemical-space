import unittest

from taxane_reaction_grammar_space.screen import _bit_vector_as_int


class ScreenTests(unittest.TestCase):
    def test_bit_vector_conversion(self):
        from rdkit import DataStructs

        vector = DataStructs.ExplicitBitVect(8)
        vector.SetBit(0)
        vector.SetBit(7)
        self.assertEqual(_bit_vector_as_int(vector), int("10000001", 2))


if __name__ == "__main__":
    unittest.main()
