import os
import sys
import unittest

from py2sterimol import py2sterimol


class Test_py2sterimol(unittest.TestCase):
    """
    Test the various functionalities of py2sterimol.
    """

    @classmethod
    def setUpClass(self):
        cwd = os.getcwd()
        self._tbu_pdb = os.path.join(cwd, "test", "data", "tbu.pdb")

    def test_tbu(self):
        """
        """
        p2s = py2sterimol(self._tbu_pdb)
        params = p2s(sterimolargs={})
        self.assertEqual(
            params["L"],
            4.11,
            "Value of L incorrect for tbu test."
        )
        self.assertEqual(
            params["B1"],
            2.77,
            "Value of B1 incorrect for tbu test."
        )
        self.assertEqual(
            params["B2"],
            2.98,
            "Value of B2 incorrect for tbu test."
        )
        self.assertEqual(
            params["B3"],
            3.16,
            "Value of B3 incorrect for tbu test."
        )
        self.assertEqual(
            params["B4"],
            3.15,
            "Value of B4 incorrect for tbu test."
        )
        self.assertEqual(
            params["B5"],
            3.17,
            "Value of B5 incorrect for tbu test."
        )


if __name__ == '__main__':
    unittest.main()
