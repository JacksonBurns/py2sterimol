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

    def test_py2sterimol(self):
        """
        """
        for molecule in molecule_names:
            with self.subTest(molecule=molecule):
                p2s = py2sterimol(self._tbu_pdb)
                results = p2s(sterimolargs={})

                for parameter in parameters:
                    with self.subTest(parameter=parameter):
                        self.assertEqual(
                            answer_lookup[molecule][parameter],
                            results[parameter],
                            "Value of {} incorrect for {}; expected {} got {}".format(
                                parameter,
                                molecule,
                                answer_lookup[molecule][parameter],
                                results[parameter],
                            )
                        )


molecule_names = [
    '1Nap',
    '35diMePh',
    '4ClPh',
    '4MeOPh',
    '4MePh',
    'Ad',
    'Bn',
    'CEt3',
    'CH2iPr',
    'CH2tBu',
    'CHEt2',
    'cHex',
    'CHiPr2',
    'CHPh2',
    'CHPr2',
    'Et',
    'H',
    'iPr',
    'Me',
    'nBu',
    'nPr',
    'Ph',
    'tBu',
]

parameters = [
    'L',
    'B1',
    'B2',
    'B3',
    'B4',
    'B5',
]

answer_lookup = {
    '1Nap': {"L": 6.28, "B1": 5.50, "B2": 3.11, "B3": 1.72, "B4": 1.71, "B5": 5.50},
    '35diMePh': {"L": 6.28, "B1": 4.26, "B2": 4.26, "B3": 2.05, "B4": 1.71, "B5": 4.30},
    '4ClPh': {"L": 7.74, "B1": 3.11, "B2": 3.11, "B3": 1.80, "B4": 1.80, "B5": 3.11},
    '4MeOPh': {"L": 8.20, "B1": 3.11, "B2": 3.11, "B3": 1.78, "B4": 2.02, "B5": 3.11},
    '4MePh': {"L": 7.22, "B1": 3.11, "B2": 3.11, "B3": 2.04, "B4": 1.71, "B5": 3.11},
    'Ad': {"L": 6.17, "B1": 3.49, "B2": 3.49, "B3": 3.16, "B4": 3.16, "B5": 3.49},
    'Bn': {"L": 4.62, "B1": 1.52, "B2": 6.02, "B3": 3.10, "B4": 3.13, "B5": 6.02},
    'CEt3': {"L": 5.05, "B1": 4.41, "B2": 3.17, "B3": 4.23, "B4": 2.77, "B5": 4.45},
    'CH2iPr': {"L": 5.05, "B1": 1.52, "B2": 4.21, "B3": 1.90, "B4": 3.16, "B5": 4.45},
    'CH2tBu': {"L": 5.05, "B1": 1.52, "B2": 4.22, "B3": 3.15, "B4": 3.16, "B5": 4.45},
    'CHEt2': {"L": 5.05, "B1": 4.41, "B2": 1.90, "B3": 3.49, "B4": 2.77, "B5": 4.45},
    'cHex': {"L": 6.17, "B1": 2.76, "B2": 3.49, "B3": 3.16, "B4": 1.91, "B5": 3.49},
    'CHiPr2': {"L": 5.05, "B1": 3.69, "B2": 2.05, "B3": 4.42, "B4": 4.41, "B5": 4.45},
    'CHPh2': {"L": 4.62, "B1": 1.99, "B2": 5.26, "B3": 1.50, "B4": 5.19, "B5": 6.02},
    'CHPr2': {"L": 6.17, "B1": 5.67, "B2": 1.90, "B3": 4.43, "B4": 2.79, "B5": 5.67},
    'Et': {"L": 4.11, "B1": 1.52, "B2": 2.98, "B3": 1.90, "B4": 1.91, "B5": 3.17},
    'H': {"L": 2.06, "B1": 1.00, "B2": 1.00, "B3": 1.00, "B4": 1.00, "B5": 1.00},
    'iPr': {"L": 4.11, "B1":  3.15, "B2": 1.90, "B3": 2.98, "B4": 2.77, "B5": 3.17},
    'Me': {"L": 3.00, "B1": 1.52, "B2": 2.04, "B3": 1.90, "B4": 1.90, "B5": 2.04},
    'nBu': {"L": 6.17, "B1":  1.52, "B2": 4.43, "B3": 1.90, "B4": 1.92, "B5": 4.54},
    'nPr': {"L": 5.05, "B1": 1.52, "B2": 3.49, "B3": 1.90, "B4": 1.91, "B5": 3.49},
    'Ph': {"L": 6.28, "B1": 3.11, "B2": 3.11, "B3": 1.71, "B4": 1.71, "B5": 3.11},
    'tBu': {"L": 4.11, "B1": 2.77, "B2": 2.98, "B3": 3.16, "B4": 3.15, "B5": 3.17},
}

if __name__ == '__main__':
    unittest.main()
