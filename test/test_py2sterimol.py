import os
import unittest

from py2sterimol import py2sterimol
from py2sterimol.utils.functions import _smiles_to_verloop


class Test_py2sterimol(unittest.TestCase):
    """
    Test py2sterimol against the original Fortran Sterimol program.
    """

    @classmethod
    def setUpClass(cls):
        cls._data_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "data"
        )
        cls._sterimol_exe = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "sterimol_exe",
        )

    def test_from_inp_file(self):
        """Test all .inp molecules against expected values."""
        for molecule in molecule_names:
            with self.subTest(molecule=molecule):
                inp_path = os.path.join(self._data_dir, f"{molecule}.inp")
                p2s = py2sterimol.from_inp_file(
                    inp_path, path_to_sterimol=self._sterimol_exe
                )
                results = p2s()

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

    def test_attached_atom_idx(self):
        """Test attached_atom_idx: different attachment points give different results."""
        # n-propyl: terminal carbons (0, 2) should be symmetric,
        # middle carbon (1) should be different
        r0 = py2sterimol(
            "CCC", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=0,
        )()
        r1 = py2sterimol(
            "CCC", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=1,
        )()
        r2 = py2sterimol(
            "CCC", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=2,
        )()

        # Terminals are symmetric
        self.assertAlmostEqual(r0["L"], r2["L"], places=1)
        self.assertAlmostEqual(r0["B1"], r2["B1"], places=1)
        self.assertAlmostEqual(r0["B5"], r2["B5"], places=1)

        # Middle is shorter (only methyl group on other side)
        self.assertAlmostEqual(r1["L"], 3.00, places=1)
        self.assertLess(r1["L"], r0["L"])

    def test_attached_atom_idx_tbutyl(self):
        """Test t-butyl: terminal vs quaternary carbon attachment."""
        # Terminal C (idx 0): full t-butyl shape
        r_term = py2sterimol(
            "CC(C)(C)", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=0,
        )()
        # Quaternary C (idx 1): three methyls radiating, shorter L
        r_quad = py2sterimol(
            "CC(C)(C)", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=1,
        )()

        # Terminal should have longer L
        self.assertAlmostEqual(r_term["L"], 4.11, places=1)
        self.assertAlmostEqual(r_quad["L"], 3.00, places=1)
        self.assertGreater(r_term["L"], r_quad["L"])

    def test_attached_atom_idx_out_of_range(self):
        """Test that out-of-range attached_atom_idx raises ValueError."""
        with self.assertRaises(ValueError):
            _smiles_to_verloop("CC", attached_atom_idx=2)

        with self.assertRaises(ValueError):
            _smiles_to_verloop("CC", attached_atom_idx=-1)

    def test_attached_atom_idx_default(self):
        """Test that None attached_atom_idx uses first heavy atom (idx 0)."""
        r_default = py2sterimol(
            "CCC", "", path_to_sterimol=self._sterimol_exe,
        )()
        r_zero = py2sterimol(
            "CCC", "", path_to_sterimol=self._sterimol_exe, attached_atom_idx=0,
        )()

        self.assertAlmostEqual(r_default["L"], r_zero["L"], places=2)
        self.assertAlmostEqual(r_default["B1"], r_zero["B1"], places=2)
        self.assertAlmostEqual(r_default["B5"], r_zero["B5"], places=2)


molecule_names = [
    '1Nap',
    '2Nap',
    '35diMePh',
    '4ClPh',
    '4MeOPh',
    '4MePh',
    'Acetophenone',
    'Ad',
    'Anisole',
    'Bn',
    'BnNH2',
    'BrPh',
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
    'IodPh',
    'Me',
    'nBu',
    'nPr',
    'pMeSPh',
    'Ph',
    'tBu',
    'tBuS',
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
    '2Nap': {"L": 6.28, "B1": 5.50, "B2": 3.11, "B3": 1.72, "B4": 1.71, "B5": 5.50},
    '35diMePh': {"L": 6.28, "B1": 4.26, "B2": 4.26, "B3": 2.05, "B4": 1.71, "B5": 4.30},
    '4ClPh': {"L": 7.74, "B1": 3.11, "B2": 3.11, "B3": 1.80, "B4": 1.80, "B5": 3.11},
    '4MeOPh': {"L": 8.20, "B1": 3.11, "B2": 3.11, "B3": 1.78, "B4": 2.02, "B5": 3.11},
    '4MePh': {"L": 7.22, "B1": 3.11, "B2": 3.11, "B3": 2.04, "B4": 1.71, "B5": 3.11},
    'Acetophenone': {"L": 5.92, "B1": 3.11, "B2": 3.11, "B3": 1.71, "B4": 1.71, "B5": 3.11},
    'Ad': {"L": 6.17, "B1": 3.49, "B2": 3.49, "B3": 3.16, "B4": 3.16, "B5": 3.49},
    'Anisole': {"L": 8.20, "B1": 3.11, "B2": 3.11, "B3": 1.78, "B4": 2.02, "B5": 3.11},
    'Bn': {"L": 4.62, "B1": 1.52, "B2": 6.02, "B3": 3.10, "B4": 3.13, "B5": 6.02},
    'BnNH2': {"L": 5.30, "B1": 1.52, "B2": 6.86, "B3": 3.10, "B4": 3.13, "B5": 6.86},
    'BrPh': {"L": 8.04, "B1": 3.07, "B2": 3.07, "B3": 1.95, "B4": 1.95, "B5": 3.11},
    'CEt3': {"L": 5.05, "B1": 4.41, "B2": 3.17, "B3": 4.23, "B4": 2.77, "B5": 4.45},
    'CH2iPr': {"L": 5.05, "B1": 1.52, "B2": 4.21, "B3": 1.90, "B4": 3.16, "B5": 4.45},
    'CH2tBu': {"L": 5.05, "B1": 1.52, "B2": 4.22, "B3": 3.15, "B4": 3.16, "B5": 4.45},
    'CHEt2': {"L": 5.05, "B1": 4.41, "B2": 1.90, "B3": 3.49, "B4": 2.77, "B5": 4.45},
    'cHex': {"L": 6.17, "B1": 2.76, "B2": 3.49, "B3": 3.16, "B4": 1.91, "B5": 3.49},
    'CHiPr2': {"L": 5.05, "B1": 3.69, "B2": 2.05, "B3": 4.42, "B4": 4.41, "B5": 4.45},
    'CHPh2': {"L": 4.62, "B1": 1.50, "B2": 6.01, "B3": 3.21, "B4": 3.02, "B5": 6.02},
    'CHPr2': {"L": 6.17, "B1": 5.67, "B2": 1.90, "B3": 4.43, "B4": 2.79, "B5": 5.67},
    'Et': {"L": 4.11, "B1": 1.52, "B2": 2.98, "B3": 1.90, "B4": 1.91, "B5": 3.17},
    'H': {"L": 2.06, "B1": 1.00, "B2": 1.00, "B3": 1.00, "B4": 1.00, "B5": 1.00},
    'iPr': {"L": 4.11, "B1": 3.15, "B2": 1.90, "B3": 2.98, "B4": 2.77, "B5": 3.17},
    'IodPh': {"L": 8.45, "B1": 2.96, "B2": 2.96, "B3": 2.15, "B4": 2.15, "B5": 3.11},
    'Me': {"L": 3.00, "B1": 1.52, "B2": 2.04, "B3": 1.90, "B4": 1.90, "B5": 2.04},
    'nBu': {"L": 6.17, "B1": 1.52, "B2": 4.43, "B3": 1.90, "B4": 1.92, "B5": 4.54},
    'nPr': {"L": 5.05, "B1": 1.52, "B2": 3.49, "B3": 1.90, "B4": 1.91, "B5": 3.49},
    'pMeSPh': {"L": 8.02, "B1": 3.11, "B2": 3.36, "B3": 1.78, "B4": 2.03, "B5": 3.36},
    'Ph': {"L": 6.28, "B1": 3.11, "B2": 3.11, "B3": 1.71, "B4": 1.71, "B5": 3.11},
    'tBu': {"L": 4.11, "B1": 2.77, "B2": 2.98, "B3": 3.16, "B4": 3.15, "B5": 3.17},
    'tBuS': {"L": 4.95, "B1": 1.70, "B2": 4.15, "B3": 1.70, "B4": 4.07, "B5": 4.67},
}

if __name__ == '__main__':
    unittest.main()
