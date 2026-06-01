import subprocess
import re
from typing import Dict, Optional

from py2sterimol.utils.functions import smiles_and_pdb_to_verloop
from py2sterimol.utils.strings import result_regex


class py2sterimol:
    """Wrapper around the original Fortran Sterimol program.

    Can be initialized from SMILES + PDB, or directly from a Verloop-formatted
    .inp file or verloop string.
    """

    def __init__(
        self,
        smiles: str,
        pdb_fpath: str,
        path_to_sterimol: str = 'sterimol.exe',
        attached_atom_idx: Optional[int] = None,
    ) -> None:
        """Create a py2sterimol instance from SMILES and a PDB file.

        Args:
            smiles: SMILES string for the molecule.
            pdb_fpath: Path to PDB file with 3D coordinates (or empty string
                       to generate coordinates from SMILES).
            path_to_sterimol: Path to the compiled Fortran sterimol executable.
            attached_atom_idx: 0-based index of the heavy atom to use as the
                               attachment point (i.e., the atom closest to the
                               bond axis origin). If None, defaults to the
                               first heavy atom in SMILES order. This is the
                               equivalent of morfeus's ``attached_index``.
        """
        self._smiles = smiles
        self._pdb_fpath = pdb_fpath
        self._path_to_sterimol = path_to_sterimol
        self._verloop_input = None
        self._attached_atom_idx = attached_atom_idx

    @classmethod
    def from_verloop_string(
        cls,
        verloop_string: str,
        path_to_sterimol: str = 'sterimol.exe',
    ) -> 'py2sterimol':
        """Create instance from a raw Verloop-format input string.

        Args:
            verloop_string: Complete Sterimol input deck (blanks, IPR, formula,
                           torsion count, torsion angles, trailing blanks).
            path_to_sterimol: Path to compiled sterimol executable.
        """
        instance = cls.__new__(cls)
        instance._smiles = None
        instance._pdb_fpath = None
        instance._path_to_sterimol = path_to_sterimol
        instance._verloop_input = verloop_string
        instance._attached_atom_idx = None
        return instance

    @classmethod
    def from_inp_file(
        cls,
        inp_fpath: str,
        path_to_sterimol: str = 'sterimol.exe',
    ) -> 'py2sterimol':
        """Create instance by reading a Verloop-format .inp file.

        Args:
            inp_fpath: Path to .inp file.
            path_to_sterimol: Path to compiled sterimol executable.
        """
        with open(inp_fpath, 'r') as f:
            verloop_string = f.read()
        return cls.from_verloop_string(verloop_string, path_to_sterimol)

    def __call__(self, sterimolargs: dict = {}) -> Dict[str, float]:
        """Run the Fortran Sterimol program and return parsed results.

        Returns:
            Dict with keys L, B1, B2, B3, B4, B5.
        """
        results = {
            'L': None,
            'B1': None,
            'B2': None,
            'B3': None,
            'B4': None,
            'B5': None,
        }

        if self._verloop_input:
            encoded = self._verloop_input
        else:
            encoded = smiles_and_pdb_to_verloop(
                self._smiles, self._pdb_fpath,
                attached_atom_idx=self._attached_atom_idx,
            )

        encoded_bytes = encoded.encode('ascii')

        sterimol_proc = subprocess.Popen(
            [self._path_to_sterimol],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = sterimol_proc.communicate(input=encoded_bytes)

        file_result = re.search(result_regex, out)
        if not file_result:
            raise RuntimeError(
                f"Failed to parse Sterimol output.\n"
                f"Stdout: {out.decode('ascii', errors='replace')}\n"
                f"Stderr: {err.decode('ascii', errors='replace')}"
            )
        for param, i in zip(results.keys(), range(1, 7)):
            results[param] = float(file_result.group(i))

        return results
