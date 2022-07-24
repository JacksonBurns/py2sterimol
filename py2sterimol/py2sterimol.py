import os
import subprocess
import sys
import re
from typing import Dict

from py2sterimol.utils.functions import pdb_to_verloop
from py2sterimol.utils.strings import result_regex


class py2sterimol:
    def __init__(
        self,
        pdb_fpath: str,
        path_to_sterimol: str = 'sterimol.exe',
    ) -> None:
        self._pdb_fpath = pdb_fpath
        self._path_to_sterimol = path_to_sterimol

    def __call__(self, sterimolargs: dict = {}) -> Dict[str, int or None]:
        """Encode and write input before calling sterimol and parsing results.

        Args:
            sterimolargs (dict, optional): Optional args for sterimol. Defaults to {}.

        Returns:
            Dict[str, int or None]: Sterimol parameters.
        """
        results = {
            'L': None,
            'B1': None,
            'B2': None,
            'B3': None,
            'B4': None,
            'B5': None,
        }
        encoded = pdb_to_verloop(self._pdb_fpath)

        # send encoded to stdout, piping
        cat_proc = subprocess.Popen(
            ['echo', encoded],
            stdout=subprocess.PIPE,
            shell=True,
        )

        # call sterimol with stdout piped to stdin, save new stdout to variable
        sterimol_proc = subprocess.Popen(
            self._path_to_sterimol,
            stdin=cat_proc.stdout,
            stdout=subprocess.PIPE,
        )
        out, _ = sterimol_proc.communicate()

        # update results dict
        file_result = re.search(result_regex, out)
        for param, i in zip(results.keys(), range(1, 7)):
            results[param] = float(file_result.group(i))

        return results
