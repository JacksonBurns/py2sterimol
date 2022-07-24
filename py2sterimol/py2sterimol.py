import os
import subprocess
import sys
import re
from typing import Dict

from py2sterimol.utils.tools import (
    pdb_to_verloop,
    result_regex,
)


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
        self._write_input()
        self._call_sterimol(sterimolargs)
        self._parse_output(results)
        return results

    def _write_input(self):
        encoded = pdb_to_verloop(self._pdb_fpath)
        with open('py2sterimol_input.inp', 'w') as file:
            file.writelines(encoded)

    def _call_sterimol(self, sterimolargs):
        """Uses subprocess to call sterimol

        Args:
            sterimolargs (dict): Extra arguments for sterimol behavior.
        """
        # find the appropriate command based on the platform
        if sys.platform == "win32":
            cmd = 'type'
        else:
            cmd = 'cat'
        # send file to stdout, piping
        cat_proc = subprocess.Popen(
            [cmd, 'py2sterimol_input.inp'],
            stdout=subprocess.PIPE,
            shell=True,
        )
        # call sterimol with stdout piped to stdin, save new stdout to file
        with open("py2sterimol_input.out", 'w') as file:
            subprocess.run(
                self._path_to_sterimol,
                stdin=cat_proc.stdout,
                stdout=file,
            )

    def _parse_output(self, results):
        """Open the output file, retrieve parameters with regex and update results dictionary.

        Args:
            results (dict): Dictionary with sterimol parameters as keys.

        Returns:
            dict: Updated results
        """
        with open("py2sterimol_input.out", 'r') as file:
            file_result = re.search(result_regex, file.read())
            for param, i in zip(results.keys(), range(1, 7)):
                results[param] = float(file_result.group(i))
        return results
