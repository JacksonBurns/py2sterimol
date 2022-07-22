from typing import Dict


class py2sterimol:
    def __init__(
        self,
        pdb,
    ) -> None:
        self._pdb = pdb

    def __call__(self, sterimolargs={}) -> Dict[str, int or None]:
        results = {
            'L': None,
            'B1': None,
            'B2': None,
            'B3': None,
            'B4': None,
            'B5': None,
        }
        self._write_input()
        self._call_sterimol()
        self._parse_output(results)
        return results

    def _write_input(self):
        pass

    def _call_sterimol(self):
        pass

    def _parse_output(self, results):
        # regex to find results
        # update results for each
        return results
