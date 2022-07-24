import re

from py2sterimol.utils import lookuptables

from rdkit import Chem


def smiles_and_pdb_to_verloop(smiles: str, pdb: str) -> str:
    """
    Take a smiles string and pdb file and return the equivalent
    Verloop-encoded input card.

    Args:
        smiles (str): SMILES string
        pdb (str): Filepath to pdb

    Returns:
        str: Verloop-encoded string.

    Example:    

    (       "\n"
            "1 \n"
            "HC(RC(H,H,H),RC(H,H,H),RC(H,H,H))\n"
            "3\n"
            "180.0 180.0 180.0 \n\n")    
    """
    encoded_string = ""

    # start with a newline
    encoded_string += '\n'

    # output mode - always do 1
    encoded_string += '1\n'

    # formula for the molecule
    encoded_string += _smiles_to_verloop(smiles) + '\n'

    # number of dihedral (torsion angles) on next line, torsion angles on next
    encoded_string += _get_torsion_string(pdb)

    # two blank lines
    encoded_string += '\n\n'

    return encoded_string


def _smiles_to_verloop(smiles: str):
    """Turn a valid SMILES string into the corresponding Verloop encoding

    Args:
        smiles (str): SMILES string
    """
    verloop = ""

    print(smiles)

    # add explicit H's
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    smiles = Chem.MolToSmiles(mol)
    smiles = smiles.replace("[", "").replace("]", "")
    print(smiles)

    # replace H's with Verloop notation
    smiles = smiles.replace(
        "HC(H)(H)", "C(H,H,H)"
    ).replace(
        "C(H)(H)H", "C(H,H,H)"
    ).replace(
        "C(H)(H)", "C(H,H)"
    )

    print(smiles)

    # for character in smiles:
    #     if character in ('(', ')'):  # side chain
    #         continue
    #     else:  # atom
    #         if character == 'H':
    #             continue
    #         elif character == 'C':
    #             # find what types of bonds it has
    #             continue
    verloop += smiles
    return verloop


def _get_torsion_string(pdb: str):
    """Take pdb and find the torsion angles needed for sterimol

    Args:
        pdb (str): Filepath to PDB

    Returns:
        String formatted as shown below:

        # of torsion angles\n
        torsion angles

        3
        180.0 80.0 90.0
    """
    torsion_string = ""

    return torsion_string
