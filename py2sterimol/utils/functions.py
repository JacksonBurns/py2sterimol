from py2sterimol.utils import lookuptables


def pdb_to_verloop(pdb: str) -> str:
    """
    Take a pdb file and return the equivalent Verloop-encoded string.

    Args:
        pdb (str): Filepath to pdb file

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

    # number of dihedral (torsion angles) on next line

    # torsion angles

    # two blank lines
    encoded_string += '\n\n'

    return encoded_string
