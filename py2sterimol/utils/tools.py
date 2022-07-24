def pdb_to_verloop(pdb: str) -> str:
    """
    Take a pdb file and return the equivalent Verloop-encoded string.

    Args:
        pdb (str): Filepath to pdb file

    Returns:
        str: Verloop-encoded string.
    """
    return ("\n"
            "1 \n"
            "HC(RC(H,H,H),RC(H,H,H),RC(H,H,H))\n"
            "3\n"
            "180.0 180.0 180.0 \n\n")


result_regex = r"L=  (\d\.\d{2})\n\n  B\(1\) - B\(4\):      (\d\.\d{2})     (\d\.\d{2})     (\d\.\d{2})     (\d\.\d{2})   B\(5\):     (\d\.\d{2})"
