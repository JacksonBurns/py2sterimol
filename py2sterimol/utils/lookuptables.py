
"""
Usage:
First key is the type of atom, second key is the bonds on the atom.
Final resulting symbol is what sterimol uses to refer to that atom.

atom_and_bonds_to_symbol["C"][(2, 1, 1)] -> "C2" - ethylene carbon

"""
atom_and_bonds_to_symbol = {
    "C": {
        (1, 1, 1, 1): "C",  # tetrahedral carbon
        (4, 2, 1): "C2",  # ethylene carbon
        (3, 1): "C3",  # acetylene carbon
        (4, 2, 1): "C4",  # amide carbon
        (3, 3, 1): "C5",  # C in aromatic 5-ring
        (2, 2, 1): "C6",  # C in aromatic 6-ring
        (3, 2, 1): "C7",  # carbon linking 5- and 6-ring
        (2, 2, 1, 1): "C8",  # cyclopropane carbon
        (2, 2, 2): "C66",  # carbon linking two 6-rings
    },
    "N": {
        (1, 1, 1, 1): "N",  # tetrahedral nitrogen
        (5, 4, 1): "N4",  # amide nitrogen
        (3, 3, 1): "N5",  # N in aromatic 5-ring
        (2, 2, 1): "N6",  # N in aromatic 6-ring
    },
    "H": {
        (1, 5): "H",  # hydrogen
    },
    "O": {
        (1, 1): "O",  # normal oxygen
        (2): "O2",  # double bonded oxygen
    },
    "P": {
        (1, 1, 1, 1): "P",  # tetrahedral phosphorus
    },
    "S": {
        (1, 1): "S",  # normal sulfur
        (1, 1, 1, 1): "S4",  # tetrahedral sulfur
        (1, 1, 1, 1, 1, 1): "S1",  # octahedral sulfur
    },
    "F": {
        (1): "F",  # fluorine
    },
    "Cl": {
        (1): "C1",  # chlorine
    },
    "Br": {
        (1): "B1",  # bromine
    },
    "I": {
        (1): "I",  # iodine
    }
}

"""
Usage:
simple lookup for bond value to symbol for sterimol

"""
bond_value_to_bond_symbol = {
    1: "",
    2: "D",
    3: "T",
    4: "A",
    5: "E",
}
