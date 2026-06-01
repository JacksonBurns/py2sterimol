import os
import re
import math

from py2sterimol.utils import lookuptables

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry


def smiles_and_pdb_to_verloop(
    smiles: str,
    pdb: str,
    attached_atom_idx=None,
) -> str:
    """
    Take a SMILES string and PDB file and return the equivalent
    Verloop-encoded input card.

    Args:
        smiles: SMILES string.
        pdb: Filepath to PDB file with 3D coordinates (or empty string).
        attached_atom_idx: 0-based index of the heavy atom to use as the
                           attachment point. If None, defaults to the first
                           heavy atom in SMILES order.

    Returns:
        Verloop-encoded string ready for the Fortran Sterimol program.
    """
    encoded_string = ""

    # Blank line (old symbol table position, now embedded in code)
    encoded_string += '\n'

    # Output mode - always 1 (numeric output only)
    encoded_string += '1\n'

    # Formula for the molecule
    encoded_string += _smiles_to_verloop(
        smiles, attached_atom_idx=attached_atom_idx,
    ) + '\n'

    # Number of dihedral (torsion angles) on next line, torsion angles on next
    encoded_string += _get_torsion_string(smiles, pdb)

    # Two blank lines to terminate input
    encoded_string += '\n\n'

    return encoded_string


def _get_sterimol_atom_symbol(atom, mol):
    """Determine the Sterimol atom symbol for an RDKit atom.

    Maps RDKit atom properties to Verloop's atomic symbol encoding.
    """
    sym = atom.GetSymbol()
    atomic_num = atom.GetAtomicNum()

    # Count non-H neighbors and their bond orders
    neighbors = []
    for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        if nbr.GetAtomicNum() != 1:  # Skip hydrogens for connectivity
            neighbors.append((nbr, bond))

    # Determine specific atom type
    if sym == 'C':
        if atom.GetIsAromatic():
            ring_info = mol.GetRingInfo()
            if ring_info.IsAtomInRingOfSize(atom.GetIdx(), 5):
                return 'C5'
            elif ring_info.IsAtomInRingOfSize(atom.GetIdx(), 6):
                return 'C6'
            else:
                return 'C6'  # Default for aromatic C
        elif len(neighbors) == 2:
            bond_orders = [b.GetBondTypeAsDouble() for _, b in neighbors]
            if any(o == 3.0 for o in bond_orders):
                return 'C3'  # acetylene
            if any(o == 2.0 for o in bond_orders):
                return 'C2'  # ethylene (double bond)
            return 'C'  # sp3 with 2 heavy neighbors (e.g., ethane middle C)
        elif len(neighbors) == 3:
            if any(b.GetBondTypeAsDouble() == 2.0 for _, b in neighbors):
                return 'C4'  # has double bond (amide/carbonyl)
            return 'C'  # sp3 carbon
        elif len(neighbors) == 4:
            return 'C'  # quaternary carbon
        else:
            return 'C'
    elif sym == 'N':
        if atom.GetIsAromatic():
            return 'N6'
        elif any(b.GetBondTypeAsDouble() == 2.0 for _, b in neighbors):
            return 'N4'  # amide-like
        return 'N'
    elif sym == 'O':
        if any(b.GetBondTypeAsDouble() == 2.0 for _, b in neighbors):
            return 'O2'
        return 'O'
    elif sym == 'S':
        bond_count = atom.GetTotalValence()
        if bond_count >= 6:
            return 'S1'
        elif bond_count >= 4:
            return 'S4'
        return 'S'
    elif sym == 'P':
        return 'P'
    elif sym == 'F':
        return 'F'
    elif sym == 'Cl':
        return 'C1'
    elif sym == 'Br':
        return 'B1'
    elif sym == 'I':
        return 'I'
    elif sym == 'H':
        return 'H'
    else:
        return sym  # Fallback: use element symbol


def _bond_symbol(bond):
    """Return the Verloop bond type symbol for a bond."""
    bt = bond.GetBondTypeAsDouble()
    if bt == 2.0:
        return 'D'
    elif bt == 3.0:
        return 'T'
    return ''  # Single bond: no symbol


def _smiles_to_verloop(smiles: str, attached_atom_idx=None) -> str:
    """Turn a valid SMILES string into the corresponding Verloop encoding.

    Uses RDKit to parse the molecular graph and generate the Verloop notation
    used by the original Fortran Sterimol program.

    Args:
        smiles: SMILES string.
        attached_atom_idx: 0-based index into the list of heavy atoms to use
                           as the attachment (root) atom. If None, defaults
                           to 0 (first heavy atom in SMILES order).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # Collect heavy atom indices in SMILES order
    heavy_indices = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1
    ]
    if not heavy_indices:
        raise ValueError("SMILES contains no heavy atoms")

    # Determine attachment atom
    if attached_atom_idx is not None:
        if attached_atom_idx < 0 or attached_atom_idx >= len(heavy_indices):
            raise ValueError(
                f"attached_atom_idx={attached_atom_idx} out of range "
                f"(0..{len(heavy_indices) - 1}) for {len(heavy_indices)} "
                f"heavy atoms"
            )
        root_idx = heavy_indices[attached_atom_idx]
    else:
        root_idx = heavy_indices[0]

    # Build adjacency for traversal
    visited = set()
    visited.add(root_idx)

    root_atom = mol.GetAtomWithIdx(root_idx)
    atom_sym = _get_sterimol_atom_symbol(root_atom, mol)

    # Count H neighbors (excluding one for the root attachment H)
    h_count = 0
    for bond in root_atom.GetBonds():
        nbr = bond.GetOtherAtom(root_atom)
        if nbr.GetAtomicNum() == 1:
            h_count += 1
    h_count = max(0, h_count - 1)

    # Build root atom with H children
    root_str = atom_sym
    if h_count > 0:
        root_str += f"({','.join(['H'] * h_count)})"

    verloop = "H" + root_str

    # Process heavy atom neighbors as main chain extensions
    for bond in root_atom.GetBonds():
        nbr = bond.GetOtherAtom(root_atom)
        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() not in visited:
            child_str = _traverse_atom(nbr, bond, mol, visited)
            if child_str:
                verloop += "R" + child_str

    return verloop


def _traverse_atom(atom, parent_bond, mol, visited):
    """Recursively traverse the molecular graph, building Verloop notation.

    Returns the Verloop string for this atom and all its descendants.
    """
    if atom.GetIdx() in visited:
        return None  # Already visited (ring closure)

    visited.add(atom.GetIdx())

    sym = _get_sterimol_atom_symbol(atom, mol)
    bond_sym = _bond_symbol(parent_bond)

    # Separate heavy atom children and hydrogen children
    heavy_children = []
    h_children = []

    for bond in atom.GetBonds():
        nbr = bond.GetOtherAtom(atom)
        if nbr.GetAtomicNum() == 1:
            # Hydrogen child
            if nbr.GetIdx() not in visited:
                h_children.append("H")
                visited.add(nbr.GetIdx())
        elif nbr.GetIdx() not in visited:
            # Heavy atom child - becomes a side chain in parentheses
            child_str = _traverse_atom(nbr, bond, mol, visited)
            if child_str:
                heavy_children.append(child_str)

    # Build children list: H's first, then heavy atoms (as side chains)
    children = h_children + heavy_children

    # Build the string
    result = bond_sym + sym
    if children:
        result += f"({','.join(children)})"

    return result


def _get_torsion_string(smiles: str, pdb: str) -> str:
    """Compute torsion angles from a 3D structure and return formatted string.

    Args:
        smiles: SMILES string for the molecule
        pdb: Filepath to PDB file with 3D coordinates

    Returns:
        Formatted string with torsion count and values:
        e.g., "3\\n180.0 90.0 180.0 "
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "0\n"

    mol = Chem.AddHs(mol)

    # Try to load coordinates from PDB
    if pdb and os.path.exists(pdb):
        try:
            mol_with_coords = Chem.MolFromPDBFile(pdb)
            if mol_with_coords is not None and mol_with_coords.GetNumConformers() > 0:
                mol = mol_with_coords
        except Exception:
            pass

    # If no valid 3D coords from PDB, generate from SMILES
    if mol.GetNumConformers() == 0:
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            if mol.GetNumConformers() > 0:
                AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass

    if mol.GetNumConformers() == 0:
        # If 3D embedding fails, use default torsion angles
        return _default_torsion_string(smiles)

    conformer = mol.GetConformer()

    # Find torsion angles: for each sequence of 4 consecutive heavy atoms
    heavy_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]

    torsions = []
    for i in range(len(heavy_indices) - 3):
        a1 = heavy_indices[i]
        a2 = heavy_indices[i + 1]
        a3 = heavy_indices[i + 2]
        a4 = heavy_indices[i + 3]

        try:
            p1 = conformer.GetAtomPosition(a1)
            p2 = conformer.GetAtomPosition(a2)
            p3 = conformer.GetAtomPosition(a3)
            p4 = conformer.GetAtomPosition(a4)

            torsion = _compute_dihedral(p1, p2, p3, p4)
            torsions.append(torsion)
        except Exception:
            torsions.append(180.0)  # Default

    if not torsions:
        return "0\n"

    torsion_str = f"{len(torsions)}\n"
    torsion_str += " ".join(f"{t:.1f}" for t in torsions) + " "
    return torsion_str


def _compute_dihedral(p1, p2, p3, p4):
    """Compute the dihedral angle between four points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = b1.Cross(b2)
    n2 = b2.Cross(b3)

    n1.Normalize()
    n2.Normalize()

    m1 = n1.Cross(b2.Normalize())
    m1.Normalize()

    x = n1.Dot(n2)
    y = m1.Dot(n2)

    angle = math.degrees(math.atan2(y, x))
    return angle


def _default_torsion_string(smiles: str) -> str:
    """Generate default torsion angles based on SMILES analysis.

    For simple molecules, uses 180.0 (anti) as default.
    For branched molecules, uses staggered conformations.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "0\n"

    # Count heavy atoms that need torsion angles
    heavy_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    n_torsions = max(0, heavy_count - 3)

    if n_torsions == 0:
        return "0\n"

    # Use 180.0 as default (anti conformation) for all torsions
    return f"{n_torsions}\n" + " ".join(["180.0"] * n_torsions) + " "
