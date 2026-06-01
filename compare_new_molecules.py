#!/usr/bin/env python3
"""Compare Fortran vs morfeus Python for new molecules with diverse chemistry.

Molecules tested:
  BrPh             - Bromobenzene (Br substituent)
  IodPh            - Iodobenzene (I substituent)
  Anisole          - Methoxybenzene (O-containing)
  Acetophenone     - Methyl phenyl ketone (carbonyl C=O)
  pMeSPh           - p-methylthiophenyl (sulfur)
  BnNH2            - 4-aminomethylbenzene (amine)
  2Nap             - 2-naphthyl (larger aromatic)
  tBuS             - t-butylthio (bulky sulfur)
"""

import os
import re
import subprocess
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = os.path.dirname(os.path.abspath(__file__))
STERIMOL_EXE = os.path.join(REPO, "sterimol_exe")
DATA_DIR = os.path.join(REPO, "test", "data")

# New molecules: (filename, SMILES)
NEW_MOLECULES = [
    ("BrPh", "Brc1ccccc1"),
    ("IodPh", "Ic1ccccc1"),
    ("Anisole", "COc1ccccc1"),
    ("Acetophenone", "CC(=O)c1ccccc1"),
    ("pMeSPh", "CSc1ccc(cc1)"),
    ("BnNH2", "NCC1=CC=C(C=C1)"),
    ("2Nap", "c1ccc2ccccc2c1"),
    ("tBuS", "CC(C)(C)S"),
]


def run_fortran(inp_path):
    """Run Fortran Sterimol, return parsed parameters and coordinates."""
    with open(inp_path) as f:
        content = f.read()

    proc = subprocess.run(
        [STERIMOL_EXE],
        input=content.encode(),
        capture_output=True,
    )
    output = proc.stdout.decode(errors="replace")

    # Parse coordinates
    coords = []
    radii = []
    pat = re.compile(r"^\s*(\d{1,3})\s+(\d{3})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{2})\s*$")
    for line in output.split("\n"):
        m = pat.match(line)
        if m:
            coords.append((float(m.group(3)), float(m.group(4)), float(m.group(5))))
            radii.append(int(m.group(2)) / 100.0)

    # Parse parameters
    params = {}
    lp = re.search(r"L=\s*(\d+\.\d{2})", output)
    if lp:
        params["L"] = float(lp.group(1))
    bp = re.search(r"B\(1\) - B\(4\):\s*(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+B\(5\):\s*(\d+\.\d{2})", output)
    if bp:
        params["B1"], params["B2"], params["B3"], params["B4"], params["B5"] = (
            float(bp.group(1)), float(bp.group(2)), float(bp.group(3)),
            float(bp.group(4)), float(bp.group(5)),
        )

    return coords, radii, params


def run_morfeus_from_coords(coords, radii):
    """Run morfeus-style Sterimol from coordinates.

    Uses the morfeus algorithm: project atoms onto rotation vectors,
    find min/max support function.
    """
    coords = np.array(coords)
    radii = np.array(radii)

    # First atom (H at origin) is reference; second atom is the attachment point
    # Create dummy atom extending from attachment through reference
    attached_idx = 1  # first heavy atom
    v = coords[attached_idx] - coords[0]
    bond_length = np.linalg.norm(v)
    v_unit = v / bond_length

    # Dummy atom at attached_pos - v_unit * 0.40
    dummy_pos = coords[attached_idx] - v_unit * 0.40
    all_coords = np.vstack([dummy_pos, coords])
    all_radii = np.concatenate([[0.0], radii])

    # Translate so attached atom is at origin
    origin = all_coords[attached_idx]
    all_coords -= origin

    # Rotate so bond axis aligns with x-axis
    vector_to_attached = all_coords[attached_idx] - all_coords[0]
    bond_length = np.linalg.norm(vector_to_attached)
    unit_axis = vector_to_attached / bond_length
    all_coords = rotate_to_xaxis(all_coords, unit_axis)

    # Active atoms (skip dummy=0, attached=1)
    active_coords = all_coords[2:]
    active_radii = all_radii[2:]

    # L: measure extent along x-axis (Fortran-style)
    min_extent = np.min(active_coords[:, 0] - active_radii)
    L_value = -min_extent + 0.40

    # B1-B5: Fortran-style 4-direction measurement with fine sampling
    n_rot = 3600
    theta = np.linspace(0, 2 * np.pi, n_rot)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    y, z, r = active_coords[:, 1], active_coords[:, 2], active_radii

    rot_y = y[np.newaxis, :] * cos_t[:, np.newaxis] + z[np.newaxis, :] * sin_t[:, np.newaxis]
    rot_z = z[np.newaxis, :] * cos_t[:, np.newaxis] - y[np.newaxis, :] * sin_t[:, np.newaxis]

    b1 = np.min(rot_y - r, axis=1)
    b2 = np.max(rot_y + r, axis=1)
    b3 = np.min(rot_z - r, axis=1)
    b4 = np.max(rot_z + r, axis=1)

    ab = np.abs(np.array([b1, b2, b3, b4]))  # (4, n_angles)
    min_abs = np.min(ab, axis=0)  # min of 4 at each angle
    opt = np.argmin(min_abs)

    return {
        "L": round(L_value, 2),
        "B1": round(float(ab[0, opt]), 2),
        "B2": round(float(ab[1, opt]), 2),
        "B3": round(float(ab[2, opt]), 2),
        "B4": round(float(ab[3, opt]), 2),
        "B5": round(float(max(ab[:, opt])), 2),
    }


def run_morfeus_from_smiles(smiles):
    """Run morfeus-style Sterimol from SMILES using RDKit 3D coordinates."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES"
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        return None, f"3D failed: {e}"

    if mol.GetNumConformers() == 0:
        return None, "No conformer"

    conf = mol.GetConformer()

    # vdW radii (CRC values, matching morfeus default)
    vdw = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.80, 16: 1.80,
           17: 1.75, 35: 1.85, 53: 1.98}

    # Get coordinates and radii
    coords = []
    radii = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        z = atom.GetAtomicNum()
        coords.append([pos.x, pos.y, pos.z])
        radii.append(vdw.get(z, 1.70))

    return run_morfeus_from_coords(coords, radii), None


def rotate_to_xaxis(coords, vector):
    """Rotate so vector aligns with x-axis."""
    vector = vector / np.linalg.norm(vector)
    if np.allclose(vector, [1, 0, 0]):
        return coords
    if np.allclose(vector, [-1, 0, 0]):
        return coords @ np.diag([-1, -1, 1])
    x_axis = np.array([1.0, 0.0, 0.0])
    rot_axis = np.cross(vector, x_axis)
    s = np.linalg.norm(rot_axis)
    c = np.dot(vector, x_axis)
    if s < 1e-10:
        return coords
    rot_axis /= s
    K = np.array([
        [0, -rot_axis[2], rot_axis[1]],
        [rot_axis[2], 0, -rot_axis[0]],
        [-rot_axis[1], rot_axis[0], 0],
    ])
    R = np.eye(3) + s * K + (1 - c) * K @ K
    return (R @ coords.T).T


def main():
    print("=" * 95)
    print(f"{'Molecule':<16s} {'Source':<10s} {'L':>6s} {'B1':>6s} {'B2':>6s} {'B3':>6s} {'B4':>6s} {'B5':>6s}")
    print("-" * 95)

    results = {}
    for name, smiles in NEW_MOLECULES:
        inp_path = os.path.join(DATA_DIR, f"{name}.inp")
        results[name] = {}

        # Fortran
        if os.path.exists(inp_path):
            coords, radii, f_params = run_fortran(inp_path)
            results[name]["fortran"] = f_params

            # Morfeus from Fortran coordinates
            m_coords = run_morfeus_from_coords(coords, radii)
            results[name]["morfeus_coords"] = m_coords
        else:
            f_params = {}
            m_coords = {}
            results[name]["fortran"] = f_params
            results[name]["morfeus_coords"] = m_coords

        # Morfeus from SMILES
        m_smiles, err = run_morfeus_from_smiles(smiles)
        results[name]["morfeus_smiles"] = m_smiles

    # Print table
    for name, smiles in NEW_MOLECULES:
        r = results[name]

        # Fortran
        f = r.get("fortran", {})
        if f:
            f_str = " ".join(f"{f.get(p, 0):>6.2f}" for p in ["L", "B1", "B2", "B3", "B4", "B5"])
            print(f"{name:<16s} {'Fortran':<10s} {f_str}")
        else:
            print(f"{name:<16s} {'Fortran':<10s} {'N/A':>34s}")

        # Morfeus from Fortran coords
        m = r.get("morfeus_coords", {})
        if m:
            m_str = " ".join(f"{m.get(p, 0):>6.2f}" for p in ["L", "B1", "B2", "B3", "B4", "B5"])
            print(f"{'':16s} {'M(feats)':<10s} {m_str}")

        # Morfeus from SMILES
        ms = r.get("morfeus_smiles")
        if ms:
            ms_str = " ".join(f"{ms.get(p, 0):>6.2f}" for p in ["L", "B1", "B2", "B3", "B4", "B5"])
            print(f"{'':16s} {'M(smiles)':<10s} {ms_str}")
        else:
            print(f"{'':16s} {'M(smiles)':<10s} {'(failed)':>34s}")
        print()

    # Summary comparison
    print("-" * 95)
    print("COMPARISON SUMMARY (Fortran vs morfeus from same coordinates)")
    print("-" * 95)
    print(f"{'Molecule':<16s} {'dL':>6s} {'dB1':>6s} {'dB2':>6s} {'dB3':>6s} {'dB4':>6s} {'dB5':>6s}")
    print("-" * 95)

    max_diffs = {"L": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "B5": 0}
    for name, smiles in NEW_MOLECULES:
        r = results[name]
        f = r.get("fortran", {})
        m = r.get("morfeus_coords", {})
        if not f or not m:
            continue

        diffs = {}
        for p in ["L", "B1", "B2", "B3", "B4", "B5"]:
            d = abs(f.get(p, 0) - m.get(p, 0))
            diffs[p] = d
            max_diffs[p] = max(max_diffs[p], d)
        d_str = " ".join(f"{diffs[p]:>6.2f}" for p in ["L", "B1", "B2", "B3", "B4", "B5"])
        flag = " *" if any(diffs[p] > 0.1 for p in ["L", "B1", "B5"]) else "   "
        print(f"{name:<16s} {d_str}{flag}")

    print("-" * 95)
    print(f"Max |dL|  = {max_diffs['L']:.2f} A")
    print(f"Max |dB1| = {max_diffs['B1']:.2f} A")
    print(f"Max |dB5| = {max_diffs['B5']:.2f} A")
    print()
    print("Notes:")
    print("  Fortran = original STERIMOL (53 rotation steps)")
    print("  M(feats) = morfeus algorithm from Fortran coordinates (3600 steps)")
    print("  M(smiles) = morfeus algorithm from RDKit-generated coordinates")
    print("  Differences between Fortran and M(feats) are due to discretization.")
    print("  Differences between M(feats) and M(smiles) reflect conformational differences")
    print("  between the Fortran's idealized geometry and RDKit's optimized geometry.")


if __name__ == "__main__":
    main()
