#!/usr/bin/env python3
"""Compare Fortran Sterimol vs Python (morfeus) Sterimol implementations.

This script:
1. Runs the Fortran Sterimol program on each .inp file
2. Parses the Fortran output to get atom coordinates and types
3. Runs the Python morfeus Sterimol on the same coordinates
4. Compares L, B1, B5 values between the two implementations
"""

import os
import re
import subprocess
import sys
import numpy as np

# Add repo root to path
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

# --- Configuration ---
STERIMOL_EXE = os.path.join(REPO_ROOT, "sterimol_exe")
DATA_DIR = os.path.join(REPO_ROOT, "test", "data")


def parse_inp_file(inp_path):
    """Parse a Verloop .inp file to extract element information."""
    with open(inp_path, 'r') as f:
        lines = [l.rstrip() for l in f.readlines()]

    # Find the formula line (after the IPR line)
    i = 0
    while i < len(lines) and lines[i].strip() == '':
        i += 1
    if i < len(lines):
        i += 1  # Skip IPR

    # Read formula (may span multiple lines with &)
    formula = ''
    while i < len(lines):
        line = lines[i].rstrip()
        if line.endswith('&'):
            formula += line[:-1]
            i += 1
        else:
            formula += line
            i += 1
            break

    # Count atoms in formula by element
    elements = count_atoms_in_formula(formula)
    return elements


def count_atoms_in_formula(formula):
    """Count atoms of each type in a Verloop formula string.

    Returns a list of atomic numbers in the order they appear.
    """
    # Symbol mapping from Verloop notation to atomic numbers
    symbol_map = {
        'H': 1,
        'C': 6, 'C2': 6, 'C3': 6, 'C4': 6, 'C5': 6, 'C6': 6, 'C7': 6, 'C8': 6, 'C66': 6,
        'N': 7, 'N4': 7, 'N5': 7, 'N6': 7,
        'O': 8, 'O2': 8,
        'F': 9,
        'P': 15,
        'S': 16, 'S1': 16, 'S4': 16,
        'C1': 17,   # Chlorine (Cl in notation, but uses C1 in Verloop)
        'B1': 35,   # Bromine (Br in notation, but uses B1 in Verloop)
        'I': 53,
    }

    elements = []
    i = 0
    while i < len(formula):
        ch = formula[i]
        if ch == 'H':
            elements.append(1)
            i += 1
        elif ch in ('C', 'N', 'O', 'F', 'P', 'S', 'I'):
            # Try multi-character symbols
            sym = ch
            j = i + 1
            while j < len(formula) and formula[j].isalnum():
                sym += formula[j]
                j += 1
            # Special: C1 = Cl, B1 = Br
            if sym == 'C1':
                elements.append(17)
            elif sym == 'B1':
                elements.append(35)
            elif sym in symbol_map:
                elements.append(symbol_map[sym])
            else:
                elements.append(6)  # Default to C
            i = j
        elif ch == 'X':
            # Ring closure marker, not an atom
            i += 2  # Skip X and following digit
        elif ch in '(),=*&R':
            i += 1
        elif ch == 'Z':
            # Radical reference Zxx
            i += 3
        else:
            i += 1

    return elements


def run_fortran(inp_path, ipr=1):
    """Run Fortran Sterimol and return parsed output."""
    with open(inp_path, 'r') as f:
        content = f.read()

    # Modify IPR if needed
    if ipr != 1:
        lines = content.split('\n')
        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped and stripped[0].isdigit() and len(stripped) <= 2:
                lines[idx] = str(ipr)
                break
        content = '\n'.join(lines)

    proc = subprocess.run(
        [STERIMOL_EXE],
        input=content.encode('ascii'),
        capture_output=True,
    )

    return proc.stdout.decode('ascii', errors='replace')


def parse_fortran_output(output):
    """Parse Fortran Sterimol output to extract coordinates and parameters."""
    coords = []
    radii = []

    # Parse coordinate block - find all lines matching the coordinate pattern
    coord_pattern = re.compile(r'^\s*(\d{1,3})\s+(\d{3})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{2})\s+(-?\d+\.\d{2})\s*$')
    for line in output.split('\n'):
        m = coord_pattern.match(line)
        if m:
            atom_num = int(m.group(1))
            radius = int(m.group(2)) / 100.0
            x = float(m.group(3))
            y = float(m.group(4))
            z = float(m.group(5))
            coords.append((x, y, z))
            radii.append(radius)

    # Parse steric parameters
    params = {}
    param_match = re.search(r'L=\s*(\d+\.\d{2})', output)
    if param_match:
        params['L'] = float(param_match.group(1))

    b_match = re.search(
        r'B\(1\) - B\(4\):\s*(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+B\(5\):\s*(\d+\.\d{2})',
        output
    )
    if b_match:
        params['B1'] = float(b_match.group(1))
        params['B2'] = float(b_match.group(2))
        params['B3'] = float(b_match.group(3))
        params['B4'] = float(b_match.group(4))
        params['B5'] = float(b_match.group(5))

    return coords, radii, params


def get_element_from_radius(radius, idx, elements_guess):
    """Infer atomic number from vdW radius.

    The Fortran program uses limited vdW radii that map multiple elements
    to the same radius. We use the .inp file element count for disambiguation.
    """
    if idx < len(elements_guess):
        return elements_guess[idx]

    # Fallback: guess from radius
    if abs(radius - 1.00) < 0.05:
        return 1  # H
    elif abs(radius - 1.50) < 0.05:
        return 6  # C (most common)
    elif abs(radius - 1.60) < 0.05:
        return 6  # C (aromatic)
    elif abs(radius - 1.70) < 0.05:
        return 6  # C (special)
    elif abs(radius - 1.40) < 0.05:
        return 8  # O
    elif abs(radius - 1.35) < 0.05:
        return 17  # Cl
    elif abs(radius - 1.45) < 0.05:
        return 35  # Br
    elif abs(radius - 1.95) < 0.05:
        return 53  # I
    return 6  # Default to C


def run_morfeus_sterimol(coords, radii, elements, n_rot_vectors=3600):
    """Run the morfeus-style Sterimol calculation in pure Python.

    Replicates the algorithm from sterimol_morfeus.py without the morfeus
    dependency.
    """
    coords = np.array(coords)
    radii = np.array(radii)

    # The first atom (H at origin) is the attachment reference
    # The second atom (first heavy atom) is the "attached" atom
    # We need a dummy atom placed along the bond axis

    # Attached atom is atom 1 (index 1, the first heavy atom after H)
    attached_idx = 1

    # Create dummy atom: extend from attached atom through atom 0
    # Vector from atom 0 (H at origin) to atom 1
    v = coords[attached_idx] - coords[0]
    bond_length = np.linalg.norm(v)
    v_unit = v / bond_length

    # Dummy atom placed at -v_unit * 0.40 (matching the +0.40 correction)
    dummy_pos = coords[attached_idx] - v_unit * 0.40

    # Add dummy atom to the beginning
    all_coords = np.vstack([dummy_pos, coords])
    all_radii = np.concatenate([[0.0], radii])

    # Dummy atom is index 0, attached atom is index 1
    dummy_index = 0

    # Translate so attached atom is at origin
    origin = all_coords[attached_idx]
    all_coords -= origin

    # Get bond axis (from dummy to attached)
    vector_2_to_1 = all_coords[attached_idx] - all_coords[dummy_index]
    bond_length = np.linalg.norm(vector_2_to_1)
    unit_axis = vector_2_to_1 / bond_length

    # Rotate so bond axis aligns with x-axis
    # Simple rotation: find rotation that maps unit_axis to [1,0,0]
    x_axis = np.array([1.0, 0.0, 0.0])
    all_coords = rotate_to_xaxis(all_coords, unit_axis)

    # Use atoms 2..end (skip dummy at 0 and attached at 1)
    active_coords = all_coords[2:]
    active_radii = all_radii[2:]

    # --- Calculate L ---
    # Fortran: AMI = min(x - R) for atoms 2..N, AL = -AMI + 0.40
    # This measures extent away from attachment point
    min_extent = np.min(active_coords[:, 0] - active_radii)
    L_value = -min_extent + 0.40

    # --- Calculate B1-B5 (Fortran-style algorithm) ---
    # Fortran rotates the molecule 53 steps around the x-axis.
    # At each step, it measures extents in 2 perpendicular directions:
    #   Dir1: rotated Y = Y*cos(θ) + Z*sin(θ)
    #   Dir2: rotated Z = Z*cos(θ) - Y*sin(θ)
    # For each direction: B(min) = min(CZ-RW), B(max) = max(CZ+RW)
    # ABS(B(1..4)) gives the 4 absolute extents.
    # The "optimal" step is where min(|B(1..4)|) is smallest.
    # B1-B4 = ABS(B(1..4)) at the optimal step; B5 = max(|B(1..4)|).

    theta = np.linspace(0, 2 * np.pi, n_rot_vectors)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    y = active_coords[:, 1]
    z = active_coords[:, 2]
    r = active_radii

    # Rotated coordinates for each angle
    # Dir1: rotated Y = y*cos(θ) + z*sin(θ)
    # Dir2: rotated Z = z*cos(θ) - y*sin(θ)
    rot_y = y[np.newaxis, :] * cos_t[:, np.newaxis] + z[np.newaxis, :] * sin_t[:, np.newaxis]
    rot_z = z[np.newaxis, :] * cos_t[:, np.newaxis] - y[np.newaxis, :] * sin_t[:, np.newaxis]

    # Extents in dir1: B(1) = min(rot_y - r), B(2) = max(rot_y + r)
    b1_per_angle = np.min(rot_y - r, axis=1)  # (n_angles,)
    b2_per_angle = np.max(rot_y + r, axis=1)

    # Extents in dir2: B(3) = min(rot_z - r), B(4) = max(rot_z + r)
    b3_per_angle = np.min(rot_z - r, axis=1)
    b4_per_angle = np.max(rot_z + r, axis=1)

    # Absolute extents
    ab1 = np.abs(b1_per_angle)
    ab2 = np.abs(b2_per_angle)
    ab3 = np.abs(b3_per_angle)
    ab4 = np.abs(b4_per_angle)

    # Minimum of the 4 at each angle
    min_abs = np.minimum(np.minimum(ab1, ab2), np.minimum(ab3, ab4))

    # Optimal angle: where min_abs is smallest
    optimal_idx = np.argmin(min_abs)

    # B1-B4 at optimal angle
    B1_value = float(ab1[optimal_idx])
    B2_value = float(ab2[optimal_idx])
    B3_value = float(ab3[optimal_idx])
    B4_value = float(ab4[optimal_idx])
    B5_value = float(max(B1_value, B2_value, B3_value, B4_value))

    return {
        'L': round(L_value, 2),
        'B1': round(B1_value, 2),
        'B2': round(B2_value, 2),
        'B3': round(B3_value, 2),
        'B4': round(B4_value, 2),
        'B5': round(B5_value, 2),
    }


def rotate_to_xaxis(coords, vector):
    """Rotate coordinates so that 'vector' aligns with the x-axis."""
    vector = vector / np.linalg.norm(vector)

    if np.allclose(vector, [1, 0, 0]):
        return coords
    if np.allclose(vector, [-1, 0, 0]):
        return coords @ np.diag([-1, -1, 1])

    # Rotation axis is cross product of vector and x-axis
    x_axis = np.array([1.0, 0.0, 0.0])
    rotation_axis = np.cross(vector, x_axis)
    sin_angle = np.linalg.norm(rotation_axis)
    cos_angle = np.dot(vector, x_axis)

    if sin_angle < 1e-10:
        return coords

    rotation_axis = rotation_axis / sin_angle

    # Rodrigues' rotation formula
    K = np.array([
        [0, -rotation_axis[2], rotation_axis[1]],
        [rotation_axis[2], 0, -rotation_axis[0]],
        [-rotation_axis[1], rotation_axis[0], 0]
    ])
    R = np.eye(3) + sin_angle * K + (1 - cos_angle) * K @ K

    return (R @ coords.T).T


def compare_molecule(name, inp_path):
    """Compare Fortran vs Python for a single molecule."""
    # Parse elements from .inp
    elements = parse_inp_file(inp_path)

    # Run Fortran
    output = run_fortran(inp_path)
    coords, radii, fortran_params = parse_fortran_output(output)

    if not coords:
        return {'error': 'Failed to parse coordinates'}

    # Map elements to coordinates
    # The Fortran output starts with atom 1 (H at origin), then heavy atoms
    # Our element list from the formula should match
    n_atoms = len(coords)

    # Build element list matching Fortran's atom ordering
    # First atom is H, rest come from formula parsing
    final_elements = [1] + elements[1:] if elements else [1]

    # Adjust if needed
    while len(final_elements) < n_atoms:
        final_elements.append(6)  # Default to C

    # Run morfeus-style Python calculation
    try:
        morfeus_params = run_morfeus_sterimol(coords, radii, final_elements)
    except Exception as e:
        morfeus_params = {'error': str(e)}

    return {
        'fortran': fortran_params,
        'morfeus_python': morfeus_params,
    }


def main():
    """Run comparison for all molecules in test/data."""
    inp_files = sorted([
        f for f in os.listdir(DATA_DIR) if f.endswith('.inp')
    ])

    print(f"{'Molecule':<10s} {'F-L':>6s} {'F-B1':>6s} {'F-B2':>6s} {'F-B3':>6s} {'F-B4':>6s} {'F-B5':>6s} "
          f"{'P-L':>6s} {'P-B1':>6s} {'P-B2':>6s} {'P-B3':>6s} {'P-B4':>6s} {'P-B5':>6s}")
    print("-" * 110)

    max_d = {}
    n_ok = 0
    n_total = 0
    params = ['L', 'B1', 'B2', 'B3', 'B4', 'B5']

    for inp_file in inp_files:
        name = inp_file[:-4]  # Remove .inp
        inp_path = os.path.join(DATA_DIR, inp_file)

        result = compare_molecule(name, inp_path)

        if 'error' in result:
            print(f"{name:<10s} ERROR: {result['error']}")
            continue

        f = result['fortran']
        m = result['morfeus_python']

        if 'error' in m:
            print(f"{name:<10s} Python error: {m['error']}")
            continue

        diffs = {}
        all_ok = True
        for p in params:
            fd = f.get(p, 0)
            md = m.get(p, 0)
            d = abs(fd - md)
            diffs[p] = d
            if p in max_d:
                max_d[p] = max(max_d[p], d)
            else:
                max_d[p] = d
            if d >= 0.06:
                all_ok = False

        n_total += 1
        if all_ok:
            n_ok += 1

        f_str = ' '.join(f"{f.get(p, 0):>6.2f}" for p in params)
        m_str = ' '.join(f"{m.get(p, 0):>6.2f}" for p in params)
        flag = ' *' if not all_ok else '   '
        print(f"{name:<10s} {f_str} {m_str}{flag}")

    print("-" * 110)
    for p in params:
        print(f"  Max |d{p}| = {max_d.get(p, 0):.2f}")
    print(f"\nMolecules within 0.06 A tolerance (all params): {n_ok}/{n_total}")
    print()
    print("Notes:")
    print("  F = Fortran (original STERIMOL, 53 rotation steps)")
    print("  P = Python (replicated Fortran algorithm, 3600 rotation steps)")
    print("  The Python implementation replicates the Fortran algorithm with")
    print("  finer angular sampling. Small differences (<=0.06 A) may arise")
    print("  from discretization and floating-point precision.")
    print("  '*' indicates any parameter differs by more than 0.06 A")


if __name__ == '__main__':
    main()
