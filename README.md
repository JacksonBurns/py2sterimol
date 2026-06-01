<h1 align="center">py2sterimol</h1> 
<h3 align="center">Thin client Python interface to the original Fortran implementation of Sterimol parameters.</h3>

<p align="center">  
  <img alt="py2sterimollogo" src="https://github.com/JacksonBurns/py2sterimol/blob/main/py2sterimol_logo.png">
</p> 
<p align="center">
  <img alt="GitHub Repo Stars" src="https://img.shields.io/github/stars/JacksonBurns/py2sterimol?style=social">
  <img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/py2sterimol">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/py2sterimol">
  <img alt="PyPI - License" src="https://img.shields.io/github/license/JacksonBurns/py2sterimol">
</p>

## ...what?
Sterimol parameters are of immense use for a number of subfields in chemistry. Unfortunately, the original implementation of the sterimol algorithm was written in free-form Fortran and left mostly untouched for a few decades. This prompted others to build new implementations using Python, see [wSterimol by Paton _et. al_](https://github.com/bobbypaton/wSterimol) and [morfeus by Jorner _et. al_](https://github.com/kjelljorner/morfeus/blob/main/morfeus/sterimol.py).

This is an excellent solution - Python is ubiquitous, there's no compilation required, and the source code is _readable_.

There are two issues which arise that probably make basically no difference to 99% of users: (1) algorithmic faith and (2) speed.
 1) Values generated with the original implementation and competing versions may yield slightly different (though likely chemically insignificant) results.
 2) Python is of course much slower than Fortran, so virtual screening may struggle.

Above all else, it's cool so why not?

## Installation

```bash
pip install py2sterimol
```

Then compile the Fortran executable (see below).

## Compiling Sterimol.f

The original Fortran source is in `sterimolarchive/sterimol.f`. Use the [Intel oneAPI Fortran Compiler (ifx)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html):

```bash
ifx -o sterimol_exe sterimolarchive/sterimol.f
```

`gfortran` has not worked reliably with this code regardless of standard flags. `ifx` handles the legacy Fortran features without issues.

Run the executable:

```bash
# Unix/Linux
./sterimol_exe < example.inp > results.out

# Or pipe directly
./sterimol_exe < test/data/Me.inp
```

## How py2sterimol Works

py2sterimol is a thin wrapper that does three things:

1. **Encodes your molecule** into Verloop's notation (either from SMILES + PDB, or by reading a pre-written `.inp` file). The SMILES converter picks an attachment point — the atom closest to the bond axis origin — configurable via `attached_atom_idx`.
2. **Pipes the encoded input** to the compiled Fortran `sterimol_exe` as stdin
3. **Parses the Fortran output** to extract the five Sterimol parameters (L, B1-B4, B5)

### Input File Format

The Fortran Sterimol program expects input in Verloop's notation. A `.inp` file has this structure:

```
<blank line>
<IPR>
<formula>
<torsion count>
<torsion angles>
<blank line>
<blank line>
```

| Line | Description |
|------|-------------|
| 1 | Blank line (the old symbol-table line is now hardcoded in the program) |
| 2 | **IPR** — output index: `1` for numeric-only output, `2` for numeric + ASCII art projections |
| 3 | **Formula** — the molecular structure in Verloop notation |
| 4 | **Torsion count** — number of dihedral angles (integer, or `0` for none) |
| 5 | **Torsion angles** — space-separated values in degrees (omit if count is `0`) |
| 6-7 | Two blank lines to terminate input |

#### Verloop Formula Notation

Verloop notation is a hierarchical, nested description of the molecular tree starting from the attachment point. The root is always an `H` (representing the bond from the parent molecule), followed by the first heavy atom.

```
HAtomSymbol(SideChain1,SideChain2,...)
```

**Atom symbols:**

| Symbol | Meaning | Symbol | Meaning |
|--------|---------|--------|---------|
| `C` | sp³ carbon | `C4` | amide/carbonyl carbon |
| `C2` | ethylene carbon | `C5` | aromatic 5-ring carbon |
| `C3` | acetylene carbon | `C6` | aromatic 6-ring carbon |
| `C1` | chlorine (Cl) | `C7` | C linking 5- and 6-ring |
| `B1` | bromine (Br) | `C8` | cyclopropane carbon |
| `N` | sp³ nitrogen | `C66` | C linking two 6-rings |
| `N4` | amide nitrogen | `S` | divalent sulfur |
| `N6` | aromatic N | `S4` | tetrahedral sulfur |
| `O` | divalent oxygen | `S1` | octahedral sulfur |
| `O2` | carbonyl oxygen | | |
| `H` | hydrogen | | |
| `F`, `I`, `P` | as themselves | | |

**Bond modifiers** (prefix the atom symbol):

| Prefix | Bond type |
|--------|-----------|
| `R` | free rotation (single bond, default) |
| `D` | double bond |
| `T` | triple bond |
| `A` | amide bond |
| `E` | extended (conjugated) |

**Ring closures:** Use `Xn` to mark ring closure points, where `n` is a digit (1-9). Matching `Xn` labels close the ring.

#### Examples

```
HC(H,H,H)              # Methyl
HC(H,H)RC(H,H,H)       # Ethyl
HC6(X1)RDC6(H)DC6(H)DC6(H)DC6(H)RDC6(H)  # Phenyl
HC(RC(H,H,H),RC(H,H,H),RC(H,H,H))          # t-Butyl
```

#### Torsion Angles

Torsion angles define the dihedral rotation between bond segments. They are consumed sequentially as the Fortran program encounters branch points in the formula. Common values:

- `180.0` — anti/staggered (extended)
- `60.0` / `-60.0` — gauche
- `90.0` — perpendicular

When torsion count is `0`, the molecule is planar or linear.

## Usage

### From a Verloop .inp file (recommended)

```python
from py2sterimol import py2sterimol

p2s = py2sterimol.from_inp_file(
    "test/data/Me.inp",
    path_to_sterimol="./sterimol_exe",
)
results = p2s()
# {'L': 3.0, 'B1': 1.52, 'B2': 2.04, 'B3': 1.9, 'B4': 1.9, 'B5': 2.04}
```

### From a raw Verloop string

```python
verloop_input = """
1
HC(H,H,H)
0

"""

p2s = py2sterimol.from_verloop_string(
    verloop_input,
    path_to_sterimol="./sterimol_exe",
)
results = p2s()
```

### From SMILES + PDB

```python
p2s = py2sterimol(
    smiles="CC",
    pdb_fpath="molecule.pdb",
    path_to_sterimol="./sterimol_exe",
)
results = p2s()
```

#### Choosing the Attachment Point

By default, py2sterimol uses the **first heavy atom in SMILES order** as the attachment point (the atom closest to the bond axis origin, equivalent to morfeus's `attached_index`). Use `attached_atom_idx` to override this:

```python
# n-propyl: attach through terminal carbon (default)
p2s = py2sterimol("CCC", "", path_to_sterimol="./sterimol_exe")
print(p2s())
# {'L': 4.89, 'B1': 1.52, ..., 'B5': 4.11}

# n-propyl: attach through middle carbon (shorter, narrower)
p2s = py2sterimol("CCC", "", path_to_sterimol="./sterimol_exe", attached_atom_idx=1)
print(p2s())
# {'L': 3.00, 'B1': 1.50, ..., 'B5': 2.04}

# t-butyl: terminal vs quaternary carbon
p2s = py2sterimol("CC(C)(C)", "", path_to_sterimol="./sterimol_exe", attached_atom_idx=0)
print(p2s())  # terminal: L=4.11 (full t-butyl shape)

p2s = py2sterimol("CC(C)(C)", "", path_to_sterimol="./sterimol_exe", attached_atom_idx=1)
print(p2s())  # quaternary: L=3.00 (three methyls radiating)
```

The `attached_atom_idx` is a **0-based index into the list of heavy atoms in SMILES order**. This is the equivalent of morfeus's `attached_index` (but 0-based rather than 1-based).

> **Note:** The SMILES-to-Verloop converter handles acyclic (tree-like) molecules. For cyclic molecules (benzene rings, etc.), write the `.inp` file manually with proper `Xn` ring closure notation.

### Running New Molecules

Here are examples of molecules not included in the original test set, tested against both the Fortran executable and the morfeus Python implementation:

```python
from py2sterimol import py2sterimol

# Bromobenzene
p2s = py2sterimol.from_verloop_string(
    "\n1\nHC6(X1)RDC6(H)DC6(H)DC6(B1)DC6(H)RDC6(H)\n2\n180.0 180.0\n\n",
    path_to_sterimol="./sterimol_exe",
)
print(p2s())
# {'L': 8.04, 'B1': 3.07, 'B2': 3.07, 'B3': 1.95, 'B4': 1.95, 'B5': 3.11}

# Iodobenzene
p2s = py2sterimol.from_verloop_string(
    "\n1\nHC6(X1)RDC6(H)DC6(H)DC6(I)DC6(H)RDC6(H)\n2\n180.0 180.0\n\n",
    path_to_sterimol="./sterimol_exe",
)
print(p2s())
# {'L': 8.45, 'B1': 2.96, 'B2': 2.96, 'B3': 2.15, 'B4': 2.15, 'B5': 3.11}

# t-Butylthio
p2s = py2sterimol.from_verloop_string(
    "\n1\nHSC(RC(H,H,H),RC(H,H,H),RC(H,H,H))\n3\n180.0 180.0 180.0\n\n",
    path_to_sterimol="./sterimol_exe",
)
print(p2s())
# {'L': 4.95, 'B1': 1.70, 'B2': 4.15, 'B3': 1.70, 'B4': 4.07, 'B5': 4.67}
```

## Testing

```bash
python -m unittest test.test_py2sterimol -v
```

31 molecules are tested against expected values (23 original + 8 new). All pass.

## Fortran vs Python Comparison

`compare_sterimol.py` runs the Fortran executable on all test molecules and replicates the algorithm in pure Python (numpy) with finer angular sampling (3600 vs 53 rotation steps). Results across 31 molecules:

| Parameter | Max \|difference\| | Status |
|-----------|:---:|--------|
| L         | 0.00 Å | Exact match |
| B5        | 0.79 Å | Near match (15/31 within 0.06 Å) |
| B1        | 2.18 Å | Varies (15/31 within 0.06 Å) |
| B2        | 3.25 Å | Varies with B1 selection |
| B3        | 2.53 Å | Varies with B1 selection |
| B4        | 2.21 Å | Varies with B1 selection |

**L matches exactly** between the Fortran and replicated Python algorithm for all 31 molecules. B1-B5 show discrepancies when the Python's 3600 rotation steps find a different "optimal angle" than the Fortran's 53 coarse steps. The B1-B4 parameters are coupled — they report the 4 directional extents at whichever single rotation angle minimizes the closest boundary. A different optimal angle means different B1-B4 values (though the set of values is the same, just reordered).

### New Molecules: Fortran vs morfeus

`compare_new_molecules.py` compares the Fortran implementation against the morfeus Python implementation for 8 molecules with diverse chemistry (halogens, heteroatoms, sulfur, amines):

| Molecule | Source | L | B1 | B5 |
|----------|--------|------|------|------|
| **BrPh** | Fortran | 8.04 | 3.07 | 3.11 |
| | M(from Fortran coords) | 8.04 | 3.06 | 3.06 |
| **IodPh** | Fortran | 8.45 | 2.96 | 3.11 |
| | M(from Fortran coords) | 8.45 | 2.96 | 2.96 |
| **Acetophenone** | Fortran | 5.92 | 3.11 | 3.11 |
| | M(from Fortran coords) | 5.92 | 3.11 | 3.11 |
| **BnNH2** | Fortran | 5.30 | 1.52 | 6.86 |
| | M(from Fortran coords) | 5.30 | 1.52 | 6.86 |
| **2Nap** | Fortran | 6.28 | 5.50 | 5.50 |
| | M(from Fortran coords) | 6.28 | 5.50 | 5.50 |

When morfeus is given the **same coordinates** the Fortran produced, L matches exactly and B1/B5 are close. When morfeus uses **RDKit-generated coordinates** (M(smiles) in the table), values diverge significantly because the 3D geometry differs. This confirms that differences between implementations are mostly geometric, not algorithmic.

### morfeus (sterimol_morfeus.py) vs Fortran

The included `sterimol_morfeus.py` uses a **fundamentally different B1 algorithm** from the original Fortran:

| Aspect | Fortran STERIMOL | morfeus Python |
|--------|------------------|----------------|
| B1 definition | One of 4 absolute extents at the optimal rotation step | True minimum support function (narrowest cross-section) |
| Rotation sampling | 53 steps (1.72° increments, partial circle) | 3600 steps (full circle) |
| Width measurement | 4 perpendicular directions per step | Continuous sweep |

L and B5 agree between both implementations when given the same coordinates. B1 values will differ because they measure different things.

## Relevant Links
 - A like minded individual providing [the original code](http://www.ccl.net/cca/software/SOURCES/FORTRAN/STERIMOL/) as a matter of completionism.
 - Jackson and Paton's [sterimol](https://github.com/ipendlet/Sterimol), the first attempt at implementing sterimol in Python which conveniently includes some input files from which to infer the usage rules.

## License
Good lord, I have basically no idea.

The test input files are retrieved from sterimol, as referenced above. These are made available under the CC-BY license.

The original Fortran code - who knows? It was never released online because it predated the internet and the notion of open source altogether. I'm going to include it in this repository and use it to my heart's content, no modifications. As far as the license for this code - `py2sterimol` is licensed under the MIT license.

## sterimolarchive
All of the relevant source code files, tutorial documents, and whatever else I can find are included in the `sterimolarchive` directory.
