# Inverse design of heterodeformations for strain soliton networks

This directory contains example scripts and data used to construct strain soliton networks from prescribed Burgers vector–line vector pairs, as described in the manuscript.

The workflow implements an **inverse design framework**, where a target soliton network is specified and the corresponding heterodeformation and atomic configuration are constructed.

---

## Repository structure

The implementation is organized as follows:

- `inverse_design_core.py`  
  Computes the deformation gradient \(F\) from Burgers vectors and line directions using exact arithmetic.

- `layers.py`  
  Defines 2D multilattices (Bravais lattice + basis atoms).

- `bilayers.py`  
  Defines bilayer structures and applies heterodeformations.

- `lammps_writer.py`  
  Writes LAMMPS data files from generated atomic configurations.

- `example_inverse_design.py`  
  A minimal working example demonstrating the full workflow.

- `relax_graphene.in`  
  LAMMPS input script for relaxation of bilayer graphene.

- `relax_mos2.in`  
  LAMMPS input script for relaxation of bilayer MoS$_2$.

---

## Running the example

First, build the oILAB library and Python bindings from the repository root:

```bash
cmake -B build -S .
cmake --build build
```

Ensure the Python bindings are available:

```bash
export PYTHONPATH=build/bindings
```

Then run the example script:

```bash
python example_inverse_design.py
```

The script will:
1. Compute the deformation gradient \(F\) from prescribed network data
2. Construct a heterodeformed bilayer
3. Compute the CSL using oILAB
4. Generate atomic coordinates
5. Write a LAMMPS data file (`initial.data`)

---

## LAMMPS relaxation scripts

This repository includes LAMMPS input scripts used to relax the generated bilayer structures:

- `relax_graphene.in` — relaxation of bilayer graphene using REBO + Kolmogorov–Crespi potentials  
- `relax_mos2.in` — relaxation of bilayer MoS$_2$ using SW + interlayer potential (ILP)

These scripts take the generated `initial.data` file as input and perform energy minimization to obtain relaxed configurations.

---

### Usage

After generating the configuration (`initial.data`) using the Python script, run LAMMPS:

```bash
mpirun -np 4 lmp_mpi -in relax_graphene.in
```

---

## Modifying the example

The script `example_inverse_design.py` contains **a single example**.

To reproduce other cases from the manuscript, edit the **example block**:

- Choose the bilayer:
  ```python
  bilayer = make_ab_graphene_bilayer()
  # or
  bilayer = make_aa_prime_mos2_bilayer()
  ```

- Specify:
  - `b1u`, `b2u`, `b3u` (Burgers vectors)
  - `l1u`, `l2u` (line directions)
  - `beta`

The third line direction is automatically set as:
```python
l3u = -(l1u + l2u)
```

Use:
- `beta = 1` for 1D and triangular networks
- `beta = 3` for honeycomb networks

---

## CSL construction

The coincidence site lattice (CSL) is computed using the oILAB Python bindings.

The script constructs a `BiCrystal2D` object and obtains the CSL primitive vectors:

```python
bicrystal = hetero_bilayer.bicrystal
```

The simulation cell is defined using CSL lattice vectors:

```python
ell1 = gb.LatticeVector2D([1, 0], bicrystal.csl)
ell2 = gb.LatticeVector2D([0, 1], bicrystal.csl)
```

The user may replace these with any integer combination of CSL basis vectors.

---

## Bilayer conventions

The atomic configurations are constructed to be consistent with the LAMMPS scripts used in this work.

### Bilayer graphene
- Top layer: atom type `1`
- Bottom layer: atom type `2`

### Bilayer MoS$_2$
- Top layer: atom types `1, 2, 3`
- Bottom layer: atom types `4, 5, 6`

These conventions are required for compatibility with the provided LAMMPS input scripts.

---

## Examples from the manuscript

Below are the Burgers vector–line vector pairs used in the paper.

---

## Bilayer graphene

Structure matrix:

```python
A = (2.46 / 2.0) * np.array([[0.0, -np.sqrt(3.0)],
                             [2.0, -1.0]], dtype=float)
```

### 1D network (Fig. 6b)

```python
b1u = np.array([Fraction(-1, 1), Fraction(1, 1)], dtype=object)
b2u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)

l1u = np.array([Fraction(-1, 1), Fraction(1, 1)], dtype=object)
l2u = np.array([Fraction(569, 1), Fraction(567, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D simple twist (Figs. 2a, 7)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(-110, 1), Fraction(-221, 1)], dtype=object)
l2u = np.array([Fraction(-111, 1), Fraction(110, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D complex twist (Fig. 9a)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(-151, 3), Fraction(-305, 3)], dtype=object)
l2u = np.array([Fraction(-154, 3), Fraction(151, 3)], dtype=object)

beta = Fraction(1, 1)
```

### 2D moir\'e degeneracy 1 (Fig. 8a)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(1, 1), Fraction(-280, 1)], dtype=object)
l2u = np.array([Fraction(-281, 1), Fraction(-1, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D moir\'e degeneracy 2 (Fig. 8b)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(1, 1), Fraction(-280, 1)], dtype=object)
l2u = np.array([Fraction(-280, 1), Fraction(-281, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D moir\'e degeneracy 3 (Fig. 8c)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(1, 1), Fraction(-280, 1)], dtype=object)
l2u = np.array([Fraction(-279, 1), Fraction(-561, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D complex strain (Fig. 9b)

```python
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

l1u = np.array([Fraction(1, 1), Fraction(-72, 1)], dtype=object)
l2u = np.array([Fraction(-218, 3), Fraction(-25, 1)], dtype=object)

beta = Fraction(1, 1)
```

---

## Bilayer MoS$_2$

Structure matrix:

```python
A = (3.19702 / 2.0) * np.array([[0.0, -np.sqrt(3.0)],
                                [2.0, -1.0]], dtype=float)
```

### 1D network (Fig. 6a)

```python
b1u = np.array([Fraction(0, 1), Fraction(1, 1)], dtype=object)
b2u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)

l1u = np.array([Fraction(-1, 1), Fraction(1, 1)], dtype=object)
l2u = np.array([Fraction(121, 1), Fraction(119, 1)], dtype=object)

beta = Fraction(1, 1)
```

### 2D twist (Fig. 2b)

```python
b1u = np.array([Fraction(-1, 1), Fraction(0, 1)], dtype=object)
b2u = np.array([Fraction(1, 1), Fraction(1, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(-1, 1)], dtype=object)

l1u = np.array([Fraction(-332, 3), Fraction(-1, 3)], dtype=object)
l2u = np.array([Fraction(331, 3), Fraction(332, 3)], dtype=object)

beta = Fraction(3, 1)
```

### 2D moir\'e degeneracy 1 (Fig. 10a)

```python
b1u = np.array([Fraction(-1, 1), Fraction(0, 1)], dtype=object)
b2u = np.array([Fraction(1, 1), Fraction(1, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(-1, 1)], dtype=object)

l1u = np.array([Fraction(-64, 1), Fraction(1, 1)], dtype=object)
l2u = np.array([Fraction(65, 1), Fraction(64, 1)], dtype=object)

beta = Fraction(3, 1)
```

### 2D moir\'e degeneracy 2 (Fig. 10b)

```python
b1u = np.array([Fraction(-1, 1), Fraction(0, 1)], dtype=object)
b2u = np.array([Fraction(1, 1), Fraction(1, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(-1, 1)], dtype=object)

l1u = np.array([Fraction(-108, 1), Fraction(-85, 1)], dtype=object)
l2u = np.array([Fraction(87, 1), Fraction(107, 1)], dtype=object)

beta = Fraction(3, 1)
```

### 2D complex strain (Fig. 11)

```python
b1u = np.array([Fraction(-1, 1), Fraction(0, 1)], dtype=object)
b2u = np.array([Fraction(1, 1), Fraction(1, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(-1, 1)], dtype=object)

l1u = np.array([Fraction(-216, 5), Fraction(1, 5)], dtype=object)
l2u = np.array([Fraction(217, 5), Fraction(216, 5)], dtype=object)

beta = Fraction(3, 1)
```

---

## Notes

1. The script computes \(l^3 = -(l^1 + l^2)\) automatically.
2. The deformation gradient is computed using exact arithmetic before conversion to floating point.
3. CSL vectors are obtained through oILAB and define the simulation cell.
4. Atomic configurations are written in LAMMPS-compatible format.
