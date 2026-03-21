This directory contains example scripts used to construct strain soliton networks from prescribed Burgers-vector / line-vector pairs, as described in the manuscript.

## Running the examples

First, configure and build oILAB from the repository root:

```bash
cmake -B build -S .
cmake --build build
```

Next, navigate to the heterostrain example directory in the build tree:

```bash
cd build/examples/python/heterostrain
```

Run the script using:

```bash
PYTHONPATH=../../../bindings python heterostrain.py
```

The compiled Python extension module `pyoilab` is located in:

```
build/bindings/
```

and is made accessible via `PYTHONPATH`.

All output files are written to the current working directory:

```
build/examples/python/heterostrain/
```

ensuring that the source tree remains unchanged.

---

# Inverse design examples for the manuscript

This README lists the Burgers-vector/line-vector pairs used in the paper and explains how to run the accompanying script.

The script `heterodeformation.py` contains **one example only**. To reproduce another case from the paper, edit the following quantities in the **example block**:

- `A` (choose graphene or MoS$_2$ lattice constant)
- `b1u`, `b2u`, `b3u`
- `l1u`, `l2u` (the script sets `l3u = -(l1u + l2u)`)
- `beta`

Use:
- `beta = 1` for 1D networks and 2D triangular networks
- `beta = 3` for 2D honeycomb networks

If **all** line vectors have integer lattice coordinates, the script computes the simulation cell vectors directly as
`ell^1 = A l^1` and `ell^2 = A l^2`.
If **any** line vector has non-integer lattice coordinates, the script uses the **oILAB Python bindings** (`pyoilab`) to compute the CSL primitive vectors.

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

1. The script automatically sets `l3u = -(l1u + l2u)`.
2. If all line vectors are integral in lattice coordinates, the simulation cell vectors are the line vectors themselves.
3. If at least one line vector has non-integer lattice coordinates, install the oILAB Python bindings and rerun the script to obtain the CSL primitive vectors.
