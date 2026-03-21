import numpy as np
from fractions import Fraction
from math import gcd
from functools import reduce

# ============================================================
# Inverse design of heterodeformation from Burgers/line pairs
# ============================================================
# Edit ONLY the example block below. Additional examples are listed
# in the accompanying README. If not all line vectors have integer
# lattice coordinates, the script will optionally use pyoilab to
# compute the CSL primitive vectors.


# ------------------------------------------------------------
# Exact rational arithmetic helpers
# ------------------------------------------------------------
def as_fraction(x):
    """Convert x to Fraction unless it is already a Fraction."""
    return x if isinstance(x, Fraction) else Fraction(x)


def gcd_fraction(x, y):
    """Greatest common divisor of two rational numbers.

    If x = a/b and y = c/d in lowest terms, define
        gcd(x, y) = gcd(ad, cb)/(bd).
    """
    x = as_fraction(x)
    y = as_fraction(y)
    a, b = x.numerator, x.denominator
    c, d = y.numerator, y.denominator
    return Fraction(gcd(a * d, c * b), b * d)


def gcd_fraction_list(vals):
    """Greatest common divisor of a list of rational numbers."""
    vals = [as_fraction(v) for v in vals]
    nz = [v for v in vals if v != 0]
    if len(nz) == 0:
        return Fraction(1, 1)
    return reduce(gcd_fraction, nz)


def dot_frac(u, v):
    """Exact dot product of two rational vectors."""
    return sum(as_fraction(ui) * as_fraction(vi) for ui, vi in zip(u, v))


def outer_frac(u, v):
    """Exact outer product of two rational 2-vectors."""
    return np.array(
        [[as_fraction(u[i]) * as_fraction(v[j]) for j in range(2)] for i in range(2)],
        dtype=object,
    )


def frac_abs(x):
    """Absolute value for Fraction."""
    x = as_fraction(x)
    return x if x >= 0 else -x


def primitive_orthogonal(l):
    """Smallest reciprocal vector orthogonal to l.

    If l = (l1, l2), then
        mbar = (l2, -l1) / gcd(l1, l2).
    """
    l1 = as_fraction(l[0])
    l2 = as_fraction(l[1])
    g = gcd_fraction_list([l1, l2])
    return np.array([l2 / g, -l1 / g], dtype=object)


def frac_matrix_to_float(M):
    """Convert a matrix with Fraction entries to a NumPy float array."""
    return np.array(
        [[float(M[i, j]) for j in range(M.shape[1])] for i in range(M.shape[0])],
        dtype=float,
    )


def inv2_frac(M):
    """Exact inverse of a 2x2 matrix with Fraction entries."""
    a, b = as_fraction(M[0, 0]), as_fraction(M[0, 1])
    c, d = as_fraction(M[1, 0]), as_fraction(M[1, 1])
    det = a * d - b * c
    if det == 0:
        raise ValueError("Q is singular, so Q^{-1} does not exist.")
    return np.array([[d / det, -b / det], [-c / det, a / det]], dtype=object)


def all_integer_lattice_coords(vecs):
    """Check whether every component of every vector is an integer."""
    for vec in vecs:
        for x in vec:
            if as_fraction(x).denominator != 1:
                return False
    return True


def print_fraction_matrix(name, M):
    """Pretty-print a 2x2 matrix with Fraction entries."""
    print(f"\n{name} =")
    print(np.array(M, dtype=object))


# ------------------------------------------------------------
# Example block
# ------------------------------------------------------------
# DEFAULT EXAMPLE:
# 2D simple twist in bilayer graphene (Fig. 7 in the manuscript)
# Replace A, b1u, b2u, b3u, l1u, l2u with any case listed in README.

# Structure matrix of the 2D triangular Bravais lattice
# Graphene: a = 2.46 A -> A = (a/2) [[0, -sqrt(3)], [2, -1]]
# MoS2:     a = 3.19702 A -> same normalized matrix scaled by a
A = (2.46 / 2.0) * np.array([[0.0, -np.sqrt(3.0)],
                             [2.0, -1.0]], dtype=float)

# Burgers vectors in lattice coordinates
b1u = np.array([Fraction(-1, 3), Fraction(-2, 3)], dtype=object)
b2u = np.array([Fraction(-1, 3), Fraction(1, 3)], dtype=object)
b3u = np.array([Fraction(2, 3), Fraction(1, 3)], dtype=object)

# Line vectors in lattice coordinates
l1u = np.array([Fraction(-110, 1), Fraction(-221, 1)], dtype=object)
l2u = np.array([Fraction(-111, 1), Fraction(110, 1)], dtype=object)

# beta = 1 for 1D networks and 2D triangular networks,
# beta = 3 for 2D honeycomb networks.
# Set this explicitly after choosing an example from the README.
beta = Fraction(1, 1)


# ------------------------------------------------------------
# Main calculation
# ------------------------------------------------------------
l3u = -(l1u + l2u)
B = [b1u, b2u, b3u]
L = [l1u, l2u, l3u]

print("Cartesian Burgers vectors:")
for i, b in enumerate(B, start=1):
    print(f"  b{i} = {A @ frac_matrix_to_float(b.reshape(2,1)).ravel()}")

print("\nCartesian line vectors:")
for i, l in enumerate(L, start=1):
    print(f"  l{i} = {A @ frac_matrix_to_float(l.reshape(2,1)).ravel()}")

mbar = [primitive_orthogonal(li) for li in L]
z1 = frac_abs(dot_frac(L[1], mbar[0]))
z2 = frac_abs(dot_frac(L[2], mbar[1]))
z3 = frac_abs(dot_frac(L[0], mbar[2]))
Z = [z1, z2, z3]

if any(z == 0 for z in Z):
    raise ValueError(
        "Q is undefined because some z^i = 0. "
        "This target network is degenerate (the line vectors are collinear)."
    )

Q = np.array([[Fraction(1, 1), Fraction(0, 1)],
              [Fraction(0, 1), Fraction(1, 1)]], dtype=object)
for i in range(3):
    Q -= outer_frac(B[i], mbar[i]) / (beta * Z[i])

Qinv = inv2_frac(Q)
Qf = frac_matrix_to_float(Q)
Qinv_f = frac_matrix_to_float(Qinv)

Ainv = np.linalg.solve(A, np.eye(2))
Finv = A @ Qf @ Ainv
F = A @ Qinv_f @ Ainv

print_fraction_matrix("Q", Q)
print_fraction_matrix("Q^{-1}", Qinv)
print("\nF^{-1} =")
print(Finv)
print("\nF =")
print(F)

# ------------------------------------------------------------
# Primitive simulation-cell vectors
# ------------------------------------------------------------
import pyoilab as gb
latticeA = gb.Lattice2D(A)
latticeB = gb.Lattice2D(A, F)
bicrystal = gb.BiCrystal2D(latticeA, latticeB, True)

print("\nCSL primitive vectors/box vectors from oILAB:")
ell1= gb.LatticeVector2D(np.array([1,0],dtype=np.int64),bicrystal.csl)
ell2= gb.LatticeVector2D(np.array([0,1],dtype=np.int64),bicrystal.csl)
print("Box vector 1 = ", bicrystal.csl.latticeBasis[:, 0])
print("Box vector 2 = ", bicrystal.csl.latticeBasis[:, 1])

# Output lattice A
print("Expected number of atoms in lattice A = ", bicrystal.sigmaA)
latticeA_boxVectors = [bicrystal.getLatticeVectorInA(ell1), bicrystal.getLatticeVectorInA(ell2)]
boxObjA=latticeA.box(latticeA_boxVectors,"latticeA")
assert len(boxObjA)==bicrystal.sigmaA, "Number of atoms is not equal to the expected number of atoms in lattice A"

# Output lattice B
print("Expected Number of atoms in lattice B = ", bicrystal.sigmaB)
latticeB_boxVectors = [bicrystal.getLatticeVectorInB(ell1), bicrystal.getLatticeVectorInB(ell2)]
boxObjB=latticeB.box(latticeB_boxVectors,"latticeB")
assert len(boxObjB)==bicrystal.sigmaB, "Number of atoms is not equal to the expected number of atoms in lattice A"

