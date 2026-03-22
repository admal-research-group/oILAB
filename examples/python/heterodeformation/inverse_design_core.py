from dataclasses import dataclass
import numpy as np
from fractions import Fraction


# ------------------------------------------------------------
# Utilities for Fraction-based linear algebra
# ------------------------------------------------------------
def frac_array(x):
    return np.array(x, dtype=object)


def to_fraction(x):
    if isinstance(x, Fraction):
        return x
    return Fraction(x).limit_denominator()


def frac_dot(a, b):
    return sum(to_fraction(ai) * to_fraction(bi) for ai, bi in zip(a, b))


def frac_outer(a, b):
    return np.array([[to_fraction(ai) * to_fraction(bj) for bj in b] for ai in a], dtype=object)


def frac_abs(x):
    return abs(to_fraction(x))


def frac_to_float_matrix(M):
    return np.array([[float(to_fraction(x)) for x in row] for row in M], dtype=float)


def inv2_frac(M):
    """
    Exact inverse of 2x2 matrix with Fraction entries.
    """
    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]

    det = a * d - b * c
    if det == 0:
        raise ValueError("Matrix is singular.")

    inv = frac_array([
        [ d / det, -b / det],
        [-c / det,  a / det],
    ])
    return inv


# ------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------
def primitive_orthogonal(l):
    """
    Compute primitive integer vector orthogonal to l = (l1, l2):
        m = (l2, -l1)
    """
    l1, l2 = l
    return frac_array([l2, -l1])


# ------------------------------------------------------------
# Main result container
# ------------------------------------------------------------
@dataclass
class InverseDesignResult:
    Q: np.ndarray
    Qinv: np.ndarray
    F: np.ndarray
    l3u: np.ndarray


# ------------------------------------------------------------
# Main routine
# ------------------------------------------------------------
def run_inverse_design(A, b1u, b2u, b3u, l1u, l2u, beta):
    """
    Compute deformation gradient F from Burgers vectors and line directions.

    Parameters
    ----------
    A : np.ndarray (2x2)
        Structure matrix of reference lattice.
    b1u, b2u, b3u : arrays (Fraction)
        Burgers vectors in lattice coordinates.
    l1u, l2u : arrays (Fraction)
        Line vectors in lattice coordinates.
    beta : Fraction
        Network factor (1 or 3).

    Returns
    -------
    InverseDesignResult
    """

    # Ensure Fraction arrays
    b1u = frac_array(b1u)
    b2u = frac_array(b2u)
    b3u = frac_array(b3u)
    l1u = frac_array(l1u)
    l2u = frac_array(l2u)
    beta = to_fraction(beta)

    # Third line direction
    l3u = -(l1u + l2u)

    B = [b1u, b2u, b3u]
    L = [l1u, l2u, l3u]

    # Orthogonal vectors
    mbar = [primitive_orthogonal(li) for li in L]

    # z^i values
    z1 = frac_abs(frac_dot(L[1], mbar[0]))
    z2 = frac_abs(frac_dot(L[2], mbar[1]))
    z3 = frac_abs(frac_dot(L[0], mbar[2]))
    Z = [z1, z2, z3]

    if any(z == 0 for z in Z):
        raise ValueError("Degenerate configuration: line vectors are collinear.")

    # Build Q
    Q = frac_array([
        [Fraction(1, 1), Fraction(0, 1)],
        [Fraction(0, 1), Fraction(1, 1)],
    ])

    for i in range(3):
        Q -= frac_outer(B[i], mbar[i]) / (beta * Z[i])

    # Inverse in exact arithmetic
    Qinv = inv2_frac(Q)

    # Convert to float for F
    Qinv_f = frac_to_float_matrix(Qinv)

    # Compute F
    Ainv = np.linalg.inv(A)
    F = A @ Qinv_f @ Ainv

    return InverseDesignResult(
        Q=Q,
        Qinv=Qinv,
        F=F,
        l3u=l3u,
    )


# ------------------------------------------------------------
# Debug / printing utilities (optional)
# ------------------------------------------------------------
def print_fraction_matrix(name, M):
    print(f"\n{name} =")
    for row in M:
        print("  ", [str(x) for x in row])
