import numpy as np
from fractions import Fraction
import pyoilab as gb

from bilayers import (
    make_ab_graphene_bilayer,
    make_aa_prime_mos2_bilayer,
    heterodeform_bilayer,
)
from inverse_design_core import run_inverse_design, print_fraction_matrix
from lammps_writer import write_lammps_data, summarize_configuration


# ------------------------------------------------------------
# Step 1: Choose reference bilayer
# ------------------------------------------------------------
bilayer = make_ab_graphene_bilayer(6.0)
# bilayer = make_aa_prime_mos2_bilayer(6.5)


# ------------------------------------------------------------
# Step 2: Specify inverse-design input
# ------------------------------------------------------------
# Example shown here: 1D simple network
# Replace these using the combinations listed in README.md.
b1u = np.array([Fraction(-1, 1), Fraction(1, 1)], dtype=object)
b2u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)
b3u = np.array([Fraction(0, 1), Fraction(0, 1)], dtype=object)

l1u = np.array([Fraction(-1, 1), Fraction(1, 1)], dtype=object)
l2u = np.array([Fraction(569, 1), Fraction(567, 1)], dtype=object)

beta = Fraction(1, 1)

# ------------------------------------------------------------
# Step 3: Compute deformation gradient F
# ------------------------------------------------------------
# By default we use the bottom-layer lattice as the reference lattice in the
# inverse-design formulas.
result = run_inverse_design(
    A=bilayer.bottom_layer.A,
    b1u=b1u,
    b2u=b2u,
    b3u=b3u,
    l1u=l1u,
    l2u=l2u,
    beta=beta,
)

F = result.F

print_fraction_matrix("Q", result.Q)

print_fraction_matrix("Q^{-1}", result.Qinv)

print("\nF = ")
print(F)


# ------------------------------------------------------------
# Step 4: Heterodeform the bilayer
# ------------------------------------------------------------
# Default convention in this example:
#   bottom layer is deformed, top layer is fixed
# The user may freely change this.
hetero_bilayer = heterodeform_bilayer(
    bilayer,
    F_bottom=F,
    F_top=np.eye(2),
)

sigmaA = abs(hetero_bilayer.bicrystal.sigmaA)
sigmaB = abs(hetero_bilayer.bicrystal.sigmaB)
exp_number_of_atoms_top_layer = sigmaA * len(hetero_bilayer.top_layer.basis_atoms)
exp_number_of_atoms_bottom_layer = sigmaB * len(hetero_bilayer.bottom_layer.basis_atoms)
print("\nExpected number of atoms in the top layer = ", exp_number_of_atoms_top_layer)
print("Expected number of atoms in the bottom layer = ", exp_number_of_atoms_bottom_layer)
print("Expected total number of atoms in the bilayer = ", exp_number_of_atoms_top_layer+exp_number_of_atoms_bottom_layer)

# ------------------------------------------------------------
# Step 5: Choose CSL vectors
# ------------------------------------------------------------
bicrystal = hetero_bilayer.bicrystal

# Primitive CSL vectors: The user may instead choose any other integer combination of CSL basis vectors.
ell1 = gb.LatticeVector2D(np.array([1, 0], dtype=np.int64), bicrystal.csl)
ell2 = gb.LatticeVector2D(np.array([0, 1], dtype=np.int64), bicrystal.csl)



# ------------------------------------------------------------
# Step 6: Build atomic positions in the chosen CSL cell
# ------------------------------------------------------------
positions, atom_types, charges, labels = hetero_bilayer.box(ell1, ell2)

print(f"\nNumber of atoms in box = {len(positions)}")


# ------------------------------------------------------------
# Step 7: Optional rotation so ell1 aligns with the x-axis
# ------------------------------------------------------------
do_rotate = True

ell1_cart = np.array(ell1.cartesian(), dtype=float)
ell2_cart = np.array(ell2.cartesian(), dtype=float)

if do_rotate:
    theta = np.arctan2(ell1_cart[1], ell1_cart[0])

    R = np.array([
        [np.cos(-theta), -np.sin(-theta)],
        [np.sin(-theta),  np.cos(-theta)],
    ])

    positions = [
        tuple(R @ np.array([x, y])) + (z,)
        for (x, y, z) in positions
    ]

    ell1_cart = R @ ell1_cart
    ell2_cart = R @ ell2_cart


# ------------------------------------------------------------
# Step 8: Write LAMMPS data file
# ------------------------------------------------------------
summarize_configuration(positions, atom_types, ell1_cart, ell2_cart)
write_lammps_data(
    filename="initial.data",
    positions=positions,
    atom_types=atom_types,
    charges=charges,
    ell1=ell1_cart,
    ell2=ell2_cart,
    atom_style=bilayer.atom_style,
)

print("\nWrote LAMMPS data file: initial.data")
