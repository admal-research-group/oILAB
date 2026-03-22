from dataclasses import dataclass
import numpy as np
import pyoilab as gb


@dataclass(frozen=True)
class BasisAtom:
    label: str
    frac_xy: np.ndarray
    z_rel: float
    atom_type: int
    charge: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "frac_xy", np.asarray(self.frac_xy, dtype=float).reshape(2))


@dataclass(frozen=True)
class Layer2D:
    name: str
    lattice: gb.Lattice2D
    basis_atoms: list[BasisAtom]
    atom_style: str = "atomic"

    @property
    def A(self):
        return np.array(self.lattice.latticeBasis, dtype=float)

    def deform(self, F):
        """
        Apply in-plane deformation to the layer.
        """
        return Layer2D(
            name=f"{self.name}_def",
            lattice=gb.Lattice2D(self.A, F),
            basis_atoms=self.basis_atoms,
            atom_style=self.atom_style,
        )

    def box(self, lattice_vectors):
        """
        Build Cartesian atomic positions for this layer.
        """
        box_obj = self.lattice.box(lattice_vectors)

        positions = []
        atom_types = []
        charges = []
        labels = []

        for site in box_obj:
            site_xy = np.array(site.cartesian(), dtype=float)

            for atom in self.basis_atoms:
                basis_xy = self.A @ atom.frac_xy
                xy = site_xy + basis_xy

                positions.append((xy[0], xy[1], atom.z_rel))
                atom_types.append(atom.atom_type)
                charges.append(atom.charge)
                labels.append(atom.label)

        return positions, atom_types, charges, labels

def triangular_lattice_structure_matrix(a: float) -> np.ndarray:
    """
    Structure matrix for the 2D triangular Bravais lattice used in the manuscript:
        A = (a/2) [[0, -sqrt(3)],
                   [2, -1]]
    """
    return (a / 2.0) * np.array(
        [[0.0, -np.sqrt(3.0)],
         [2.0, -1.0]],
        dtype=float,
    )


def make_graphene_layer_from_basis(name: str, basis_frac_xy: list[tuple[float, float]]) -> Layer2D:
    """
    Construct a graphene layer from a specified two-atom in-plane basis.
    """
    a = 2.46
    A = triangular_lattice_structure_matrix(a)
    lattice = gb.Lattice2D(A)

    basis_atoms = [
        BasisAtom(label="C", frac_xy=np.array(basis_frac_xy[0]), z_rel=0.0, atom_type=1, charge=0.0),
        BasisAtom(label="C", frac_xy=np.array(basis_frac_xy[1]), z_rel=0.0, atom_type=1, charge=0.0),
    ]

    return Layer2D(
        name=name,
        lattice=lattice,
        basis_atoms=basis_atoms,
        atom_style="atomic",
    )


def make_mos2_layer_from_basis(
    name: str,
    basis_data: list[tuple[str, tuple[float, float], float, int, float]],
) -> Layer2D:
    """
    Construct a MoS2 layer from a specified basis.

    Each basis entry is:
        (label, frac_xy, z_rel, atom_type, charge)
    """
    a = 3.19702
    A = triangular_lattice_structure_matrix(a)
    lattice = gb.Lattice2D(A)

    basis_atoms = [
        BasisAtom(
            label=label,
            frac_xy=np.array(frac_xy, dtype=float),
            z_rel=float(z_rel),
            atom_type=int(atom_type),
            charge=float(charge),
        )
        for (label, frac_xy, z_rel, atom_type, charge) in basis_data
    ]

    return Layer2D(
        name=name,
        lattice=lattice,
        basis_atoms=basis_atoms,
        atom_style="full",
    )
