from dataclasses import dataclass, field
import numpy as np
import pyoilab as gb

from layers import (
    Layer2D,
    make_graphene_layer_from_basis,
    make_mos2_layer_from_basis,
)


@dataclass(frozen=True)
class Bilayer2D:
    """
    Bilayer consisting of two 2D multilattices.

    Parameters
    ----------
    name : str
        Name of the bilayer.
    top_layer : Layer2D
        Top layer.
    bottom_layer : Layer2D
        Bottom layer.
    interlayer_spacing : float
        Distance between the Bravais lattice planes of the two layers.
    """
    name: str
    top_layer: Layer2D
    bottom_layer: Layer2D
    interlayer_spacing: float

    # Constructed once when the bilayer is created
    bicrystal: gb.BiCrystal2D = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(
            self,
            "bicrystal",
            gb.BiCrystal2D(
                self.bottom_layer.lattice,
                self.top_layer.lattice,
                True,
            ),
        )

    @property
    def atom_style(self) -> str:
        if self.top_layer.atom_style != self.bottom_layer.atom_style:
            raise ValueError("Top and bottom layers must use the same atom_style.")
        return self.top_layer.atom_style

    def box(self, ell1, ell2):
        """
        Build Cartesian atomic positions using CSL lattice vectors.

        Parameters
        ----------
        ell1, ell2 : gb.LatticeVector2D
            CSL lattice vectors belonging to self.bicrystal.csl.
            These may be primitive vectors or any integer combination of them.

        Returns
        -------
        positions : list[(x, y, z)]
        atom_types : list[int]
        charges : list[float]
        labels : list[str]
        molecule_ids : list[int]
        """
        # Convert CSL vectors into each layer's lattice vectors
        ell1_bottom = self.bicrystal.getLatticeVectorInA(ell1)
        ell2_bottom = self.bicrystal.getLatticeVectorInA(ell2)

        ell1_top = self.bicrystal.getLatticeVectorInB(ell1)
        ell2_top = self.bicrystal.getLatticeVectorInB(ell2)

        # Build each layer independently
        pos_bottom, types_bottom, charges_bottom, labels_bottom = \
            self.bottom_layer.box([ell1_bottom, ell2_bottom])

        pos_top, types_top, charges_top, labels_top = \
            self.top_layer.box([ell1_top, ell2_top])

        # Shift top layer by interlayer spacing
        pos_top = [
            (x, y, z + self.interlayer_spacing)
            for (x, y, z) in pos_top
        ]

        # Offset atom types in top layer
        type_offset = max(types_bottom)
        types_top = [t + type_offset for t in types_top]

        # Molecule IDs: bottom = 1, top = 2
        molecule_ids_bottom = [1] * len(pos_bottom)
        molecule_ids_top = [2] * len(pos_top)

        # Concatenate
        positions = pos_bottom + pos_top
        atom_types = types_bottom + types_top
        charges = charges_bottom + charges_top
        labels = labels_bottom + labels_top
        molecule_ids = molecule_ids_bottom + molecule_ids_top

        return positions, atom_types, charges, labels, molecule_ids


def heterodeform_bilayer(
    bilayer: Bilayer2D,
    F_bottom=None,
    F_top=None,
    new_name=None,
):
    """
    Apply in-plane deformation gradients to the two layers of a bilayer.

    Default behavior: no deformation unless provided explicitly.

    Parameters
    ----------
    bilayer : Bilayer2D
        Input bilayer.
    F_bottom : np.ndarray or None
        Deformation gradient for the bottom layer.
    F_top : np.ndarray or None
        Deformation gradient for the top layer.
    new_name : str or None
        Optional name for the deformed bilayer.

    Returns
    -------
    Bilayer2D
        New bilayer with deformed layers. Its bicrystal is constructed once
        during Bilayer2D initialization.
    """
    if F_bottom is None:
        F_bottom = np.eye(2)
    if F_top is None:
        F_top = np.eye(2)

    bottom_layer = bilayer.bottom_layer.deform(F_bottom)
    top_layer = bilayer.top_layer.deform(F_top)

    return Bilayer2D(
        name=bilayer.name if new_name is None else new_name,
        top_layer=top_layer,
        bottom_layer=bottom_layer,
        interlayer_spacing=bilayer.interlayer_spacing,
    )

def make_ab_graphene_bilayer(interlayer_spacing) -> Bilayer2D:
    """
    AB-stacked bilayer graphene reference configuration.

    Top layer basis:
        s1 = (0, 0)
        s2 = (1/3, 2/3)

    Bottom layer basis:
        t1 = (0, 0)
        t2 = (2/3, 1/3)
    """
    top_layer = make_graphene_layer_from_basis(
        name="graphene_top",
        basis_frac_xy=[
            (0.0, 0.0),
            (1.0 / 3.0, 2.0 / 3.0),
        ],
    )

    bottom_layer = make_graphene_layer_from_basis(
        name="graphene_bottom",
        basis_frac_xy=[
            (0.0, 0.0),
            (2.0 / 3.0, 1.0 / 3.0),
        ],
    )

    return Bilayer2D(
        name="graphene_ab",
        top_layer=top_layer,
        bottom_layer=bottom_layer,
        interlayer_spacing=interlayer_spacing,
    )


def make_aa_prime_mos2_bilayer(interlayer_spacing) -> Bilayer2D:
    """
    AA'-stacked bilayer MoS2 reference configuration.

    Top layer basis:
        s1 = (0, 0, 0)
        s2 = (0, 0, 1)
        s3 = (1/3, 2/3, 1/2)

    Bottom layer basis:
        t1 = (1/3, 2/3, 0)
        t2 = (1/3, 2/3, 1)
        t3 = (0, 0, 1/2)
    """
    t = 3.1902

    top_layer = make_mos2_layer_from_basis(
        name="mos2_top",
        basis_data=[
            ("S",  (0.0, 0.0),             0.0 * t, 1, -0.42),
            ("S",  (0.0, 0.0),             1.0 * t, 2, -0.42),
            ("Mo", (1.0 / 3.0, 2.0 / 3.0), 0.5 * t, 3,  0.84),
        ],
    )

    bottom_layer = make_mos2_layer_from_basis(
        name="mos2_bottom",
        basis_data=[
            ("S",  (1.0 / 3.0, 2.0 / 3.0), 0.0 * t, 1, -0.42),
            ("S",  (1.0 / 3.0, 2.0 / 3.0), 1.0 * t, 2, -0.42),
            ("Mo", (0.0, 0.0),             0.5 * t, 3,  0.84),
        ],
    )

    return Bilayer2D(
        name="mos2_aa_prime",
        top_layer=top_layer,
        bottom_layer=bottom_layer,
        interlayer_spacing=interlayer_spacing,
    )
