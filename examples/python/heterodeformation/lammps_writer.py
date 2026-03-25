import numpy as np


def compute_triclinic_box(ell1, ell2, zlo=-20.0, zhi=20.0):
    """
    Convert two 2D Cartesian box vectors into LAMMPS triclinic box parameters.

    Assumes the user has already rotated the system, if desired.

    Parameters
    ----------
    ell1, ell2 : array-like
        2D Cartesian simulation cell vectors.
    zlo, zhi : float
        Lower and upper z bounds.

    Returns
    -------
    dict
        Dictionary with xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz.
    """
    ell1 = np.asarray(ell1, dtype=float).reshape(2)
    ell2 = np.asarray(ell2, dtype=float).reshape(2)

    # Standard 2D triclinic reduction used by LAMMPS
    lx = np.linalg.norm(ell1)
    if lx <= 0.0:
        raise ValueError("ell1 has zero length.")

    ex = ell1 / lx
    xy = np.dot(ell2, ex)

    ell2_perp = ell2 - xy * ex
    ly = np.linalg.norm(ell2_perp)
    if ly <= 0.0:
        raise ValueError("ell1 and ell2 are collinear.")

    return {
        "xlo": 0.0,
        "xhi": float(lx),
        "ylo": 0.0,
        "yhi": float(ly),
        "zlo": float(zlo),
        "zhi": float(zhi),
        "xy": float(xy),
        "xz": 0.0,
        "yz": 0.0,
    }


def write_lammps_data(
    filename,
    positions,
    atom_types,
    charges,
    ell1,
    ell2,
    atom_style="atomic",
    molecule_ids=None,
    z_padding=20.0,
    comment="LAMMPS data file",
):
    """
    Write a LAMMPS data file for atomistic simulations.

    This function generates a triclinic simulation box from the in-plane
    lattice vectors (ell1, ell2) and writes atomic data in either
    'atomic' or 'full' atom_style format.

    Parameters
    ----------
    filename : str
        Output file name.

    positions : list[(x, y, z)]
        Cartesian atomic positions.

    atom_types : list[int]
        Integer atom types (1-based), consistent with LAMMPS conventions.

    charges : list[float]
        Atomic charges. Used only when atom_style='full'. Ignored otherwise.

    ell1, ell2 : array-like
        2D Cartesian lattice vectors defining the in-plane periodic cell.
        These are converted internally into a triclinic LAMMPS box.

    atom_style : str, optional
        LAMMPS atom style. Supported options:
        - 'atomic' : writes (id, type, x, y, z)
        - 'full'   : writes (id, molecule-ID, type, charge, x, y, z)

    molecule_ids : list[int] or None, optional
        Molecule (or layer) IDs for each atom, required for certain
        many-body or interlayer potentials (e.g., ilp/tmd).
        If None and atom_style='full', all atoms are assigned molecule ID = 1.

    z_padding : float, optional
        Half-width of the simulation box in the z-direction. The box is
        constructed as [ -z_padding, z_padding ]. This should be large
        enough to avoid spurious interactions between periodic images
        in z when using boundary p p p.

    comment : str, optional
        Header comment written at the top of the data file.
    """
    if len(positions) != len(atom_types):
        raise ValueError("positions and atom_types must have the same length.")

    if atom_style == "full" and len(positions) != len(charges):
        raise ValueError("positions and charges must have the same length for atom_style='full'.")

    if atom_style == "full" and molecule_ids is not None and len(positions) != len(molecule_ids):
        raise ValueError("positions and molecule_ids must have the same length for atom_style='full'.")

    box = compute_triclinic_box(ell1, ell2, zlo=-z_padding, zhi=z_padding)

    n_atoms = len(positions)
    n_atom_types = max(atom_types)

    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"{comment}\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_atom_types} atom types\n\n")

        f.write(f"{box['xlo']:.16f} {box['xhi']:.16f} xlo xhi\n")
        f.write(f"{box['ylo']:.16f} {box['yhi']:.16f} ylo yhi\n")
        f.write(f"{box['zlo']:.16f} {box['zhi']:.16f} zlo zhi\n")
        f.write(f"{box['xy']:.16f} {box['xz']:.16f} {box['yz']:.16f} xy xz yz\n\n")

        if atom_style == "atomic":
            f.write("Atoms # atomic\n\n")
            for i, ((x, y, z), atom_type) in enumerate(zip(positions, atom_types), start=1):
                f.write(f"{i:d} {atom_type:d} {x:.16f} {y:.16f} {z:.16f}\n")

        elif atom_style == "full":
            f.write("Atoms # full\n\n")
            if molecule_ids is None:
                molecule_ids = [1] * n_atoms

            for i, ((x, y, z), mol_id, atom_type, charge) in enumerate(
                zip(positions, molecule_ids, atom_types, charges), start=1
            ):
                f.write(
                    f"{i:d} {mol_id:d} {atom_type:d} "
                    f"{charge:.16f} {x:.16f} {y:.16f} {z:.16f}\n"
                )
        else:
            raise ValueError("atom_style must be either 'atomic' or 'full'.")


def summarize_configuration(positions, atom_types, ell1, ell2):
    """
    Print a compact summary of the configuration.
    """
    print("\nConfiguration summary")
    print("---------------------")
    print(f"Number of atoms      : {len(positions)}")
    print(f"Number of atom types : {max(atom_types)}")
    print(f"ell1                 : {np.asarray(ell1, dtype=float)}")
    print(f"ell2                 : {np.asarray(ell2, dtype=float)}")
