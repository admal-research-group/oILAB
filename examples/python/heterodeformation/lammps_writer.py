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
    z_padding=20.0,
    comment="LAMMPS data file",
):
    """
    Write a LAMMPS data file.

    Parameters
    ----------
    filename : str
        Output file name.
    positions : list[(x, y, z)]
        Cartesian atomic positions.
    atom_types : list[int]
        Atom types.
    charges : list[float]
        Atomic charges. Used only for atom_style='full'.
    ell1, ell2 : array-like
        2D Cartesian simulation cell vectors.
    atom_style : str
        'atomic' or 'full'.
    z_padding : float
        Half-width of z box: [ -z_padding, z_padding ].
    comment : str
        Header comment.
    """
    if len(positions) != len(atom_types):
        raise ValueError("positions and atom_types must have the same length.")

    if atom_style == "full" and len(positions) != len(charges):
        raise ValueError("positions and charges must have the same length for atom_style='full'.")

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
            for i, ((x, y, z), atom_type, charge) in enumerate(
                zip(positions, atom_types, charges), start=1
            ):
                mol_id = 1
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
