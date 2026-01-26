#!/usr/bin/env python3

"""CG amorphous silica script (MDAnalysis)

Three main functions:
  1) build_sio2_molecules(): map each Si to its two nearest O atoms and write a SiO2-molecule GRO
  2) compute_sio2_com(): compute per SiO2 molecule COM and write a COM-only GRO (Si atoms at COMs)
  3) set_all_resids_to_one(): rewrite GRO with all residue ids set to 1

"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances


def build_sio2_molecules(
    amorph_gro: str = "AmorphSilica.gro",
    out_gro: str = "ASiO2.gro",
    n_residues: int | None = None,
) -> None:
    """New Universe where each residue is one SiO2 (Si + 2 nearest O) is created.

    amorph_gro : str
        Input GRO file with atoms named 'Si' and 'O' (amorphous silica).
    out_gro : str
        Output GRO file with SiO2 per residue.
    n_residues : int | None
        Number of Si atoms (residues) to build. If None, inferred from input.
    """
    u = mda.Universe(amorph_gro)

    Si = u.select_atoms("name Si")
    O = u.select_atoms("name O")

    n_Si = len(Si)
    if n_residues is None:
        n_residues = n_Si
    if n_residues != n_Si:
        raise ValueError(f"n_residues={n_residues} does not match number of Si atoms in input ({n_Si}).")

    # Distance matrix: (n_Si, n_O)
    dist_arr = distances.distance_array(Si.positions, O.positions, box=u.dimensions)

    n_atoms = n_residues * 3

    # Create SiO2 Universe: 3 atoms per residue
    resindices = np.repeat(np.arange(n_residues), 3)
    segindices = np.zeros(n_residues, dtype=int)

    sio2 = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=resindices,
        residue_segindex=segindices,
        trajectory=True,
    )

    sio2.add_TopologyAttr("name", ["Si", "O1", "O2"] * n_residues)
    sio2.add_TopologyAttr("type", ["Si", "O", "O"] * n_residues)
    sio2.add_TopologyAttr("resname", ["SIO2"] * n_residues)
    sio2.add_TopologyAttr("resid", list(range(1, n_residues + 1)))
    sio2.add_TopologyAttr("segid", ["SIO2"])

    coords = np.zeros((n_atoms, 3), dtype=float)

    # For each Si: pick two closest O atoms
    for i in range(n_residues):
        os_sorted = np.argsort(dist_arr[i])
        o1 = os_sorted[0]
        o2 = os_sorted[1]

        coords[i * 3 + 0, :] = Si[i].position
        coords[i * 3 + 1, :] = O[o1].position
        coords[i * 3 + 2, :] = O[o2].position

    sio2.atoms.positions = coords
    # Preserve the original box (important for later replication)
    sio2.dimensions = u.dimensions.copy()

    sio2.atoms.write(out_gro)
    print(f"[build_sio2_molecules] Wrote: {out_gro}  (residues={n_residues}, atoms={n_atoms})")


def compute_sio2_com(
    sio2_gro: str = "ASiO2.gro",
    out_com_gro: str = "ASiO2-COM.gro",
    masses: np.ndarray | None = None,
) -> None:
    
    if masses is None:
        masses = np.array([28.0, 16.0, 16.0], dtype=float)

    u = mda.Universe(sio2_gro)
    Si = u.select_atoms("name Si")

    n_Si = len(Si)
    if u.atoms.n_atoms != n_Si * 3:
        raise ValueError(
            f"Expected 3 atoms per SiO2 (total={n_Si*3}), but found {u.atoms.n_atoms} atoms in {sio2_gro}."
        )

    all_centers = np.zeros((n_Si, 3), dtype=float)

    for i in range(n_Si):
        start = i * 3
        stop = start + 3
        trio = u.atoms[start:stop]
        all_centers[i, :] = trio.center(weights=masses)

    # Overwrite Si positions with COMs (keeps single-Si per residue output style from your notebook)
    Si.positions = all_centers

    # Keep box
    Si.dimensions = u.dimensions.copy()

    Si.write(out_com_gro)
    print(f"[compute_sio2_com] Wrote: {out_com_gro}  (Si/COM sites={n_Si})")


def set_all_resids_to_one(
    in_gro: str = "ASiO2-COM.gro",
    out_gro: str = "ASiO2-COM-ResID1.gro",
) -> None:
    
    u = mda.Universe(in_gro)
    for atom in u.atoms:
        atom.residue.resid = 1
    u.atoms.write(out_gro)
    print(f"[set_all_resids_to_one] Wrote: {out_gro}")



def main() -> None:
    
    build_sio2_molecules("AmorphSilica.gro", "ASiO2.gro")
    compute_sio2_com("ASiO2.gro", "ASiO2-COM.gro")
    set_all_resids_to_one("ASiO2-COM.gro", "ASiO2-COM-ResID1.gro")


if __name__ == "__main__":
    main()
