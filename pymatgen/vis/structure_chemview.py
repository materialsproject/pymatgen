"""Visualization for structures using chemview."""

from __future__ import annotations

import numpy as np
from monty.dev import requires

from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    from chemview import MolecularViewer
    from chemview.utils import get_atom_color

    chemview_loaded = True
except ImportError:
    chemview_loaded = False


@requires(chemview_loaded, "To use quick_view, you need to have chemview installed.")
def quick_view(
    structure,
    bonds=True,
    conventional=False,
    transform=None,
    show_box=True,
    bond_tol=0.2,
    stick_radius=0.1,
):
    """
    A function to visualize pymatgen Structure objects in jupyter notebook using chemview package.

    Args:
        structure: pymatgen Structure
        bonds: (bool) visualize bonds. Bonds are found by comparing distances
            to added covalent radii of pairs. Defaults to True.
        conventional: (bool) use conventional cell. Defaults to False.
        transform: (list) can be used to make supercells with pymatgen.Structure.make_supercell method
        show_box: (bool) unit cell is shown. Defaults to True.
        bond_tol: (float) used if bonds=True. Sets the extra distance tolerance when finding bonds.
        stick_radius: (float) radius of bonds.

    Returns:
        A chemview.MolecularViewer object
    """
    struct = structure.copy()
    if conventional:
        struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()

    if transform:
        struct.make_supercell(transform)
    atom_types = [specie.symbol for specie in struct.species]

    if bonds:
        bonds = []
        for idx in range(len(struct) - 1):
            sym_i = struct[idx].specie.symbol
            for j in range(idx + 1, len(struct)):
                sym_j = struct[j].specie.symbol
                max_d = CovalentRadius.radius[sym_i] + CovalentRadius.radius[sym_j] + bond_tol
                if struct.get_distance(idx, j, np.array([0, 0, 0])) < max_d:
                    bonds.append((idx, j))
    bonds = bonds or None

    mv = MolecularViewer(struct.cart_coords, topology={"atom_types": atom_types, "bonds": bonds})

    if bonds:
        mv.ball_and_sticks(stick_radius=stick_radius)
    for site in struct:
        el = site.specie.symbol
        coord = site.coords
        rad = CovalentRadius.radius[el]
        mv.add_representation(
            "spheres",
            {"coordinates": coord.astype(float), "colors": [get_atom_color(el)], "radii": [rad * 0.5], "opacity": 1},
        )
    if show_box:
        o = np.array([0, 0, 0])
        a, b, c = struct.lattice.matrix[0], struct.lattice.matrix[1], struct.lattice.matrix[2]
        starts = [o, o, o, a, a, b, b, c, c, a + b, a + c, b + c]
        ends = [
            a,
            b,
            c,
            a + b,
            a + c,
            b + a,
            b + c,
            c + a,
            c + b,
            a + b + c,
            a + b + c,
            a + b + c,
        ]
        colors = [0xFFFFFF for i in range(12)]
        mv.add_representation(
            "lines",
            {
                "startCoords": np.array(starts),
                "endCoords": np.array(ends),
                "startColors": colors,
                "endColors": colors,
            },
        )
    return mv
